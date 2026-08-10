#!/bin/bash
set -e
# pipefail so the feeder pipeline below reports a failure in `bash "$GEN"` and not just in the
# samtools stage that terminates it.
set -o pipefail
GEN=$1
K=$2
T=$3
L=$4
HASH_SIZE=$5

: "${RUFUS_ROOT:=/opt/RUFUS}"
RDIR="$RUFUS_ROOT"
JELLYFISH="$RDIR/bin/externals/jellyfish/src/jellyfish_project/bin/jellyfish"

trap 'rm -f "$FIFO_FQ"' EXIT

# If we're using a region-specific hash, adjust size accordingly (1MB hashes made w/ 1G)
if [ -e "$GEN.Jhash" ]
then
	echo "Skipping jelly, $GEN.Jhash alreads exists"
else
	# Extra thread safety
	PID=$$
	FIFO_FQ="${GEN}.fq.${PID}"

	rm -f "$FIFO_FQ"
	mkfifo "$FIFO_FQ"

	# bash "$GEN" | "$RDIR/bin/PassThroughSamCheck" "$GEN.Jelly.chr" > "$FIFO_FQ" &
	# samtools fastq validated bit-identical to PassThroughSamCheck for counting (job 16682513); the
	# generator's first command emits -h so the header is present -- and only its first, since
	# samtools aborts on a second @HD mid-stream.
	bash "$GEN" | samtools fastq -@ "$T" - > "$FIFO_FQ" &
	FEEDER=$!

	# -C is canonical ("Count both strand, canonical representation")
	# -L is filtering out low frequency kmers ("Don't output k-mer with count < lower-count")
	# For subject, we keep kmers with 2+ counts
	# For controls, we keep kmers with 2+ counts OR the provided argument to rufus (_argParLowK)
	# These arguments combined, I interpret this as keeping kmers with a single read
	# --disk means hash will be written to disk if entire thing can't be held in memory
	# -s (intial hash size) is G + Gcek (genome size * coverage * error * kmer length) ~228G for 300x, 22.8G for 30x, etc
	# guessing this starting number is far too low and there's a lot of memory swapping happening here
	# good area of parallelization and possible merging after - will neeed to think through
	
	# Capture jellyfish's status explicitly rather than letting `set -e` abort with it.
	# Callers must be able to tell "the tool failed" (OOM, disk full, crash) apart from
	# "this region legitimately has no k-mers" -- both otherwise leave an empty histogram
	# and were previously indistinguishable. In a sharded whole-genome run that turns a
	# lost shard into a silent "no variants here". See exit-code contract below.
	set +e
	"$JELLYFISH" count --disk -m "$K" -L "$L" -s "$HASH_SIZE" -t "$T" -o "$GEN.Jhash" -C "$FIFO_FQ"
	jf_rc=$?

	if [ "$jf_rc" -ne 0 ]; then
		# jellyfish is gone, so nothing is draining the FIFO and the feeder is blocked
		# mid-write. Tear it down before reaping, or `wait` never returns.
		rm -f "$FIFO_FQ"
		kill "$FEEDER" 2>/dev/null
		wait "$FEEDER" 2>/dev/null
		set -e
		echo "ERROR: jellyfish count failed (exit $jf_rc) for $GEN" >&2
		exit 2
	fi

	# The feeder's status is not jellyfish's. If the feeder dies partway -- a truncated stream,
	# an unreadable input, a malformed generator -- jellyfish sees a clean EOF on the FIFO and
	# reports success over however many reads happened to arrive, so an under-counted hash
	# looks identical to a complete one. Reap it explicitly and treat a failure as a tool
	# failure (2), discarding the partial hash so a rerun cannot pick it up via the
	# skip-if-exists check at the top.
	wait "$FEEDER"
	feeder_rc=$?
	set -e

	if [ "$feeder_rc" -ne 0 ]; then
		rm -f "$GEN.Jhash"
		echo "ERROR: read feeder failed (exit $feeder_rc) for $GEN; k-mer counts would be incomplete" >&2
		exit 2
	fi

	# A zero exit with no output file means jellyfish died without reporting it.
	if [ ! -s "$GEN.Jhash" ]; then
		echo "ERROR: jellyfish count reported success but produced no $GEN.Jhash" >&2
		exit 2
	fi
fi

if [ ! -s "$GEN.Jhash.histo" ]; then
	set +e
	"$JELLYFISH" histo -f -o "$GEN.Jhash.histo" "$GEN.Jhash"
	histo_rc=$?
	set -e
	if [ "$histo_rc" -ne 0 ] || [ ! -s "$GEN.Jhash.histo" ]; then
		echo "ERROR: jellyfish histo failed (exit $histo_rc) for $GEN" >&2
		exit 2
	fi
fi

# Exit-code contract for callers (runRufus.sh check_empty_hashes depends on this):
#   0 - counted OK, k-mers found
#   1 - ran successfully, but the region genuinely contains no k-mers
#   2 - the counting tool itself failed; the result says nothing about coverage
if [ $(awk '$2 > 0' "$GEN.Jhash.histo" | wc -l ) -eq "0" ]; then
	exit 1
fi
exit 0