#!/bin/bash
# F4 — SPECIFICITY (the canonical-hash regression guard).
#
# subjectA and subjectB are simulated from the SAME genome (fixtures/ref/parent.fa) with
# different wgsim seeds. There are genuinely ZERO real differences between them, so RUFUS
# must call NOTHING. Anything it emits is a false positive by construction.
#
# This is the test that catches the class of bug where non-canonical control hashes leak
# germline k-mers into the call set — a regression that recall-only tests cannot see.
#
# Run directly:  bash f4_specificity.sh
# Or submit:     sbatch f4_specificity.sh
#SBATCH --account=marth-rw
#SBATCH --partition=marth-rw
#SBATCH --job-name=rufus_f4_specificity
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
set -euo pipefail

# sbatch COPIES this script into the slurmd spool dir, so BASH_SOURCE points at the copy and
# any path derived from it resolves against the spool, not the repo. Prefer SLURM_SUBMIT_DIR
# (set by sbatch to the dir you submitted from) and fall back to script location for direct runs.
HERE="${FUNCTIONAL_CASES_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}}"
FIX="$(cd "$HERE/../fixtures" 2>/dev/null && pwd)" || {
  echo "ERROR: cannot locate fixtures from HERE=$HERE"
  echo "       submit from tests/functional/cases/ or set FUNCTIONAL_CASES_DIR=<abs path to cases/>"
  exit 1; }
DATA=/uufs/chpc.utah.edu/common/HIPAA/u0746015/marth_software/RUFUS/resources/reg_test_files
SIF=${SIF:-/uufs/chpc.utah.edu/common/HIPAA/u0746015/marth_software/RUFUS/zenodo_images/rufus_dev.sif}
OUT=${OUT:-$DATA/runs/f4_specificity}
THREADS=${SLURM_CPUS_PER_TASK:-8}

REF=$FIX/ref/tiny.fa
SUBJ=$FIX/specificity/subjectA.bam
CTRL=$FIX/specificity/subjectB.bam

rm -rf "$OUT"; mkdir -p "$OUT"; cd "$OUT"

module load apptainer 2>/dev/null || true
command -v apptainer >/dev/null || { echo "ERROR: apptainer unavailable"; exit 1; }
[ -f "$SIF" ] || { echo "ERROR: SIF not found: $SIF"; exit 1; }
for f in "$REF" "$SUBJ" "$CTRL"; do
  [ -f "$f" ] || { echo "ERROR: fixture missing: $f (run make_fixtures.sh)"; exit 1; }
done

echo "=== F4 specificity | $(date) ==="
echo "  subject : $(basename "$SUBJ")"
echo "  control : $(basename "$CTRL")   (SAME genome, different seed)"
echo "  expect  : ZERO variants called"
echo

# NOTE: -R chr20 keeps RUFUS in windowed mode (1G hash default, ~4GB). Omitting -R would select
# whole-genome mode and its much larger hash default — needless for a 200kb fixture. No -hs:
# the shipped default is what we want to test. See run 17378733 for why overriding it is fatal.
set +e
apptainer exec --bind "$FIX,$OUT,$DATA" "$SIF" \
  bash /opt/RUFUS/runRufus.sh \
    -s "$SUBJ" -c "$CTRL" -r "$REF" \
    -k 25 -m 5 -R chr20 -t "$THREADS" -z
RC=$?
set -e
echo "runRufus exit: $RC"

# ---- assertions ----
STATUS=$(cat "$OUT/region_status.log" 2>/dev/null || echo "<no region_status.log>")
echo "region_status: $STATUS"
WORK="$OUT/rufus_chr20"
fail() { echo "RESULT: FAIL — $*"; exit 1; }

# PRECONDITION: prove the pipeline actually RAN before trusting an empty result.
# A negative test is vacuous if it also passes when nothing executed. Run 17378733 produced
# NO_VARIANTS too — reason no_control_kmers — because jellyfish had been OOM-killed. So assert
# real k-mer counting happened for BOTH samples first.
echo "-- preconditions --"
for s in subjectA subjectB; do
  JH="$WORK/$s.bam.chr20.generator.Jhash"; HI="$JH.histo"
  [ -s "$JH" ] || fail "$s: Jhash missing/empty — jellyfish never produced a hash (OOM signature)"
  [ -s "$HI" ] || fail "$s: histo missing/empty"
  TOT=$(awk '{s+=$2} END{print s+0}' "$HI")
  [ "$TOT" -ge 100000 ] || fail "$s: only $TOT distinct kmers (<100000) — counting did not really run"
  # NB: `sort -k2,2nr | head -1` would SIGPIPE sort, and with `set -euo pipefail` that kills the
  # script silently before fail() can report. Single-process awk avoids the pipe entirely.
  MODE=$(awk 'BEGIN{m=-1;c=0} {if($2+0>m){m=$2+0;c=$1+0}} END{print c}' "$HI")
  { [ "$MODE" -ge 12 ] && [ "$MODE" -le 45 ]; } || fail "$s: modal kmer depth $MODE outside 12-45 (expect ~24 at 30x)"
  echo "   $s: distinct=$TOT modal_depth=$MODE  OK"
done

# The reason must be no_unique_hashes specifically. no_control_kmers / ERROR mean a stage DIED.
grep -q 'no_unique_hashes' <<<"$STATUS" \
  || fail "status reason is not 'no_unique_hashes' (a dead stage can also yield NO_VARIANTS): $STATUS"

# The specificity result at k-mer level: nothing survived subject-minus-control.
HL=$(ls "$WORK"/subjectA*.HashList 2>/dev/null | head -1 || true)
[ -n "$HL" ] && [ -f "$HL" ] || fail "HashList missing — subtraction stage did not run"
NHL=$(wc -l < "$HL")
[ "$NHL" -eq 0 ] || fail "$NHL unique subject kmers survived subtraction on identical genomes (expected 0)"
echo "   HashList: 0 unique kmers  OK"

# And no VCF records, if a VCF was emitted at all.
VCF=$(ls "$OUT"/temp.RUFUS.Final.*.vcf.gz "$OUT"/RUFUS.Final.*.vcf.gz 2>/dev/null | head -1 || true)
if [ -n "$VCF" ]; then
  N=$(zcat "$VCF" | grep -vc '^#' || true)
  echo "   VCF: $VCF records=$N"
  [ "$N" -eq 0 ] || { echo "RESULT: FAIL — $N false positive(s) on identical genomes:";
                      zcat "$VCF" | grep -v '^#' | head -20; exit 1; }
fi

echo "RESULT: PASS — pipeline ran on both samples, 0 unique kmers, 0 variants called"
exit 0
