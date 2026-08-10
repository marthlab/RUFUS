#!/bin/bash
# F7 — MULTI-SUBJECT INPUT (-s repeated).
#
# A sample sequenced across several files (different centers, flowcells, or split lanes) is passed
# as repeated -s arguments; setup_slurm.sh emits exactly that form from its comma-separated -s list.
# runRufus.sh concatenates one read-emitting command per subject into a single generator, and the
# generator body is run as ONE stream (`bash generator | samtools ...`). So every subject file must
# contribute its reads while exactly ONE command emits a SAM header -- samtools aborts on a second
# @HD mid-stream, and because `collate` buffers before emitting, it discards the WHOLE stream rather
# than truncating it. The visible symptom is a clean-looking run that calls nothing at all, in every
# region ("no_reads_passed_filter"), with the real error buried in stderr.
#
# The test splits the somatic tumor fixture into two files by read name -- a true partition, so the
# union is exactly the original single file -- and asserts a multi-subject run reproduces the
# single-file result:
#   PARTITION    the two parts are non-empty and their read counts sum to the original
#   HEADER       the multi-subject generator has 2 read commands but emits exactly 1 @HD
#   RECALL       every planted somatic is recovered from the split input
#   CONCORDANCE  split call set == single-file call set (bam), and the cram arm matches too
#
# The HEADER assertion is what pinpoints a regression: without it, a reintroduced duplicate header
# shows up only as "no calls", which is indistinguishable from a dozen unrelated failures. F3 already
# establishes bam/cram call-set equality, so the cram arm here isolates the -T/-cr generator line
# rather than re-testing format concordance.
#
# Run directly:  bash f7_multi_subject.sh
# Or submit:     sbatch f7_multi_subject.sh
# EXTRA_BIND=host:container[,...] overlays uncommitted code (dev loop) before the image is rebuilt.
#SBATCH --account=marth-rw
#SBATCH --partition=marth-rw
#SBATCH --job-name=rufus_f7_multi_subject
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
set -euo pipefail

# sbatch copies this script to the spool dir, so BASH_SOURCE-derived paths break. See F4.
HERE="${FUNCTIONAL_CASES_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}}"
FIX="$(cd "$HERE/../fixtures" 2>/dev/null && pwd)" || { echo "ERROR: cannot locate fixtures from HERE=$HERE"; exit 1; }
DATA=/uufs/chpc.utah.edu/common/HIPAA/u0746015/marth_software/RUFUS/resources/reg_test_files
SIF=${SIF:-/uufs/chpc.utah.edu/common/HIPAA/u0746015/marth_software/RUFUS/zenodo_images/rufus_dev.sif}
OUTROOT=${OUTROOT:-$DATA/runs/f7_multi_subject}
THREADS=${SLURM_CPUS_PER_TASK:-8}
WINDOW=10
EXTRA_BIND=${EXTRA_BIND:-}

REF=$FIX/ref/tiny.fa
SOMATIC=$FIX/designed/somatic.vcf     # expected somatic loci
SPLIT=$OUTROOT/_split                 # derived inputs, rebuilt every run

module load apptainer 2>/dev/null || true
command -v apptainer >/dev/null || { echo "ERROR: apptainer unavailable"; exit 1; }
[ -f "$SIF" ] || { echo "ERROR: SIF not found: $SIF"; exit 1; }
for f in "$REF" "$FIX/somatic/tumor.bam" "$FIX/somatic/normal.bam" "$FIX/somatic/normal.cram" "$SOMATIC"; do
  [ -f "$f" ] || { echo "ERROR: fixture missing: $f (run make_fixtures.sh)"; exit 1; }
done

rm -rf "$OUTROOT"; mkdir -p "$SPLIT"
BINDS="$FIX,$OUTROOT,$DATA"; [ -n "$EXTRA_BIND" ] && BINDS="$BINDS,$EXTRA_BIND"

fail() { echo "RESULT: FAIL — $*"; exit 1; }
insif() { apptainer exec --bind "$BINDS" "$SIF" "$@"; }

echo "=== F7 multi-subject input (-s repeated) | $(date) ==="
[ -n "$EXTRA_BIND" ] && echo "  OVERLAY : $EXTRA_BIND  (testing an uncommitted local fix)"
echo

# ---------------------------------------------------------------------------------------------
# Build the split. Reads are assigned by QNAME so both mates land in the same part: the parts are
# a true partition of the original, which is what makes "split == single file" a fair comparison.
# Coordinate order is preserved within each part, but the concatenated stream is no longer globally
# sorted -- RUFUS name-sorts the filtered mates before assembly, so that is expected to wash out.
# ---------------------------------------------------------------------------------------------
echo "-- building split fixtures from tumor.bam --"
# All of the setup runs in ONE container invocation. Container startup dominates this step -- it has
# been seen at ~5 min per exec on a busy node -- so issuing eleven of them made a 9-minute test take
# well over an hour. The work itself is seconds. Written as a file rather than an inline `bash -c`
# so the awk program needs no second level of shell quoting.
cat > "$SPLIT/build_split.sh" <<'SPLITEOF'
#!/bin/bash
# Partition SRC into two BAMs by QNAME (both mates land in the same part), add CRAM copies of each,
# and record read counts so the caller can verify the parts really are a partition.
set -euo pipefail
SRC=$1; REF=$2; OUT=$3
for p in 1 2; do
  samtools view -h "$SRC" \
    | awk -v part="$p" 'BEGIN{OFS="\t"} /^@/{print; next} { if (!($1 in m)) m[$1]=(++n % 2); if (m[$1]==part%2) print }' \
    | samtools view -b -o "$OUT/tumor.part$p.bam" -
  samtools index "$OUT/tumor.part$p.bam"
  samtools view -C -T "$REF" -o "$OUT/tumor.part$p.cram" "$OUT/tumor.part$p.bam"
  samtools index "$OUT/tumor.part$p.cram"
done
{ samtools view -c "$SRC"
  samtools view -c "$OUT/tumor.part1.bam"
  samtools view -c "$OUT/tumor.part2.bam"; } > "$OUT/counts.txt"
SPLITEOF
insif bash "$SPLIT/build_split.sh" "$FIX/somatic/tumor.bam" "$REF" "$SPLIT"

# PARTITION: non-empty parts that sum to the original. A split that silently produced one empty
# file would make the multi-subject run trivially equal to a single-subject run and prove nothing.
N_ALL=$(awk 'NR==1{print $1}' "$SPLIT/counts.txt")
N_P1=$(awk 'NR==2{print $1}' "$SPLIT/counts.txt")
N_P2=$(awk 'NR==3{print $1}' "$SPLIT/counts.txt")
[ -n "$N_ALL" ] && [ -n "$N_P1" ] && [ -n "$N_P2" ] || fail "split build did not produce read counts — see $SPLIT"
echo "   reads: original=$N_ALL  part1=$N_P1  part2=$N_P2  (sum=$((N_P1 + N_P2)))"
[ "$N_P1" -gt 0 ] && [ "$N_P2" -gt 0 ] || fail "split produced an empty part (part1=$N_P1 part2=$N_P2)"
[ "$((N_P1 + N_P2))" -eq "$N_ALL" ] || fail "split is not a partition: $N_P1 + $N_P2 != $N_ALL"

# ---------------------------------------------------------------------------------------------
# run_case <label> <extra runRufus args...>   -> writes $OUTROOT/<label>, sets $CALLS_OUT
#
# Sets a global rather than echoing the path: called as `CALLS=$(run_case ...)` any diagnostic it
# printed -- including fail()'s RESULT line -- would be captured into the variable instead of the
# job log, and `exit 1` would only leave the substitution subshell, so set -e would abort the whole
# script with no message at all. A failing case must say why in the log run_all.sh greps.
# ---------------------------------------------------------------------------------------------
CALLS_OUT=""
run_case() {
  local label=$1; shift
  local out="$OUTROOT/$label" rc=0
  rm -rf "$out"; mkdir -p "$out"
  # Do not let a failing arm abort the script via set -e -- the assertions below give a far better
  # diagnosis than a bare non-zero exit, and a broken multi-subject run is exactly what this catches.
  ( cd "$out"
    apptainer exec --bind "$BINDS" "$SIF" \
      bash /opt/RUFUS/runRufus.sh "$@" -k 25 -L -m 5 -R chr20 -t "$THREADS" -z >"$out/run.log" 2>&1
  ) || rc=$?
  local vcf
  vcf=$(set +f; ls "$out"/temp.RUFUS.Final.*.vcf.gz 2>/dev/null | awk 'NR==1')
  [ -n "$vcf" ] || fail "$label: no VCF (runRufus exit=$rc, status: $(cat "$out/region_status.log" 2>/dev/null)) — see $out/run.log"
  zcat "$vcf" | awk -F'\t' '!/^#/ {print $1"\t"$2"\t"$4"\t"$5}' | sort -k1,1 -k2,2n > "$out/calls.tsv"
  CALLS_OUT="$out/calls.tsv"
}

# Locate the subject generator RUFUS built for a run (named after the FIRST -s file).
subject_generator() {
  local out=$1 first=$2
  set +f; ls "$out/rufus_chr20/$first".*.generator 2>/dev/null | awk 'NR==1'
}

echo
echo "-- run 1/3: baseline, single tumor.bam --"
run_case baseline -s "$FIX/somatic/tumor.bam" -c "$FIX/somatic/normal.bam" -r "$REF"
BASE_CALLS=$CALLS_OUT
echo "   baseline calls:   $(wc -l < "$BASE_CALLS")"

echo "-- run 2/3: multi-subject, tumor.part1.bam + tumor.part2.bam --"
run_case multi_bam \
  -s "$SPLIT/tumor.part1.bam" -s "$SPLIT/tumor.part2.bam" -c "$FIX/somatic/normal.bam" -r "$REF"
MBAM_CALLS=$CALLS_OUT
echo "   multi-bam calls:  $(wc -l < "$MBAM_CALLS")"

echo "-- run 3/3: multi-subject cram (-cr decode path) --"
run_case multi_cram \
  -s "$SPLIT/tumor.part1.cram" -s "$SPLIT/tumor.part2.cram" -c "$FIX/somatic/normal.cram" -cr "$REF"
MCRAM_CALLS=$CALLS_OUT
echo "   multi-cram calls: $(wc -l < "$MCRAM_CALLS")"

# ---- HEADER: the generator must carry every subject but emit exactly one @HD ----
echo
echo "-- generator header discipline and read completeness (root-cause assertions) --"
# Two complementary checks, because they fail on different regressions:
#   @HD count     catches a header emitted per subject -- the generator still carries every record,
#                 but samtools truncates (fastq) or discards (collate) the stream when it hits the
#                 second @HD, so the damage is downstream of the generator.
#   record count  catches a subject file being dropped from the generator outright (a loop or flag
#                 bug, or the header-strip mangling a pre-built .generator command).
# The record count matters because concordance alone would not reliably catch a dropped file HERE:
# each part is ~15x and the planted somatics sit near 0.46 VAF, so one part on its own still clears
# the shipped -m 5 and would very likely produce the same calls. Compare against the single-file
# baseline rather than a literal, so the check tracks the fixture (note the generator's -F 3328
# drops 10 secondary/supplementary records, so this is 40000, not the 40010 `samtools view -c` gives).
BASE_GEN=$(subject_generator "$OUTROOT/baseline" "tumor.bam")
[ -n "$BASE_GEN" ] && [ -f "$BASE_GEN" ] || fail "baseline subject generator not found under $OUTROOT/baseline/rufus_chr20"
BASE_RECS=$(insif bash -c "bash '$BASE_GEN' 2>/dev/null | awk '!/^@/{r++} END{print r+0}'")
[ "${BASE_RECS:-0}" -gt 0 ] || fail "baseline generator emitted no records — cannot calibrate the completeness check"
echo "   baseline    records=$BASE_RECS (single file, the target every multi-subject stream must match)"

for arm in "multi_bam:tumor.part1.bam" "multi_cram:tumor.part1.cram"; do
  label=${arm%%:*}; first=${arm##*:}
  gen=$(subject_generator "$OUTROOT/$label" "$first")
  [ -n "$gen" ] && [ -f "$gen" ] || fail "$label: subject generator not found under $OUTROOT/$label/rufus_chr20 (run used -z, so it should be retained)"
  n_cmds=$(awk 'NF{c++} END{print c+0}' "$gen")
  n_hdr=$(awk '/(^| )-h( |$)/{c++} END{print c+0}' "$gen")
  read -r n_hd n_rec < <(insif bash -c "bash '$gen' 2>/dev/null | awk '/^@HD/{h++} !/^@/{r++} END{print h+0, r+0}'")
  printf "   %-11s commands=%s  header-flags=%s  @HD=%s  records=%s\n" "$label" "$n_cmds" "$n_hdr" "$n_hd" "$n_rec"
  [ "$n_cmds" -eq 2 ] || fail "$label: generator has $n_cmds read commands, expected 2 (one per -s file)"
  [ "$n_hdr" -eq 1 ] || fail "$label: $n_hdr commands carry the header flag, expected exactly 1"
  [ "$n_hd"  -eq 1 ] || fail "$label: generator stream contains $n_hd @HD lines, expected exactly 1 — downstream samtools will reject it"
  [ "$n_rec" -eq "$BASE_RECS" ] || fail "$label: generator stream carries $n_rec records, expected $BASE_RECS — a subject file was dropped from the generator, so its reads never reach RUFUS"
done

# ---- RECALL: the split input must still recover every planted somatic ----
echo
echo "-- recall (multi-subject run must call each planted somatic) --"
recall_ok=1
while IFS=$'\t' read -r _c pos _i _ref _alt _rest; do
  b=$(awk -F'\t' -v p="$pos" -v w="$WINDOW" '$2>=p-w && $2<=p+w{f=1} END{exit !f}' "$MBAM_CALLS"  && echo 1 || echo 0)
  c=$(awk -F'\t' -v p="$pos" -v w="$WINDOW" '$2>=p-w && $2<=p+w{f=1} END{exit !f}' "$MCRAM_CALLS" && echo 1 || echo 0)
  printf "   @%-7s multi-bam=%s multi-cram=%s\n" "$pos" "$b" "$c"
  { [ "$b" = 1 ] && [ "$c" = 1 ]; } || recall_ok=0
done < <(grep -v '^#' "$SOMATIC")
[ "$recall_ok" = 1 ] || fail "a planted somatic was missed by a multi-subject run"

# ---- CONCORDANCE: the same reads, split across files, must give the same calls ----
echo
echo "-- concordance (split call set must equal single-file call set) --"
for arm in multi_bam multi_cram; do
  if diff -q "$BASE_CALLS" "$OUTROOT/$arm/calls.tsv" >/dev/null; then
    echo "   $arm: identical to baseline ($(wc -l < "$BASE_CALLS") calls)"
  else
    echo "   *** DIFF (< baseline only, > $arm only) ***"
    diff "$BASE_CALLS" "$OUTROOT/$arm/calls.tsv" | head -20
    fail "$arm produced different calls than the same reads in a single file"
  fi
done

echo
echo "RESULT: PASS — multi-subject input concordant with single-file ($(wc -l < "$BASE_CALLS") calls), single @HD per generator, all planted somatics recovered"
exit 0
