#!/bin/bash
# F3 — INPUT FORMAT CONCORDANCE (BAM vs CRAM).
#
# runRufus.sh accepts bam / cram / generator for -s and -c (FASTQ is NOT a primary input — it only
# supplements the filter via -q1/-q2). The fixture CRAM is a lossless re-encode of the BAM (same
# reads), so RUFUS must produce IDENTICAL calls from either. CRAM is what the real COLO829 accuracy
# gate uses, so this exercises the -cr decode path end to end.
#
# Runs the somatic scenario twice — once from BAM, once from CRAM — and asserts:
#   CONCORDANCE  the two call sets are byte-identical (CHROM/POS/REF/ALT)
#   RECALL       both recover the planted somatic variants (so concordance isn't "both empty")
#
# Run directly:  bash f3_input_formats.sh
# Or submit:     sbatch f3_input_formats.sh
#SBATCH --account=marth-rw
#SBATCH --partition=marth-rw
#SBATCH --job-name=rufus_f3_formats
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
set -euo pipefail

# sbatch copies this script to the spool dir, so BASH_SOURCE-derived paths break. See F4.
HERE="${FUNCTIONAL_CASES_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}}"
FIX="$(cd "$HERE/../fixtures" 2>/dev/null && pwd)" || { echo "ERROR: cannot locate fixtures from HERE=$HERE"; exit 1; }
DATA=/uufs/chpc.utah.edu/common/HIPAA/u0746015/marth_software/RUFUS/resources/reg_test_files
SIF=${SIF:-/uufs/chpc.utah.edu/common/HIPAA/u0746015/marth_software/RUFUS/zenodo_images/rufus_dev.sif}
OUTROOT=${OUTROOT:-$DATA/runs/f3_formats}
THREADS=${SLURM_CPUS_PER_TASK:-8}
WINDOW=10

REF=$FIX/ref/tiny.fa
SOMATIC=$FIX/designed/somatic.vcf     # expected somatic loci

module load apptainer 2>/dev/null || true
command -v apptainer >/dev/null || { echo "ERROR: apptainer unavailable"; exit 1; }
[ -f "$SIF" ] || { echo "ERROR: SIF not found: $SIF"; exit 1; }
for f in "$REF" "$FIX/somatic/tumor.bam" "$FIX/somatic/tumor.cram" "$FIX/somatic/normal.bam" "$FIX/somatic/normal.cram"; do
  [ -f "$f" ] || { echo "ERROR: fixture missing: $f (run make_fixtures.sh)"; exit 1; }
done

fail() { echo "RESULT: FAIL — $*"; exit 1; }

# run_format <label> <extra runRufus args...>   -> writes $OUTROOT/<label>, echoes the calls.tsv path
run_format() {
  local label=$1; shift
  local out="$OUTROOT/$label"
  rm -rf "$out"; mkdir -p "$out"
  ( cd "$out"
    apptainer exec --bind "$FIX,$out,$DATA" "$SIF" \
      bash /opt/RUFUS/runRufus.sh "$@" -k 25 -L -m 5 -R chr20 -t "$THREADS" -z >"$out/run.log" 2>&1
  )
  local vcf; vcf=$(ls "$out"/temp.RUFUS.Final.*.vcf.gz 2>/dev/null | head -1 || true)
  [ -n "$vcf" ] || fail "$label: no VCF (status: $(cat "$out/region_status.log" 2>/dev/null)) — see $out/run.log"
  zcat "$vcf" | awk -F'\t' '!/^#/ {print $1"\t"$2"\t"$4"\t"$5}' | sort -k1,1 -k2,2n > "$out/calls.tsv"
  echo "$out/calls.tsv"
}

echo "=== F3 input-format concordance (BAM vs CRAM) | $(date) ==="

echo "-- run 1/2: BAM --"
BAM_CALLS=$(run_format bam  -s "$FIX/somatic/tumor.bam"  -c "$FIX/somatic/normal.bam"  -r  "$REF")
echo "   BAM calls:  $(wc -l < "$BAM_CALLS")"

echo "-- run 2/2: CRAM (-cr decode path) --"
CRAM_CALLS=$(run_format cram -s "$FIX/somatic/tumor.cram" -c "$FIX/somatic/normal.cram" -cr "$REF")
echo "   CRAM calls: $(wc -l < "$CRAM_CALLS")"

# ---- RECALL: both formats must recover the planted somatic variants (guards against "both empty") ----
echo "-- recall (both formats must call each planted somatic) --"
recall_ok=1
while IFS=$'\t' read -r _c pos _i _ref _alt _rest; do
  b=$(awk -F'\t' -v p="$pos" -v w="$WINDOW" '$2>=p-w && $2<=p+w{f=1} END{exit !f}' "$BAM_CALLS"  && echo 1 || echo 0)
  c=$(awk -F'\t' -v p="$pos" -v w="$WINDOW" '$2>=p-w && $2<=p+w{f=1} END{exit !f}' "$CRAM_CALLS" && echo 1 || echo 0)
  printf "   @%-7s BAM=%s CRAM=%s\n" "$pos" "$b" "$c"
  { [ "$b" = 1 ] && [ "$c" = 1 ]; } || recall_ok=0
done < <(grep -v '^#' "$SOMATIC")
[ "$recall_ok" = 1 ] || fail "a planted somatic was missed by at least one format"

# ---- CONCORDANCE: identical reads -> identical call sets ----
echo "-- concordance (BAM call set must equal CRAM call set) --"
if diff -q "$BAM_CALLS" "$CRAM_CALLS" >/dev/null; then
  echo "   identical: $(wc -l < "$BAM_CALLS") calls match exactly"
else
  echo "   *** DIFF (< BAM only, > CRAM only) ***"; diff "$BAM_CALLS" "$CRAM_CALLS" | head -20
  fail "BAM and CRAM produced different call sets from identical reads"
fi

echo "RESULT: PASS — BAM and CRAM concordant ($(wc -l < "$BAM_CALLS") calls), both recover the planted somatics"
exit 0
