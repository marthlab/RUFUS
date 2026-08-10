#!/bin/bash
# F2 — TUMOR/NORMAL SOMATIC (single-control path; the COLO829-shaped test on fixtures).
#
# tumor  = germline background + somatic variants (somatic HET ~0.5, germline HOM)
# normal = germline background only
#
# Two assertions:
#   RECALL      every planted SOMATIC variant must be CALLED (incl. the 1kb deletion — the only
#               fixture variant exercising the large-SV assembly path)
#   SPECIFICITY every shared GERMLINE variant must NOT be called (present in the normal control,
#               so subtraction must remove it)
#
# Differs from F1 only in configuration: ONE control instead of a trio.
#
# Run directly:  bash f2_somatic.sh
# Or submit:     sbatch f2_somatic.sh
#SBATCH --account=marth-rw
#SBATCH --partition=marth-rw
#SBATCH --job-name=rufus_f2_somatic
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
set -euo pipefail

# sbatch copies this script to the spool dir, so BASH_SOURCE-derived paths break. See F4.
HERE="${FUNCTIONAL_CASES_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}}"
FIX="$(cd "$HERE/../fixtures" 2>/dev/null && pwd)" || {
  echo "ERROR: cannot locate fixtures from HERE=$HERE"; exit 1; }
DATA=/uufs/chpc.utah.edu/common/HIPAA/u0746015/marth_software/RUFUS/resources/reg_test_files
SIF=${SIF:-/uufs/chpc.utah.edu/common/HIPAA/u0746015/marth_software/RUFUS/zenodo_images/rufus_dev.sif}
OUT=${OUT:-$DATA/runs/f2_somatic}
THREADS=${SLURM_CPUS_PER_TASK:-8}
WINDOW=10          # positional tolerance: RUFUS may left-align / atomize a call vs the planted pos

# Dev-loop overlay: bind local file(s) over the container to validate a fix before a CI rebuild.
#   EXTRA_BIND=/path/to/repo/scripts/Foo.pl:/opt/RUFUS/scripts/Foo.pl sbatch f2_somatic.sh
EXTRA_BIND=${EXTRA_BIND:-}

REF=$FIX/ref/tiny.fa
TUMOR=$FIX/somatic/tumor.bam
NORMAL=$FIX/somatic/normal.bam
SOMATIC=$FIX/designed/somatic.vcf     # expected TO be called
GERMLINE=$FIX/designed/germline.vcf   # expected NOT to be called

rm -rf "$OUT"; mkdir -p "$OUT"; cd "$OUT"
module load apptainer 2>/dev/null || true
command -v apptainer >/dev/null || { echo "ERROR: apptainer unavailable"; exit 1; }
[ -f "$SIF" ] || { echo "ERROR: SIF not found: $SIF"; exit 1; }
for f in "$REF" "$TUMOR" "$NORMAL" "$SOMATIC" "$GERMLINE"; do
  [ -f "$f" ] || { echo "ERROR: fixture missing: $f (run make_fixtures.sh)"; exit 1; }
done

echo "=== F2 tumor/normal somatic | $(date) ==="
echo "  subject : tumor.bam        (germline HOM + somatic HET)"
echo "  control : normal.bam       (germline only)"
echo "  expect  : $(grep -vc '^#' "$SOMATIC") somatic CALLED, $(grep -vc '^#' "$GERMLINE") germline NOT called"
echo

BINDS="$FIX,$OUT,$DATA"
if [ -n "$EXTRA_BIND" ]; then
  BINDS="$BINDS,$EXTRA_BIND"
  echo "  OVERLAY : $EXTRA_BIND  (testing an uncommitted local fix)"; echo
fi

set +e
apptainer exec --bind "$BINDS" "$SIF" \
  bash /opt/RUFUS/runRufus.sh \
    -s "$TUMOR" -c "$NORMAL" -r "$REF" \
    -k 25 -L -m 5 -R chr20 -t "$THREADS" -z
    # -L is the shipped default (setup_slurm hardcodes `-k 25 -L -vs`); somatic HET calls are
    # labelled *-Mosaic and are stripped by VilterAutosomeOnly.withoutMosaic when -L is absent.
RC=$?
set -e
echo "runRufus exit: $RC"

STATUS=$(cat "$OUT/region_status.log" 2>/dev/null || echo "<none>")
echo "region_status: $STATUS"
WORK="$OUT/rufus_chr20"
fail() { echo "RESULT: FAIL — $*"; exit 1; }

# ---- preconditions: prove the pipeline actually ran on BOTH samples ----
echo "-- preconditions --"
for s in tumor normal; do
  JH="$WORK/$s.bam.chr20.generator.Jhash"; HI="$JH.histo"
  [ -s "$JH" ] || fail "$s: Jhash missing/empty — jellyfish did not run (OOM signature)"
  [ -s "$HI" ] || fail "$s: histo missing/empty"
  TOT=$(awk '{s+=$2} END{print s+0}' "$HI")
  [ "$TOT" -ge 100000 ] || fail "$s: only $TOT distinct kmers — counting did not really run"
  MODE=$(awk 'BEGIN{m=-1;c=0} {if($2+0>m){m=$2+0;c=$1+0}} END{print c}' "$HI")
  echo "   $s: distinct=$TOT modal_depth=$MODE  OK"
done

# ---- collect calls ----
VCF=$(ls "$OUT"/temp.RUFUS.Final.*.vcf.gz "$OUT"/RUFUS.Final.*.vcf.gz 2>/dev/null | head -1 || true)
[ -n "$VCF" ] || fail "no VCF emitted (status: $STATUS) — expected somatic calls"
CALLS=$OUT/calls.tsv
zcat "$VCF" | awk -F'\t' '!/^#/ {print $1"\t"$2"\t"$4"\t"$5}' > "$CALLS"
NCALL=$(wc -l < "$CALLS")
echo "-- calls: $NCALL total (VCF: $(basename "$VCF")) --"

called_near() { awk -F'\t' -v p="$1" -v w="$WINDOW" '$2>=p-w && $2<=p+w {f=1} END{exit !f}' "$CALLS"; }

# ---- RECALL: every somatic must be called ----
echo "-- recall (somatic, must be CALLED) --"
MISS=0
while IFS=$'\t' read -r _c pos _i ref alt _rest; do
  klass=$(( ${#ref} > 50 || ${#alt} > 50 ? 1 : 0 ))  # flag the large event for the log
  tag=$([ "$klass" -eq 1 ] && echo " [large SV]" || echo "")
  if called_near "$pos"; then
    echo "   somatic @$pos (${ref:0:8}>${alt:0:8})$tag  CALLED"
  else
    echo "   somatic @$pos (${ref:0:8}>${alt:0:8})$tag  *** MISSED ***"; MISS=$((MISS+1))
  fi
done < <(grep -v '^#' "$SOMATIC")

# ---- SPECIFICITY: no germline may be called ----
echo "-- specificity (germline, must NOT be called) --"
FP=0
while IFS=$'\t' read -r _c pos _i ref alt _rest; do
  if called_near "$pos"; then
    echo "   germline @$pos (${ref:0:8}>${alt:0:8})  *** FALSELY CALLED ***"; FP=$((FP+1))
  else
    echo "   germline @$pos (${ref:0:8}>${alt:0:8})  correctly absent"
  fi
done < <(grep -v '^#' "$GERMLINE")

# ---- PRECISION: every call must map to a planted somatic locus ----
# Record count > planted count is EXPECTED (bcftools norm atomizes an MNV into per-base records).
echo "-- precision (no calls outside a planted locus) --"
UNEXP=0
while IFS=$'\t' read -r _c pos _r _a; do
  ok=0
  while read -r sp; do
    if [ "$pos" -ge $((sp-WINDOW)) ] && [ "$pos" -le $((sp+WINDOW)) ]; then ok=1; break; fi
  done < <(grep -v '^#' "$SOMATIC" | cut -f2)
  if [ "$ok" -eq 0 ]; then echo "   *** UNEXPECTED call @$pos ***"; UNEXP=$((UNEXP+1)); fi
done < "$CALLS"
[ "$UNEXP" -eq 0 ] && echo "   all $NCALL call(s) map to a planted locus  OK"

# ---- LARGE-SV SIZE: the 1kb deletion is the only variant exercising SV assembly. Position alone
# would pass even if RUFUS miscalled it as a small indel, so assert the event size too. ----
echo "-- large-SV size (planted 1000bp deletion @170000) --"
SVDELTA=$(zcat "$VCF" | awk -F'\t' '!/^#/ && $2>=169990 && $2<=170010 {
  d=length($4)-length($5); if(d<0)d=-d; if(d>m)m=d} END{print m+0}')
echo "   largest event near 170000: ${SVDELTA}bp (expect ~1000)"
SVOK=$([ "$SVDELTA" -ge 900 ] && echo 1 || echo 0)

NSOM=$(grep -vc '^#' "$SOMATIC")
echo
echo "-- summary --"
echo "   unexpected calls : $UNEXP"
echo "   somatic recalled : $((NSOM-MISS))/$NSOM"
echo "   germline leaked  : $FP"
echo "   total calls      : $NCALL"

[ "$MISS" -eq 0 ]  || fail "$MISS somatic variant(s) missed"
[ "$FP" -eq 0 ]    || fail "$FP germline variant(s) leaked into the call set"
[ "$UNEXP" -eq 0 ] || fail "$UNEXP call(s) at loci where nothing was planted (false positives)"
[ "$SVOK" -eq 1 ]  || fail "large deletion mis-sized: ${SVDELTA}bp near 170000, expected ~1000 (SV assembly regressed?)"
echo "RESULT: PASS — all $NSOM somatic called (incl. ${SVDELTA}bp SV), no germline leaked"
exit 0
