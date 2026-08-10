#!/bin/bash
# F1 — TRIO DE NOVO (RUFUS's flagship path).
#
# child  = germline background + de novo variants (de novo HET, germline HOM)
# mother = germline background only
# father = germline background only
#
# Two assertions, and both matter:
#   RECALL      every planted DE NOVO variant must be CALLED
#   SPECIFICITY every shared GERMLINE variant must NOT be called (it's in both parents,
#               so control subtraction has to remove it)
#
# Unlike F4 this exercises the full pipeline through filter -> assembly -> interpret -> VCF.
#
# Run directly:  bash f1_trio_denovo.sh
# Or submit:     sbatch f1_trio_denovo.sh
#SBATCH --account=marth-rw
#SBATCH --partition=marth-rw
#SBATCH --job-name=rufus_f1_trio
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
OUT=${OUT:-$DATA/runs/f1_trio}
THREADS=${SLURM_CPUS_PER_TASK:-8}
WINDOW=10          # positional tolerance: RUFUS may left-align / represent an MNV or INS differently

# Dev-loop overlay: bind local file(s) over the container to validate a fix BEFORE rebuilding the
# image and round-tripping through CI. Same pattern as editing files mounted over the container.
#   EXTRA_BIND=/path/to/repo/scripts/Foo.pl:/opt/RUFUS/scripts/Foo.pl sbatch f1_trio_denovo.sh
EXTRA_BIND=${EXTRA_BIND:-}

REF=$FIX/ref/tiny.fa
CHILD=$FIX/trio/child.bam
MOM=$FIX/trio/mother.bam
DAD=$FIX/trio/father.bam
DENOVO=$FIX/designed/denovo.vcf       # expected TO be called
GERMLINE=$FIX/designed/germline.vcf   # expected NOT to be called

rm -rf "$OUT"; mkdir -p "$OUT"; cd "$OUT"
module load apptainer 2>/dev/null || true
command -v apptainer >/dev/null || { echo "ERROR: apptainer unavailable"; exit 1; }
[ -f "$SIF" ] || { echo "ERROR: SIF not found: $SIF"; exit 1; }
for f in "$REF" "$CHILD" "$MOM" "$DAD" "$DENOVO" "$GERMLINE"; do
  [ -f "$f" ] || { echo "ERROR: fixture missing: $f (run make_fixtures.sh)"; exit 1; }
done

echo "=== F1 trio de novo | $(date) ==="
echo "  subject  : child.bam        (germline HOM + de novo HET)"
echo "  controls : mother.bam, father.bam"
echo "  expect   : $(grep -vc '^#' "$DENOVO") de novo CALLED, $(grep -vc '^#' "$GERMLINE") germline NOT called"
echo

BINDS="$FIX,$OUT,$DATA"
if [ -n "$EXTRA_BIND" ]; then
  BINDS="$BINDS,$EXTRA_BIND"
  echo "  OVERLAY : $EXTRA_BIND  (testing an uncommitted local fix)"
  echo
fi

set +e
apptainer exec --bind "$BINDS" "$SIF" \
  bash /opt/RUFUS/runRufus.sh \
    -s "$CHILD" -c "$MOM" -c "$DAD" -r "$REF" \
    -k 25 -L -m 5 -R chr20 -t "$THREADS" -z
    # -L is REQUIRED here and is the shipped default (setup_slurm hardcodes `-k 25 -L -vs`).
    # RUFUS labels these de novo calls X-Mosaic/5I-Mosaic/3X-Mosaic; without -L the post-filter
    # VilterAutosomeOnly.withoutMosaic strips every one of them and the VCF comes back empty
    # even though RUFUS.interpret found all three at the correct positions.
RC=$?
set -e
echo "runRufus exit: $RC"

STATUS=$(cat "$OUT/region_status.log" 2>/dev/null || echo "<none>")
echo "region_status: $STATUS"
WORK="$OUT/rufus_chr20"
fail() { echo "RESULT: FAIL — $*"; exit 1; }

# ---- preconditions: prove the pipeline actually ran on ALL THREE samples ----
echo "-- preconditions --"
for s in child mother father; do
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
[ -n "$VCF" ] || fail "no VCF emitted (status: $STATUS) — expected de novo calls"
CALLS=$OUT/calls.tsv
zcat "$VCF" | awk -F'\t' '!/^#/ {print $1"\t"$2"\t"$4"\t"$5}' > "$CALLS"
NCALL=$(wc -l < "$CALLS")
echo "-- calls: $NCALL total (VCF: $(basename "$VCF")) --"

# helper: is there a call within +-WINDOW of position $1 ?
called_near() { awk -F'\t' -v p="$1" -v w="$WINDOW" '$2>=p-w && $2<=p+w {f=1} END{exit !f}' "$CALLS"; }

# ---- RECALL: every de novo must be called ----
echo "-- recall (de novo, must be CALLED) --"
MISS=0
while IFS=$'\t' read -r _c pos _i ref alt _rest; do
  if called_near "$pos"; then
    echo "   de novo @$pos ($ref>$alt)  CALLED"
  else
    echo "   de novo @$pos ($ref>$alt)  *** MISSED ***"; MISS=$((MISS+1))
  fi
done < <(grep -v '^#' "$DENOVO")

# ---- SPECIFICITY: no germline may be called ----
echo "-- specificity (germline, must NOT be called) --"
FP=0
while IFS=$'\t' read -r _c pos _i ref alt _rest; do
  if called_near "$pos"; then
    echo "   germline @$pos ($ref>$alt)  *** FALSELY CALLED ***"; FP=$((FP+1))
  else
    echo "   germline @$pos ($ref>$alt)  correctly absent"
  fi
done < <(grep -v '^#' "$GERMLINE")

# ---- PRECISION: every call must be attributable to a planted de novo locus ----
# NB: record count > planted count is EXPECTED and fine — bcftools norm atomizes a planted MNV
# into one record per base (ATC>CGA becomes A>C, T>G, C>A sharing one ID). What must not happen
# is a call landing somewhere we planted nothing.
echo "-- precision (no calls outside a planted locus) --"
UNEXP=0
while IFS=$'\t' read -r _c pos _r _a; do
  ok=0
  while read -r dp; do
    if [ "$pos" -ge $((dp-WINDOW)) ] && [ "$pos" -le $((dp+WINDOW)) ]; then ok=1; break; fi
  done < <(grep -v '^#' "$DENOVO" | cut -f2)
  if [ "$ok" -eq 0 ]; then echo "   *** UNEXPECTED call @$pos ***"; UNEXP=$((UNEXP+1)); fi
done < "$CALLS"
[ "$UNEXP" -eq 0 ] && echo "   all $NCALL call(s) map to a planted locus  OK"

NDN=$(grep -vc '^#' "$DENOVO")
echo
echo "-- summary --"
echo "   unexpected calls : $UNEXP"
echo "   de novo recalled : $((NDN-MISS))/$NDN"
echo "   germline leaked  : $FP"
echo "   total calls      : $NCALL"

[ "$MISS" -eq 0 ]  || fail "$MISS de novo variant(s) missed"
[ "$FP" -eq 0 ]    || fail "$FP germline variant(s) leaked into the call set"
[ "$UNEXP" -eq 0 ] || fail "$UNEXP call(s) at loci where nothing was planted (false positives)"
echo "RESULT: PASS — all $NDN de novo called, no germline leaked"
exit 0
