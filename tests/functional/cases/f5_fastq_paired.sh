#!/bin/bash
# F5 — PAIRED-END FASTQ INPUT (whole-genome).
#
# FASTQ is unaligned so it runs whole-genome only (no -R). Paired filtering goes through the -q1/-q2
# direct path (auto-wired from the two -s files). Same underlying reads as the somatic BAM/CRAM
# fixtures, so RUFUS must recover the same planted somatic variants.
#
#   RECALL     every planted somatic is called in the FINAL VCF
#   PRECISION  no call lands outside a planted locus
#
# EXTRA_BIND=host:container[,host:container...] overlays uncommitted code (dev loop) before the image
# is rebuilt; it is a harmless no-op once the FASTQ support is in the image.
#
# Run directly:  bash f5_fastq_paired.sh   |   Submit:  sbatch f5_fastq_paired.sh
#SBATCH --account=marth-rw
#SBATCH --partition=marth-rw
#SBATCH --job-name=rufus_f5_fastq_paired
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
set -euo pipefail

HERE="${FUNCTIONAL_CASES_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}}"
FIX="$(cd "$HERE/../fixtures" 2>/dev/null && pwd)" || { echo "ERROR: cannot locate fixtures from HERE=$HERE"; exit 1; }
DATA=/uufs/chpc.utah.edu/common/HIPAA/u0746015/marth_software/RUFUS/resources/reg_test_files
SIF=${SIF:-/uufs/chpc.utah.edu/common/HIPAA/u0746015/marth_software/RUFUS/zenodo_images/rufus_dev.sif}
OUT=${OUT:-$DATA/runs/f5_fastq_paired}
THREADS=${SLURM_CPUS_PER_TASK:-8}
WINDOW=10
EXTRA_BIND=${EXTRA_BIND:-}

REF=$FIX/ref/tiny.fa
SOMATIC=$FIX/designed/somatic.vcf
fail() { echo "RESULT: FAIL — $*"; exit 1; }

rm -rf "$OUT"; mkdir -p "$OUT"; cd "$OUT"
module load apptainer 2>/dev/null || true
command -v apptainer >/dev/null || { echo "ERROR: apptainer unavailable"; exit 1; }
[ -f "$SIF" ] || { echo "ERROR: SIF not found: $SIF"; exit 1; }
for f in "$REF" "$FIX/fastq/tumor.1.fq.gz" "$FIX/fastq/tumor.2.fq.gz" "$FIX/fastq/normal.1.fq.gz" "$FIX/fastq/normal.2.fq.gz"; do
  [ -f "$f" ] || fail "fixture missing: $f (run make_fixtures.sh)"
done

BINDS="$FIX,$OUT,$DATA"; [ -n "$EXTRA_BIND" ] && BINDS="$BINDS,$EXTRA_BIND"

echo "=== F5 paired-end FASTQ (whole-genome) | $(date) ==="
# no -R (fastq is whole-genome); -hs 1G because WG default is 64G (237GB), absurd for a 200kb fixture.
set +e
apptainer exec --bind "$BINDS" "$SIF" \
  bash /opt/RUFUS/runRufus.sh \
    -s "$FIX/fastq/tumor.1.fq.gz"  "$FIX/fastq/tumor.2.fq.gz" \
    -c "$FIX/fastq/normal.1.fq.gz" "$FIX/fastq/normal.2.fq.gz" \
    -r "$REF" -k 25 -L -m 5 -hs 1G -t "$THREADS" -z
RC=$?
set -e
echo "runRufus exit: $RC   region_status: $(cat "$OUT/region_status.log" 2>/dev/null)"

# preconditions: both samples really counted (WG generators are named *.wg.generator)
WORK="$OUT/rufus_wg"
for s in tumor normal; do
  H=$(ls "$WORK/$s".*.wg.generator.Jhash.histo 2>/dev/null | head -1)
  [ -n "$H" ] && [ -s "$H" ] || fail "$s: no k-mer histo — counting did not run"
  TOT=$(awk '{s+=$2} END{print s+0}' "$H")
  [ "$TOT" -ge 100000 ] || fail "$s: only $TOT distinct kmers — counting did not really run"
  echo "   $s counted: $TOT distinct kmers"
done

VCF=$(ls "$OUT"/temp.RUFUS.Final.*.vcf.gz "$OUT"/RUFUS.Final.*.vcf.gz 2>/dev/null | head -1 || true)
[ -n "$VCF" ] || fail "no final VCF (status: $(cat "$OUT/region_status.log" 2>/dev/null))"
zcat "$VCF" | awk -F'\t' '!/^#/{print $2}' | sort -n > "$OUT/calledpos.txt"
NCALL=$(wc -l < "$OUT/calledpos.txt")
echo "-- final VCF: $NCALL calls --"

MISS=0
while IFS=$'\t' read -r _c pos _r; do
  if awk -v p="$pos" -v w="$WINDOW" '$1>=p-w && $1<=p+w{f=1} END{exit !f}' "$OUT/calledpos.txt"; then
    echo "   somatic @$pos CALLED"
  else echo "   somatic @$pos *** MISSED ***"; MISS=1; fi
done < <(grep -v '^#' "$SOMATIC" | awk -F'\t' '{print $1"\t"$2"\t"$4}')

UNEXP=0
while read -r cp; do
  ok=0; while read -r sp; do { [ "$cp" -ge $((sp-WINDOW)) ] && [ "$cp" -le $((sp+WINDOW)) ]; } && { ok=1; break; }; done < <(grep -v '^#' "$SOMATIC" | cut -f2)
  [ "$ok" -eq 0 ] && { echo "   *** UNEXPECTED call @$cp ***"; UNEXP=1; }
done < "$OUT/calledpos.txt"

[ "$MISS" -eq 0 ]  || fail "a planted somatic was missed"
[ "$UNEXP" -eq 0 ] || fail "a call landed where nothing was planted"
echo "RESULT: PASS — paired-end FASTQ (WG) recovered all planted somatics, no spurious calls"
