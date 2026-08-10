#!/bin/bash
# F6 — SINGLE-END FASTQ INPUT (-se, whole-genome).
#
# Guards two things that were both broken until 2026-07-22:
#   1. the run COMPLETES (does not hang) — the single-end filter used to deadlock on a FIFO fd;
#   2. interpret RECOVERS every planted somatic — asserted at the RAW interpret VCF
#      (*.V2.overlap.hashcount.fastq.bam.vcf), NOT the final VCF.
#
# Why the raw VCF: single-end calls are currently tagged StrandBias (SB) by RUFUS.interpret and
# dropped by the post-filter, so the FINAL VCF is empty. That strand-bias behaviour is a SEPARATE,
# known issue with its own fix — this test deliberately asserts upstream of it so it stays meaningful.
# WHEN strand bias is fixed: tighten this to assert the final VCF (see F5) and drop the raw-VCF check.
#
# Genuine single-end input = R1 reads only (one read per fragment).
# EXTRA_BIND overlays uncommitted code before the image is rebuilt (harmless no-op afterward).
#
# Run directly:  bash f6_fastq_single_end.sh   |   Submit:  sbatch f6_fastq_single_end.sh
#SBATCH --account=marth-rw
#SBATCH --partition=marth-rw
#SBATCH --job-name=rufus_f6_fastq_se
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=00:20:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
set -euo pipefail

HERE="${FUNCTIONAL_CASES_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}}"
FIX="$(cd "$HERE/../fixtures" 2>/dev/null && pwd)" || { echo "ERROR: cannot locate fixtures from HERE=$HERE"; exit 1; }
DATA=/uufs/chpc.utah.edu/common/HIPAA/u0746015/marth_software/RUFUS/resources/reg_test_files
SIF=${SIF:-/uufs/chpc.utah.edu/common/HIPAA/u0746015/marth_software/RUFUS/zenodo_images/rufus_dev.sif}
OUT=${OUT:-$DATA/runs/f6_fastq_se}
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
for f in "$REF" "$FIX/fastq/tumor.1.fq.gz" "$FIX/fastq/normal.1.fq.gz"; do
  [ -f "$f" ] || fail "fixture missing: $f (run make_fixtures.sh)"
done

BINDS="$FIX,$OUT,$DATA"; [ -n "$EXTRA_BIND" ] && BINDS="$BINDS,$EXTRA_BIND"

echo "=== F6 single-end FASTQ (-se, whole-genome) | $(date) ==="
# Bounded walltime is part of the test: a FIFO-hang regression shows up as a TIMEOUT, not a pass.
set +e
apptainer exec --bind "$BINDS" "$SIF" \
  bash /opt/RUFUS/runRufus.sh \
    -s "$FIX/fastq/tumor.1.fq.gz" -c "$FIX/fastq/normal.1.fq.gz" \
    -r "$REF" -se -k 25 -L -m 5 -hs 1G -t "$THREADS" -z
RC=$?
set -e
echo "runRufus exit: $RC   region_status: $(cat "$OUT/region_status.log" 2>/dev/null)"

WORK="$OUT/rufus_wg"
# 1) completion + real filtering
MUT=$(ls "$WORK"/*.Mutations.fastq 2>/dev/null | head -1)
[ -n "$MUT" ] && [ -s "$MUT" ] || fail "no Mutations.fastq — single-end filter produced nothing (hang/regression?)"
echo "   mutant reads kept: $(( $(wc -l < "$MUT") / 4 ))"

# 2) recall at the RAW interpret VCF (strand-bias-independent)
RAW=$(ls "$WORK"/*.V2.overlap.hashcount.fastq.bam.vcf 2>/dev/null | head -1)
[ -n "$RAW" ] && [ -s "$RAW" ] || fail "no raw interpret VCF — pipeline did not reach interpret"
grep -v '^#' "$RAW" | awk -F'\t' '{print $2}' | sort -n > "$OUT/rawpos.txt"
echo "-- raw interpret calls: $(wc -l < "$OUT/rawpos.txt") --"

MISS=0
while IFS=$'\t' read -r _c pos _r; do
  if awk -v p="$pos" -v w="$WINDOW" '$1>=p-w && $1<=p+w{f=1} END{exit !f}' "$OUT/rawpos.txt"; then
    echo "   somatic @$pos found by interpret"
  else echo "   somatic @$pos *** MISSED by interpret ***"; MISS=1; fi
done < <(grep -v '^#' "$SOMATIC" | awk -F'\t' '{print $1"\t"$2"\t"$4}')
[ "$MISS" -eq 0 ] || fail "interpret missed a planted somatic in single-end mode"

FINAL=$(ls "$OUT"/temp.RUFUS.Final.*.vcf.gz 2>/dev/null | head -1 || true)
[ -n "$FINAL" ] && echo "note: final VCF has $(zcat "$FINAL" | grep -vc '^#') calls (empty expected until strand-bias fix)"
echo "RESULT: PASS — single-end FASTQ completes and interpret recovers all planted somatics"
