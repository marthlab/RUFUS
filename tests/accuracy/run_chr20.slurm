#!/bin/bash
# tests/accuracy/run_chr20.slurm
#
# chr20 accuracy run for the pre-release gate: full-strength COLO829 tumor vs normal in the shipped
# configuration. Direct `apptainer exec runRufus.sh` (isolates the CALLER; orchestration is separate).
# Data paths + image come from config.sh (copy config.sh.example -> config.sh first). SUBMIT FROM THIS DIR:
#
#     cd tests/accuracy && sbatch run_chr20.slurm        # then, when it finishes:  bash compare.sh
#
# EDIT --account / --partition below for your cluster — SBATCH headers cannot read config.sh.
#SBATCH --account=marth-rw
#SBATCH --partition=marth-rw
#SBATCH --job-name=rufus_chr20_accuracy
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=96G
#SBATCH --time=1-00:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
set -euo pipefail

# sbatch copies this script to the spool dir, so find config.sh via SLURM_SUBMIT_DIR (see functional/cases).
HERE="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
CFG="$HERE/config.sh"
[ -f "$CFG" ] || { echo "ERROR: missing $CFG — copy config.sh.example to config.sh and edit it (README.md)"; exit 1; }
# shellcheck source=/dev/null
source "$CFG"

THREADS=${SLURM_CPUS_PER_TASK:-32}
OUT=$RUFUS_TEST_DATA/runs/chr20_m5          # large run outputs live in the external data area, NOT the repo
# Start from a clean work dir: RUFUS resumes by skipping stages whose output exists, so a dirty dir
# silently reuses stale intermediates from a prior run (that hung a chr20 job for hours).
rm -rf "$OUT"; mkdir -p "$OUT"; cd "$OUT"

module load apptainer 2>/dev/null || true
command -v apptainer >/dev/null || { echo "ERROR: apptainer not available (module load apptainer)"; exit 1; }
for f in "$RUFUS_SIF" "$REF_FASTA" "$TUMOR_CRAM" "$NORMAL_CRAM"; do
  [ -e "$f" ] || { echo "ERROR: missing $f  (check config.sh / see README.md)"; exit 1; }
done

BINDS="$RUFUS_TEST_DATA"

# preflight: prove the reference decodes the CRAMs inside the container (fails in seconds, not hours)
for c in "$TUMOR_CRAM" "$NORMAL_CRAM"; do
  n=$(apptainer exec --bind "$BINDS" "$RUFUS_SIF" samtools view -T "$REF_FASTA" "$c" ${REGION}:5000000-5000200 2>/dev/null | wc -l)
  [ "$n" -gt 0 ] || { echo "PREFLIGHT FAILED: 0 reads decoded from $(basename "$c") — reference/bind is wrong"; exit 1; }
  echo "preflight OK: $(basename "$c") -> $n reads"
done

echo "=== RUFUS chr20 accuracy | image=$(basename "$RUFUS_SIF") | $(date) ==="
# -vs (very-fast assembly) = the shipped setup_slurm config. Full assembly DEADLOCKS on real chr20 data.
# Do NOT add -hs: the whole-genome value 64G OOMs a windowed run; the 1G windowed default is correct.
apptainer exec --bind "$BINDS" "$RUFUS_SIF" \
  bash /opt/RUFUS/runRufus.sh \
    -s "$TUMOR_CRAM" -c "$NORMAL_CRAM" -cr "$REF_FASTA" \
    -k 25 -L -vs -m 5 -z \
    -R "$REGION" -t "$THREADS"

echo "=== done -> $OUT  $(date) ==="
echo "next: cd $HERE && bash compare.sh"
