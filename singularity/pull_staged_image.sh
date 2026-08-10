#!/bin/bash
# Pull a published RUFUS image from Docker Hub into a SIF on HPC (CHPC).
#
# Replaces the old manual `sudo singularity build rufus.sif rufus.def` flow: the image is now
# built and published by CI, and apptainer converts it to a SIF directly from Docker Hub.
#
# Usage:
#   pull_staged_image.sh [tag] [dest_dir]
#     tag       Docker Hub tag to pull (default: stage). Use a version like v1.1.11 for prod.
#     dest_dir  Directory to write the SIF into (default: the CHPC zenodo_images dir).
set -euo pipefail

IMAGE="docker://stefinfection/rufus"
TAG="${1:-stage}"
DEST_DIR="${2:-/uufs/chpc.utah.edu/common/HIPAA/u0746015/marth_software/RUFUS/zenodo_images}"

command -v apptainer >/dev/null 2>&1 || { echo "ERROR: apptainer not found (try: module load apptainer)"; exit 1; }
mkdir -p "$DEST_DIR"

SIF_PATH="${DEST_DIR}/rufus_${TAG}.sif"
echo "Pulling ${IMAGE}:${TAG} -> ${SIF_PATH}"
apptainer pull --force "$SIF_PATH" "${IMAGE}:${TAG}"

echo "Verifying image..."
apptainer exec "$SIF_PATH" bash /opt/RUFUS/tests/smoke_test.sh

echo "Done: ${SIF_PATH}"
