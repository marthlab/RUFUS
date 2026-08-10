#!/bin/bash
set -euo pipefail

: "${WORK_DIR:?WORK_DIR must be set}"

# Takes in a list of regions to perform a pileup on in bam, returns a single bgzipped pileup vcf
regions=$1
bam=$2
ref=$3

DEPTH=${DEPTH:-100}
bcftools mpileup -Oz -d $DEPTH -f "$ref" -r "${regions}" "$bam"