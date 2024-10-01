#!/bin/bash

chr=$1
bam=$2
ref=$3
start_coord=$4
end_coord=$5

DEPTH=500

if [ ! -z "$start_coord" ] && [ ! -z "$end_coord" ]; then
  bcftools mpileup -Ov -d $DEPTH -f $ref -r "${chr}:${start_coord}-${end_coord}" -o mpileup_${chr}_${start_coord}_${end_coord}.vcf $bam
else
  bcftools mpileup -Ov -d $DEPTH -f $ref -r "${chr}" -o mpileup_${chr}.vcf $bam
fi