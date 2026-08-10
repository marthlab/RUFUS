#!/bin/bash
host_ref_file="$1"

# ENV override
: "${RUFUS_ROOT:=/opt/RUFUS}"

bwa="$RUFUS_ROOT/bin/externals/bwa/src/bwa_project/bwa"

if [ -z "$host_ref_file" ]; then
    echo "ERROR: usage: build_bwa_indexes.sh <reference fasta>" >&2
    exit 1
fi

# Handle .gz extension
if [[ "$host_ref_file" == *.gz ]]; then
    host_ref_file_no_gz="${host_ref_file%.gz}"
else
    host_ref_file_no_gz="$host_ref_file"
fi

fasta_idx="${host_ref_file_no_gz}"

if [ ! -f "$fasta_idx" ]; then
    echo "ERROR: reference $fasta_idx does not exist or cannot be read." >&2
    if [ "$host_ref_file" != "$fasta_idx" ]; then
        echo "       BWA cannot index a compressed reference; decompress $host_ref_file first." >&2
    fi
    exit 1
fi

# Index the file
$bwa index -a bwtsw "$fasta_idx" || { echo "ERROR: bwa index failed on $fasta_idx" >&2; exit 1; }
samtools faidx "$fasta_idx" || { echo "ERROR: samtools faidx failed on $fasta_idx" >&2; exit 1; }

echo "These indexes were created by RUFUS for an intermediate BWA step." > README.md
echo "If you wish to reuse them for the next run to save some time, copy them into the same directory as your REFERENCE_FASTA." >> README.md