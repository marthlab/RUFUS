#!/bin/bash
# Works for any number of sample columns in VCF
: "${WORK_DIR:?WORK_DIR must be set}"

input_file=$1

cat "$input_file" | awk -F'\t' '
BEGIN {
    expected_cols = 0
}

# Always print header lines
/^##/ {
    print
    next
}

# Column header: record expected column count and print
/^#CHROM/ {
    expected_cols = NF
    print
    next
}

# Data lines: only print if column count matches
{
    if (expected_cols > 0 && NF == expected_cols) {
        print
    }
}
'
