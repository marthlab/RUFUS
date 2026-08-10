#!/bin/bash
set -euo pipefail

: "${WORK_DIR:?WORK_DIR must be set}"

# Calculates and adds an HD_MED info field and allele frequency field to the first column of the vcf file (in RUFUS, this is the tumor/subject) and adds a format field to that same first column with the tag HD_AF. This value is the HD_MED / DP[0]. Prints new file to "hd_af.$IN_VCF".

IN_VCF=$1 # THIS VCF MUST HAVE A REGION SPECIFIC NAME IF RUNNING REGION MODE IN PARALLEL
SUBJECT_SAMPLE_NAME="$2"
FORMATTED_REGION="$3"

if [ -z "${3:-}" ]; then
    echo "ERROR: FORMATTED_REGION (arg 3) is required" >&2
    echo "ERROR - $FORMATTED_REGION is fmtd region and $IN_VCF is in_vcf"
    exit 1
fi

TEMP_FILE="$WORK_DIR/fields.$FORMATTED_REGION.tsv"
TEMP_HD_FILE="$WORK_DIR/hd.$FORMATTED_REGION.tsv"
TEMP_AF_FILE="$WORK_DIR/af.$FORMATTED_REGION.tsv"
TEMP_HDR_FILE="$WORK_DIR/hdr.$FORMATTED_REGION.txt"
IN_VCF_BASENAME="$(basename "$IN_VCF")"
HD_MED_VCF="$WORK_DIR/hd_med.$IN_VCF_BASENAME"
HD_AF_VCF="$WORK_DIR/hd_af.$IN_VCF_BASENAME"

cleanup() {
rm -f "$TEMP_FILE" \
      "${TEMP_AF_FILE}.gz" "${TEMP_AF_FILE}.gz.tbi" \
      "${TEMP_HD_FILE}.gz" "${TEMP_HD_FILE}.gz.tbi" \
      "$HD_MED_VCF" \
      "$TEMP_HDR_FILE"
} 
trap "cleanup" EXIT

# Add HD_MED info field if it doesn't already exist
bcftools query -s "$SUBJECT_SAMPLE_NAME" -f '%CHROM\t%POS\t%REF\t%ALT\t%HD\n' "$IN_VCF" > "$TEMP_FILE"

awk -F'\t' '
function median(arr, n) {
    # Sort the array
    asort(arr)
    if (n % 2 == 1) {
        # If odd, return the middle element
        return arr[int((n + 1) / 2)]
    } else {
        # If even, return the lower of the two middle elements
        return arr[int(n / 2)]
    }
}
{
    # Split the comma-delimited list
    if ($5 == "" || $5 == ".") {
        print $0 "\t0"
        next
    }
    split($5, values, "_")
    count = 0
    # Filter out -1 values and store remaining in new array
    for (i in values) {
        if (values[i] != -1) {
            count++
            filtered[count] = values[i]
        }
    }
    # Calculate the median if there are valid values
    if (count > 0) {
        med = median(filtered, count)
    } else {
        med = 0
    }
    # Print the original line with the median value appended
    print $0 "\t" med
}' "$TEMP_FILE" > "$TEMP_HD_FILE"

bgzip "$TEMP_HD_FILE"

# Index text file
tabix -f -s1 -b2 -e2 "${TEMP_HD_FILE}.gz"

# Make a header line to insert
echo -e '##INFO=<ID=HD_MED,Number=1,Type=Integer,Description="Median of HD array, not including -1s">' > "$TEMP_HDR_FILE"

# Write HD_MED file out
bcftools annotate -s "$SUBJECT_SAMPLE_NAME" -a "${TEMP_HD_FILE}.gz" -h "$TEMP_HDR_FILE" -Oz -c CHROM,POS,REF,ALT,-,HD_MED "$IN_VCF" > "$HD_MED_VCF"

# Pull out fields to text file
bcftools query -s "$SUBJECT_SAMPLE_NAME" -f '%CHROM\t%POS\t%REF\t%ALT\t%HD_MED\t[%DP]\n' "$HD_MED_VCF" > "$TEMP_FILE"

# Calculate AF to 4-digit precision, add as column 7
awk '{ if($6 == 0) printf "%s\t%s\t%s\t%s\t%s\t%s\t%.4f\n", $1, $2, $3, $4, $5, $6, 0; else if($6 < $5) printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $1, $2, $3, $4, $5, $6, "1.00"; else printf "%s\t%s\t%s\t%s\t%s\t%s\t%.4f\n", $1, $2, $3, $4, $5, $6, $5/$6; }' $TEMP_FILE > $TEMP_AF_FILE
bgzip "$TEMP_AF_FILE"

# Index text file
tabix -f -s1 -b2 -e2 "${TEMP_AF_FILE}.gz"

# Make a header line to insert
echo -e '##FORMAT=<ID=HD_AF,Number=1,Type=Float,Description="kMer-based allele frequency for subject sample only (HD_MED/DP)">' > "$TEMP_HDR_FILE"

bcftools annotate -s "$SUBJECT_SAMPLE_NAME" -a "${TEMP_AF_FILE}.gz" -h "$TEMP_HDR_FILE" -Oz -c CHROM,POS,REF,ALT,-,-,FORMAT/HD_AF "$HD_MED_VCF" > "$HD_AF_VCF"