#!/bin/bash
# Downloaded hashes get put in ${WORKING_DIR}/rufus_temp/downloaded_{hash_type}_hashes
# File name depends only upon matching *wg* or *fmtd_region* in the correct type/version directory

HASH_TYPE="$1"
HASH_VERSION="$2"
FMTD_REGION="$3" # Must be "wg" if not a region

s3_path="s3://rufus.marth.lab/public_access_data/rufus_resources/${HASH_TYPE}_hashes/${HASH_VERSION}/"

# Check that we only have one file to download first
matches=$(aws s3 ls --no-sign-request $s3_path | grep "$FMTD_REGION" | awk '{print $NF}')
count=$(echo "$matches" | grep -c "^")

if [ -z "$matches" ]; then
    count=0
else
    count=$(echo "$matches" | wc -l)
fi

if [ $count -eq 0 ]; then
    echo "Error: No files found matching pattern '$FMTD_REGION'" >&2
    exit 1
elif [ $count -gt 1 ]; then
    echo "Error: Multiple files found matching pattern '$FMTD_REGION':" >&2
    echo "$matches" >&2
    exit 1
else
    aws s3 cp --no-sign-request "${s3_path}${matches}" "${FMTD_REGION}_${HASH_TYPE}_${HASH_VERSION}.Jhash"
fi