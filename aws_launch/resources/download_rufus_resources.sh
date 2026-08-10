#!/bin/bash
set -e

DEST_DIR="${1:-./rufus_resources}"

# Constants
MANIFEST_URL=https://s3.us-east-1.amazonaws.com/rufus.marth.lab/public_access_data/rufus_resource_manifest.txt

echo "Downloading Rufus resources to $DEST_DIR..."

if command -v aws &> /dev/null; then
    echo "Using AWS CLI..."
    aws s3 sync s3://rufus.marth.lab/public_access_data/rufus_resources/ "$DEST_DIR" --no-sign-request
elif command -v rclone &> /dev/null; then
    echo "Using rclone..."
	rclone copy :s3:rufus.marth.lab/public_access_data/rufus_resources/ ${DEST_DIR}/rufus_resources/ --s3-provider=AWS --s3-region=us-east-1 --s3-no-check-bucket --s3-env-auth=false -P
else
    # Prompt user for path to file manifest
    echo "Neither AWS CLI nor rclone found. Trying wget."
    FILE_LIST=$(mktemp)
    wget -O "$FILE_LIST" $MANIFEST_URL

    mkdir -p ${DEST_DIR}/rufus_resources
    while IFS= read -r file; do
        wget -P ${DEST_DIR}/rufus_resources https://rufus.marth.lab.s3.us-east-1.amazonaws.com/"$file"
    done < $FILE_LIST 
fi

echo "Download complete!"
