#!/bin/bash

# Script to generate a manifest file for Rufus resources
# This creates a user-friendly file listing with sizes and download URLs

BUCKET="rufus.marth.lab"
PREFIX="public_access_data/rufus_resources/"
REGION="us-east-1"
OUTPUT_DIR="."

echo "Generating manifest for s3://${BUCKET}/${PREFIX}..."

# Generate basic file listing with sizes
echo "Creating file listing..."
aws s3 ls s3://${BUCKET}/${PREFIX} --recursive --no-sign-request --human-readable > "${OUTPUT_DIR}/file_listing_human.txt"

# Generate machine-readable listing (bytes)
aws s3 ls s3://${BUCKET}/${PREFIX} --recursive --no-sign-request > "${OUTPUT_DIR}/file_listing_raw.txt"

# Create a formatted manifest with download URLs
echo "Creating formatted manifest..."
cat > "${OUTPUT_DIR}/MANIFEST.md" << 'EOF'
# Rufus Resources File Manifest

This manifest lists all available files for download.

## Quick Download Commands

### Download a specific file:
```bash
wget https://rufus.marth.lab.s3.us-east-1.amazonaws.com/public_access_data/rufus_resources/<filename>
```

### Download all files:
```bash
aws s3 sync s3://rufus.marth.lab/public_access_data/rufus_resources/ ./rufus_resources/ --no-sign-request
```

---

## Available Files

EOF

# Parse the file listing and create markdown table
echo "| File | Size | Download URL |" >> "${OUTPUT_DIR}/MANIFEST.md"
echo "|------|------|--------------|" >> "${OUTPUT_DIR}/MANIFEST.md"

while read -r line; do
    # Parse: date time size filename
    if [[ $line =~ ^[0-9]{4}-[0-9]{2}-[0-9]{2} ]]; then
        date=$(echo "$line" | awk '{print $1}')
        time=$(echo "$line" | awk '{print $2}')
        size=$(echo "$line" | awk '{print $3}')
        filename=$(echo "$line" | awk '{$1=$2=$3=""; print $0}' | sed 's/^ *//')
        
        # Create relative filename (remove prefix)
        rel_filename="${filename#${PREFIX}}"
        
        # Create download URL
        url="https://${BUCKET}.s3.${REGION}.amazonaws.com/${filename}"
        
        # Add to markdown table
        echo "| \`${rel_filename}\` | ${size} | [Download](${url}) |" >> "${OUTPUT_DIR}/MANIFEST.md"
    fi
done < "${OUTPUT_DIR}/file_listing_human.txt"

# Create a simple text file list for wget
echo "Creating wget download list..."
awk '{print $4}' "${OUTPUT_DIR}/file_listing_raw.txt" | while read -r file; do
    echo "https://${BUCKET}.s3.${REGION}.amazonaws.com/${file}"
done > "${OUTPUT_DIR}/download_urls.txt"

# Create wget script
cat > "${OUTPUT_DIR}/download_all_wget.sh" << 'WGET_EOF'
#!/bin/bash
# Download all Rufus resources using wget

set -e

DEST_DIR="${1:-./rufus_resources}"
mkdir -p "$DEST_DIR"

echo "Downloading Rufus resources to $DEST_DIR..."
echo "This will download approximately 50GB of data."
read -p "Continue? (y/n) " -n 1 -r
echo
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    exit 1
fi

wget -P "$DEST_DIR" -i download_urls.txt -nc -c

echo "Download complete!"
WGET_EOF

chmod +x "${OUTPUT_DIR}/download_all_wget.sh"

# Generate checksums (optional, can be slow for large files)
read -p "Generate MD5 checksums? This may take a while for large files. (y/n) " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    echo "Generating MD5 checksums..."
    cat > "${OUTPUT_DIR}/checksums.md5" << 'CHECKSUM_HEADER'
# MD5 Checksums for Rufus Resources
# Verify downloads with: md5sum -c checksums.md5
CHECKSUM_HEADER
    
    awk '{print $4}' "${OUTPUT_DIR}/file_listing_raw.txt" | while read -r file; do
        echo "Computing checksum for ${file}..."
        checksum=$(aws s3 cp s3://${BUCKET}/"${file}" - --no-sign-request | md5sum | awk '{print $1}')
        rel_file="${file#${PREFIX}}"
        echo "${checksum}  ${rel_file}" >> "${OUTPUT_DIR}/checksums.md5"
    done
    echo "Checksums saved to checksums.md5"
fi

echo ""
echo "Manifest generation complete!"
echo ""
echo "Generated files:"
echo "  - MANIFEST.md           : Human-readable file listing with download links"
echo "  - file_listing_human.txt: File listing with human-readable sizes"
echo "  - file_listing_raw.txt  : File listing with sizes in bytes"
echo "  - download_urls.txt     : Plain text list of download URLs"
echo "  - download_all_wget.sh  : Executable script to download all files with wget"
if [[ $REPLY =~ ^[Yy]$ ]]; then
    echo "  - checksums.md5         : MD5 checksums for verification"
fi
echo ""
echo "You can now commit these files to your repository for users to reference."
