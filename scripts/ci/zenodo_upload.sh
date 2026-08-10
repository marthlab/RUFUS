#!/bin/bash
# Publish a new version of the RUFUS Zenodo record with the given SIF as its file.
#
# Usage: zenodo_upload.sh <sif_path> <version_tag>
#
# Creates a new version under an existing concept record (so all releases share one concept
# DOI), uploads the SIF, sets minimal metadata, and publishes. Intended to run in CI on a
# git tag push (see .github/workflows/build-publish.yml).
#
# Required environment:
#   ZENODO_TOKEN              personal access token with deposit:write + deposit:actions
#   ZENODO_CONCEPT_RECORD_ID  the concept (all-versions) record id of the existing RUFUS record
set -euo pipefail

SIF_PATH="${1:?usage: zenodo_upload.sh <sif_path> <version_tag>}"
VERSION="${2:?usage: zenodo_upload.sh <sif_path> <version_tag>}"
: "${ZENODO_TOKEN:?ZENODO_TOKEN must be set}"
: "${ZENODO_CONCEPT_RECORD_ID:?ZENODO_CONCEPT_RECORD_ID must be set}"

API="https://zenodo.org/api"
AUTH="Authorization: Bearer ${ZENODO_TOKEN}"

# jq is preinstalled on GitHub-hosted ubuntu runners.
api() {
    # api <curl-args...> ; echoes the JSON body, fails on HTTP >= 400
    local body http
    body=$(curl -sS -w $'\n%{http_code}' -H "$AUTH" "$@")
    http=$(tail -n1 <<<"$body")
    body=$(sed '$d' <<<"$body")
    if [ "$http" -ge 400 ]; then
        echo "Zenodo API error (HTTP $http):" >&2
        echo "$body" >&2
        return 1
    fi
    echo "$body"
}

echo "Resolving latest deposition for concept record ${ZENODO_CONCEPT_RECORD_ID}..."
LATEST_ID=$(api "${API}/records/${ZENODO_CONCEPT_RECORD_ID}" | jq -r '.id')
echo "Latest record id: ${LATEST_ID}"

echo "Creating new version..."
NEWVER=$(api -X POST "${API}/deposit/depositions/${LATEST_ID}/actions/newversion")
DRAFT_URL=$(jq -r '.links.latest_draft' <<<"$NEWVER")
DRAFT_ID="${DRAFT_URL##*/}"
echo "New draft deposition id: ${DRAFT_ID}"

# Remove files inherited from the previous version so only the new SIF remains.
echo "Clearing inherited files..."
DRAFT=$(api "${API}/deposit/depositions/${DRAFT_ID}")
BUCKET=$(jq -r '.links.bucket' <<<"$DRAFT")
for fid in $(jq -r '.files[].id' <<<"$DRAFT"); do
    api -X DELETE "${API}/deposit/depositions/${DRAFT_ID}/files/${fid}" >/dev/null || true
done

echo "Uploading ${SIF_PATH}..."
FNAME=$(basename "$SIF_PATH")
api -X PUT "${BUCKET}/${FNAME}" --upload-file "$SIF_PATH" >/dev/null

echo "Setting version metadata..."
api -X PUT "${API}/deposit/depositions/${DRAFT_ID}" \
    -H "Content-Type: application/json" \
    -d "{\"metadata\": $(jq -n --arg v "$VERSION" '.version=$v | .publication_date=(now|strftime("%Y-%m-%d"))' \
            <<<"$(jq '.metadata' <<<"$DRAFT")")}" >/dev/null

echo "Publishing..."
PUBLISHED=$(api -X POST "${API}/deposit/depositions/${DRAFT_ID}/actions/publish")
DOI=$(jq -r '.doi' <<<"$PUBLISHED")
echo "Published RUFUS ${VERSION} to Zenodo. DOI: ${DOI}"
