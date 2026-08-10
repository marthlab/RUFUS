#!/bin/bash

# Env override
: "${RUFUS_ROOT:=/opt/RUFUS}"

# Arguments
CONTAINER_ID="$1"
ENV_FILE="$2"
REGION="$3"

# Constants
DEFAULT_KG1_HASH_VERSION="v3.0"
DEFAULT_CONTROL_HASH_VERSION="v1.0"

# Import env file now that we're inside of container
set -a
source <(grep -v '^#' $ENV_FILE | grep -v '^[[:space:]]*$' | sed 's/\r$//')
set +a

# Resolves a hash file for a given region and hash type.
# Checks local directory first (via shared resolve_hashes.sh), falls back to S3 download.
#
# Args:
#   $1 - region: region string (e.g., chr1:1-1000000) or empty for whole-genome
#   $2 - geo_type: "local" to check local dir first, "remote" for S3 only
#   $3 - hash_type: e.g., "kg1", "control"
#
# Returns: prints resolved hash file path to stdout
get_hash() {
    local region="$1"
    local geo_type="$2"
    local hash_type="$3"
    local hash_type_upper
    hash_type_upper=$(echo "$hash_type" | tr '[:lower:]' '[:upper:]')

    local fmtd_region
    if [ -z "$region" ]; then
        fmtd_region="wg"
    else
        fmtd_region=$(echo "$region" | tr ':-' '_')
    fi

    # Try local resolution via shared resolve_hashes.sh
    if [ "$geo_type" == "local" ]; then
        local env_var="${hash_type_upper}_HASH_LOCAL_DIR"
        local host_dir="${!env_var}"

        local resolved
        resolved=$(docker exec "$CONTAINER_ID" bash -c \
            "source ${RUFUS_ROOT}/resource_helpers/resolve_hashes.sh && resolve_hash_for_region '${host_dir}' '${fmtd_region}'" 2>/dev/null)

        if [ $? -eq 0 ] && [ -n "$resolved" ]; then
            echo "Using local ${hash_type} hash for ${fmtd_region}: $resolved" >&2
            echo "$resolved"
            return 0
        fi

        echo "Could not find local ${hash_type} hash for ${fmtd_region} in ${host_dir}. Falling back to S3." >&2
    fi

    # S3 download fallback
    local version_var="${hash_type_upper}_HASH_VERSION"
    local hash_version="${!version_var}"
    if [ -z "$hash_version" ]; then
        local default_var="DEFAULT_${hash_type_upper}_HASH_VERSION"
        hash_version="${!default_var}"
    fi

    echo "Fetching version ${hash_version} ${hash_type} hash for ${fmtd_region} from S3" >&2
    docker exec "$CONTAINER_ID" bash -c \
        "cd /home/ubuntu && bash ${RUFUS_ROOT}/resource_helpers/download_hash.sh '${hash_type}' '${hash_version}' '${fmtd_region}'" >&2 \
        || { echo "ERROR: Failed to download ${hash_type} hash for ${fmtd_region} from S3" >&2; return 1; }

    echo "/home/ubuntu/${fmtd_region}_${hash_type}_${hash_version}.Jhash"
}
export -f get_hash


# Compose control argument of both or one of control hashes and paired control files
ctrl_arg=""
if [ "$CONTROL_HASH_LOCAL_DIR" != "" ]; then
    control_hash=$(get_hash $REGION "local" "control") || exit 1
    ctrl_arg+="-e $control_hash "
elif [ "$CONTROL_HASH_VERSION" != "" ]; then
    echo "Using control hash version: $CONTROL_HASH_VERSION" >&2
    control_hash=$(get_hash $REGION "remote" "control") || exit 1
    ctrl_arg+="-e $control_hash "
elif [ "${#CONTROL_FILE_ARRAY[@]}" -eq 0 ]; then
    echo "No local control hashes, paired controls, or control hash version provided, fetching default $DEFAULT_CONTROL_HASH_VERSION hashes piecemeal" >&2
    CONTROL_HASH_VERSION="$DEFAULT_CONTROL_HASH_VERSION"
    control_hash=$(get_hash $REGION "remote" "control") || exit 1
    ctrl_arg+="-e $control_hash "
fi

# If we have paired controls provided, also use those
if [ "${#CONTROL_FILE_ARRAY[@]}" -ne 0 ]; then
    # Concatenate controls into -c delimited string
    for control in "${CONTROL_FILE_ARRAY[@]}"; do
        ctrl_arg+="-c $control "
    done
fi

# Compose kg1 hash argument
if [ "$KG1_HASH_LOCAL_DIR" != "" ]; then
    kg1_hash=$(get_hash $REGION "local" "kg1") || exit 1
    kg1_hash_arg="-e $kg1_hash"
elif [ "$KG1_HASH_VERSION" != "" ]; then
    kg1_hash=$(get_hash $REGION "remote" "kg1") || exit 1
    kg1_hash_arg="-e $kg1_hash"
fi

ref_arg="-r $REFERENCE_FASTA"
# If subject_file ends with cram, need to change region_arg to -cr
if [[ "$SUBJECT_FILE" == *.cram ]]; then
    ref_arg="-cr $REFERENCE_FASTA"
fi

# Region arg
if [ "$REGION" != "" ]; then
    region_arg="-R $REGION"
    fmtd_reg=$(echo "$REGION" | tr ':-' '_')
else
    fmtd_reg="whole_genome"
fi

cd $WORKING_DIR

echo "Running RUFUS for $REGION on $SUBJECT_FILE..."

RUFUS_CMD="$RUFUS_ROOT/runRufus.sh \
  -s $SUBJECT_FILE \
  $ctrl_arg \
  $ref_arg \
  -m $KMER_DEPTH_CUTOFF \
  -k $KMER_LENGTH \
  -t $THREAD_LIMIT \
  $OTHER_FLAGS \
  $kg1_hash_arg \
  $region_arg"

docker exec "$CONTAINER_ID" bash -c \
    "$RUFUS_CMD \
    > rufus_supplementals/logs/${fmtd_reg}.out \
    2> rufus_supplementals/logs/${fmtd_reg}.err"

# Clean up hash files
if [ "$REGION" == "" ]; then
    if [ -f "/home/ubuntu/downloaded_control_hashes/*wg*.Jhash" ]; then
        docker exec ${CONTAINER_ID} rm /home/ubuntu/downloaded_control_hashes/*wg*control*.Jhash
    fi

    if [ -f "/home/ubuntu/downloaded_kg1_hashes/*wg*.Jhash" ]; then
        docker exec ${CONTAINER_ID} rm /home/ubuntu/downloaded_kg1_hashes/*wg*.Jhash
    fi

else 
    if [ -f "/home/ubuntu/downloaded_control_hashes/*$fmtd_reg*.Jhash" ]; then
        docker exec ${CONTAINER_ID} rm /home/ubuntu/downloaded_control_hashes/*$fmtd_reg*control*.Jhash
    fi

    if [ -f "/home/ubuntu/downloaded_kg1_hashes/*$fmtd_reg*.Jhash" ]; then
        docker exec ${CONTAINER_ID} rm /home/ubuntu/downloaded_kg1_hashes/*$fmtd_reg*.Jhash
    fi
fi