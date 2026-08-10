#!/bin/bash
# Shared hash resolution library for RUFUS.
# Container-runtime agnostic — works inside both Docker and Singularity containers.
# The calling workflow is responsible for getting into the container first.
#
# Dependencies (must be sourced/available before this script):
#   - chunk_utilities.sh (get_chunk_region, get_num_chunks)
#   - resource_helpers/download_hash.sh (for S3 downloads)

: "${RUFUS_ROOT:=/opt/RUFUS}"

# Source chunk utilities if not already loaded
if ! declare -f get_chunk_region &>/dev/null; then
    . "${RUFUS_ROOT}/singularity/launch_utilities/chunk_utilities.sh"
fi

# Resolve a single hash file for a given region from a directory.
# Uses glob matching: *{fmtd_region}*.Jhash
#
# Args:
#   $1 - hash_dir: directory containing Jhash files
#   $2 - fmtd_region: formatted region string (e.g., chr1_1_1000000) or "wg" for whole-genome
#
# Returns: prints the resolved file path to stdout
# Exits 1 on error (0 matches or >1 matches)
resolve_hash_for_region() {
    local hash_dir="$1"
    local fmtd_region="$2"

    local -a matches
    mapfile -t matches < <(find "$hash_dir" -maxdepth 1 -name "*${fmtd_region}*.Jhash" \( -type f -o -type l \) 2>/dev/null)

    if [ ${#matches[@]} -eq 0 ]; then
        echo "ERROR: No hash file found matching '*${fmtd_region}*.Jhash' in ${hash_dir}" >&2
        return 1
    elif [ ${#matches[@]} -gt 1 ]; then
        echo "ERROR: Multiple hash files found matching '*${fmtd_region}*.Jhash' in ${hash_dir}:" >&2
        printf '  %s\n' "${matches[@]}" >&2
        echo "Please ensure only one matching .Jhash file exists per region." >&2
        return 1
    fi

    echo "${matches[0]}"
}

# Download all region hashes from S3 for a given hash type and version.
#
# Args:
#   $1 - hash_type: e.g., "kg1", "control"
#   $2 - hash_version: e.g., "v3.0", "v1.0"
#   $3 - window_size: window size in KB (0 for whole-genome)
#   $4 - genome_build: e.g., "GRCh38"
#   $5 - dest_dir: directory to save downloaded hashes
download_hashes() {
    local hash_type="$1"
    local hash_version="$2"
    local window_size="$3"
    local genome_build="$4"
    local dest_dir="$5"

    mkdir -p "$dest_dir"

    local original_dir
    original_dir="$(pwd)"

    cd "$dest_dir" || { echo "ERROR: cannot cd to $dest_dir" >&2; return 1; }

    if [ "$window_size" -eq 0 ]; then
        echo "INFO: Downloading whole-genome ${hash_type} hash (version ${hash_version}) from S3..."
        bash "${RUFUS_ROOT}/resource_helpers/download_hash.sh" "$hash_type" "$hash_version" "wg" \
            || { echo "ERROR: Failed to download whole-genome ${hash_type} hash from S3" >&2; cd "$original_dir"; return 1; }
    else
        local num_chunks
        num_chunks=$(get_num_chunks "$window_size" "$genome_build")
        echo "INFO: Downloading ${num_chunks} ${hash_type} region hashes (version ${hash_version}) from S3..."

        local i fmtd_region region
        for ((i = 0; i < num_chunks; i++)); do
            region=$(get_chunk_region "$i" "$window_size" "$genome_build")
            fmtd_region=$(echo "$region" | tr ':-' '_')

            bash "${RUFUS_ROOT}/resource_helpers/download_hash.sh" "$hash_type" "$hash_version" "$fmtd_region" \
                || { echo "ERROR: Failed to download ${hash_type} hash for region ${region} from S3" >&2; cd "$original_dir"; return 1; }

            # Progress indicator every 100 regions
            if (( (i + 1) % 100 == 0 )); then
                echo "INFO: Downloaded ${hash_type} hashes for $((i + 1))/${num_chunks} regions"
            fi
        done
        echo "INFO: Finished downloading all ${num_chunks} ${hash_type} region hashes"
    fi

    cd "$original_dir"
}

# Validate that hash files exist for all regions in a directory.
#
# Args:
#   $1 - hash_dir: directory containing Jhash files
#   $2 - window_size: window size in KB (0 for whole-genome)
#   $3 - genome_build: e.g., "GRCh38"
#
# Returns 0 if all regions have exactly one matching hash, 1 otherwise.
validate_all_region_hashes() {
    local hash_dir="$1"
    local window_size="$2"
    local genome_build="$3"

    if [ ! -d "$hash_dir" ]; then
        echo "ERROR: Hash directory does not exist: ${hash_dir}" >&2
        return 1
    fi

    if [ "$window_size" -eq 0 ]; then
        # Whole-genome mode: look for *wg*.Jhash
        resolve_hash_for_region "$hash_dir" "wg" > /dev/null \
            || { echo "ERROR: Whole-genome hash validation failed in ${hash_dir}" >&2; return 1; }
        echo "INFO: Validated whole-genome hash in ${hash_dir}"
        return 0
    fi

    local num_chunks
    num_chunks=$(get_num_chunks "$window_size" "$genome_build")

    # Sample first 10 and last 10 indices (deduped) rather than checking all chunks.
    local sample_size=10
    local -a indices=()
    for ((i = 0; i < sample_size && i < num_chunks; i++)); do
        indices+=("$i")
    done
    for ((i = num_chunks - sample_size; i < num_chunks; i++)); do
        if (( i >= sample_size )); then   # avoid duplicates with the first window
            indices+=("$i")
        fi
    done

    echo "INFO: Spot-checking ${#indices[@]} of ${num_chunks} region hashes in ${hash_dir}..."

    local i region fmtd_region errors=0
    for i in "${indices[@]}"; do
        region=$(get_chunk_region "$i" "$window_size" "$genome_build")
        fmtd_region=$(echo "$region" | tr ':-' '_')

        if ! resolve_hash_for_region "$hash_dir" "$fmtd_region" > /dev/null; then
            errors=$((errors + 1))
            if [ $errors -ge 5 ]; then
                echo "ERROR: Too many missing hashes (showed first 5). Aborting validation." >&2
                return 1
            fi
        fi
    done

    if [ $errors -gt 0 ]; then
        echo "ERROR: ${errors} region hash(es) missing or ambiguous in ${hash_dir}" >&2
        return 1
    fi

    echo "INFO: Spot-check passed (${#indices[@]}/${num_chunks} regions) in ${hash_dir}"
    return 0
}
