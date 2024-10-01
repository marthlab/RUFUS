#!/bin/bash
set -e

#LOCAL_TESTING_UTIL_PATH=/home/ubuntu/RUFUS/singularity/launch_utilities/
#UTIL_PATH=$LOCAL_TESTING_UTIL_PATH
UTIL_PATH=/opt/RUFUS/singularity/launch_utilities/

GENOME_HELPERS_PATH=${UTIL_PATH}genome_helpers.sh
. $GENOME_HELPERS_PATH

# Returns string formatted chr:start-end for the given chunk number argument
function get_chunk_region() {
  local chunkNum=$1
  local chunkSize=$2
  local build=$3

  # Convert chunk size to bp (input as kbp)
  local adjusted_size=$((chunkSize * 1000))

  # Get chroms and lengths according to build
  local -a chrom_arr
  get_chroms chrom_arr "$build"
  local -a chrom_lengths
  get_lengths "$build" chrom_lengths

  # Create chunk table corresponding to lengths and chunk size
  local -a chunk_table
  get_chunk_table "$adjusted_size" "$build" chunk_table

  # Find the chromosome that the chunk is in
  chr_idx=0
  curr_chunk=${chunk_table[chr_idx]}
  leftover=$((chunkNum - curr_chunk))
  while [ "$leftover" -ge 0 ]; do
    chr_idx=$((chr_idx + 1))
    curr_chunk=${chunk_table[chr_idx]}
    leftover=$((leftover - curr_chunk))
  done
  chr="chr${chrom_arr[chr_idx]}"

#  # Get the leftover back to the positive range
  leftover=$((leftover + curr_chunk))

  # Calculate the start and end of the chunk
  local chunkStart=$((leftover * adjusted_size + 1))
  local chunkEnd=$((chunkStart + adjusted_size - 1))
  curr_length=${chrom_lengths[chr_idx]}
  if [ "$chunkEnd" -gt "$curr_length" ]; then
    chunkEnd=${chrom_lengths[chr_idx]}
  fi

  # Write out
  echo "$chr:${chunkStart}-${chunkEnd}"
}

# Returns the number of chunks for the given genome build
# Takes in 1) the chunk size and 2) the genome build
function get_num_chunks() {
  local chunk_size=$1
  if [ "$chunk_size" = "0" ]; then
	exit 0
  fi
 
  local adjusted_size=$((chunk_size * 1000))
  local build=$2
  local -a chrom_lengths

  get_lengths "$build" chrom_lengths

  local num_chunks=0
  for ((i=0; i<${#chrom_lengths[@]}; i++)); do

    len=${chrom_lengths[i]}
    curr_chunks=$(($len / $adjusted_size))

    # Add on last chunk if there is a remainder
    remainder=$(($len % $adjusted_size))
    if ((remainder > 0)); then
      curr_chunks=$(($curr_chunks + 1))
    fi

    num_chunks=$((num_chunks + curr_chunks))
  done
  echo "$num_chunks"
}
