#!/bin/bash

# TODO: this will need to change to singularity path
GENOME_HELPERS_PATH=/Users/genetics/Documents/code/RUFUS/singularity/launch_utilities/genome_helpers.sh

. $GENOME_HELPERS_PATH

# Returns string formatted chr:start-end for the given chunk number argument
function get_chunk_region() {
  local chunkNum=$1
  local chunkSize=$2
  local build=$3

  local -a chrom_lengths
  get_lengths "$build" chrom_lengths

  local -a chrom_arr
  get_chroms chrom_arr "$build"

  # Calculate chromosome and coordinates for the given chunk
  chunkStart=0
  chunkEnd=0

  # Calculate the absolute coordinate of the chunk
  absCoord=$((chunkNum * chunkSize))

  # Find the chromosome that the chunk is in
  idx=0
  remainder=$((absCoord - chrom_lengths[idx]))
  while ((remainder > 0)); do
    idx=$((idx + 1))
    remainder=$((remainder - chrom_lengths[idx]))
  done
  # Translate idx into chr
  chr="chr${chrom_arr[idx]}"

  # Get remainder back to positive
  remainder=$((remainder + chrom_lengths[idx]))

  # Calculate the start and end of the chunk
  multiplier=$((remainder/chunkSize))
  chunkStart=$((multiplier * chunkSize))
  if [[ $chunkStart -eq 0 ]]; then
    chunkStart=1
  fi

  chunkEnd=$((chunkStart + chunkSize - 1))

  # Adjust the end of the chunk if it exceeds the chromosome length
  if ((chunkEnd > chrom_lengths[idx])); then
    chunkEnd=$((chrom_lengths[idx]))
  fi

  echo "$chr:${chunkStart}-${chunkEnd}"
}

# Returns the number of chunks for the given genome build
# Takes in 1) the chunk size and 2) the genome build
function get_num_chunks() {
  local chunk_size=$1
  local build=$2
  local -a chrom_lengths

  get_lengths "$build" chrom_lengths

  local num_chunks=0
  for ((i=0; i<${#chrom_lengths[@]}; i++)); do

    len=${chrom_lengths[i]}
    curr_chunks=$(($len / $chunk_size))

    # Add on last chunk if there is a remainder
    remainder=$(($len % $chunk_size))
    if ((remainder > 0)); then
      curr_chunks=$(($curr_chunks + 1))
    fi

    num_chunks=$((num_chunks + curr_chunks))
  done
  echo "$num_chunks"
}