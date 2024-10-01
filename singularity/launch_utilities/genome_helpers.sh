#!/bin/bash

# Returns the chromosomes associated with the given genome build
function get_chroms() {
  local build=$2
  case $build in
    "GRCh38")
      local -n chroms=$1
      chroms=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "X" "Y")
      ;;
    *)
      echo "Genome build not yet supported"
      ;;
  esac
}

# Returns the lengths of the chromosomes associated with the given genome build
function get_lengths() {
  local build=$1
  local -n lengths=$2

  case "$build" in
    "GRCh38")
      lengths=(248956422 242193529 198295559 190214555 181538259 170805979 159345973 145138636 138394717 133797422 135086622 133275309 114364328 107043718 101991189 90338345 83257441 80373285 58617616 64444167 46709983 50818468 156040895 57227415)
      ;;
    *)
      echo "Genome $build not yet supported"
      ;;
  esac
}

# Returns the number of chunks for each chromosome given the chunk size
function get_chunk_table() {
  local chunk_size=$1
  local build=$2
  local -n local_chunk_table=$3

  # Get chrom lengths
  local -a local_lengths
  get_lengths "$build" local_lengths

  for length in "${local_lengths[@]}"; do
    local num_chunks=$((length / chunk_size))
    if ((length % chunk_size != 0)); then
      num_chunks=$((num_chunks + 1))
    fi
    local_chunk_table+=($num_chunks)
  done
}

function get_ref_path() {
	local build=$1

	case "$build" in
		"GRCh38")
			# TODO: need to actually put reference in this spot
			echo "/opt/RUFUS/resources/references/GRCh38_full_analysis_set_plus_decoy_hla.fa"
			;;
		*)
			echo "Genome $build not yet supported"
			;;
	esac
}
