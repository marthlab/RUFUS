#!/bin/bash

#LOCAL_TESTING_UTIL_PATH=/home/ubuntu/RUFUS/singularity/launch_utilities/
#UTIL_PATH=$LOCAL_TESTING_UTIL_PATH
: "${RUFUS_ROOT:=/opt/RUFUS}"

UTIL_PATH=${RUFUS_ROOT}/singularity/launch_utilities/
CHUNK_UTILITIES=${UTIL_PATH}chunk_utilities.sh
. $CHUNK_UTILITIES

KG1_FILE_PATH=${RUFUS_ROOT}/resources/1kg_window_hashes/

sub_dir=$(get_1kg_file "$1" "$2" "$3")
echo "${KG1_FILE_PATH}${sub_dir}"