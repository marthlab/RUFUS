#!/bin/bash

#LOCAL_TESTING_UTIL_PATH=/home/ubuntu/RUFUS/singularity/launch_utilities/
UTIL_PATH=/opt/RUFUS/singularity/launch_utilities/
CHUNK_UTILITIES=${UTIL_PATH}chunk_utilities.sh
. $CHUNK_UTILITIES

region=$(get_chunk_region "$1" "$2" "$3")
echo "$region"
