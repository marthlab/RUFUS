#!/bin/bash
# This script creates three files in the directory it's run in (pwd):
# 1. A rufus call slurm script
# 2. A rufus post-process slurm script
# 3. A bash script to batch submit the two above slurm scripts
# v1.0.0-gamma

#LOCAL_TESTING_UTIL_PATH=/home/ubuntu/RUFUS/singularity/launch_utilities/
#UTIL_PATH=$LOCAL_TESTING_UTIL_PATH
UTIL_PATH=/opt/RUFUS/singularity/launch_utilities/

PARSER=${UTIL_PATH}arg_parser.sh
. $PARSER "$@"

GENOME_HELPERS=${UTIL_PATH}genome_helpers.sh
. $GENOME_HELPERS

CHUNK_UTILITIES=${UTIL_PATH}chunk_utilities.sh
. $CHUNK_UTILITIES

NUM_CHUNKS=$(get_num_chunks "$WINDOW_SIZE_RUFUS_ARG" "$GENOME_BUILD_RUFUS_ARG")

WORKING_DIR=$(pwd)

echo -en "##RUFUS_callCommand=" > rufus.cmd

# Compose run script(s)
RUFUS_SLURM_SCRIPT="rufus_call.slurm"
HEADER_LINES=("#!/bin/bash"
"#SBATCH --job-name=rufus_call"
"#SBATCH --time=${SLURM_TIME_LIMIT_RUFUS_ARG}" 
"#SBATCH --account=${SLURM_ACCOUNT_RUFUS_ARG}" 
"#SBATCH --partition=${SLURM_PARTITION_RUFUS_ARG}"
"#SBATCH --cpus-per-task=${THREAD_LIMIT_RUFUS_ARG}"
"#SBATCH -o ${WORKING_DIR}/slurm_out/%A_%a.out"
"#SBATCH -e ${WORKING_DIR}/slurm_err/%A_%a.err"
)

# Helper function to avoid redundant echoes
function write_out_rest_of_rufus_args() {
    for control in "${CONTROLS_RUFUS_ARG[@]}"; do
            echo -en "-c $control " >> $RUFUS_SLURM_SCRIPT
            echo -en "-c $control " >> rufus.cmd
    done

    # Add in optional hashes if provided
    if [ -n "$REFERENCE_HASH_RUFUS_ARG" ]; then
      echo -en "-f $REFERENCE_HASH_RUFUS_ARG " >> $RUFUS_SLURM_SCRIPT
      echo -en "-f $REFERENCE_HASH_RUFUS_ARG " >> rufus.cmd
    fi
    if [ -n "$EXCLUDE_HASH_LIST_RUFUS_ARG" ]; then
      for exclude in "${EXCLUDE_HASH_LIST_RUFUS_ARG[@]}"; do
        echo -en "-e $exclude " >> $RUFUS_SLURM_SCRIPT
        echo -en "-e $exclude " >> rufus.cmd
      done
    fi
    echo -e "-r $REFERENCE_RUFUS_ARG -m $KMER_DEPTH_CUTOFF_RUFUS_ARG -k 25 -t $THREAD_LIMIT_RUFUS_ARG -L -vs \$REGION_ARG" >> $RUFUS_SLURM_SCRIPT
    echo -e "-r $REFERENCE_RUFUS_ARG -m $KMER_DEPTH_CUTOFF_RUFUS_ARG -k 25 -t $THREAD_LIMIT_RUFUS_ARG -L -vs \$REGION_ARG" >> rufus.cmd
}

# Don't overwrite a run if already exists
if [ -f "$RUFUS_SLURM_SCRIPT" ]; then
	echo "ERROR: $RUFUS_SLURM_SCRIPT already exists - are you overwriting an existing output? Please delete run_rufus.slurm and retry"
	exit 1
fi

for line in "${HEADER_LINES[@]}"
do
    echo -e "$line" >> $RUFUS_SLURM_SCRIPT
done

if [ -n "$EMAIL_RUFUS_ARG" ]; then
	echo -e "#SBATCH --mail-type=ALL" >> $RUFUS_SLURM_SCRIPT
	echo -e "#SBATCH --mail-user=${EMAIL_RUFUS_ARG}" >> $RUFUS_SLURM_SCRIPT
fi

if [ "$WINDOW_SIZE_RUFUS_ARG" = "0" ]; then
	echo "" >> $RUFUS_SLURM_SCRIPT
	echo -e "REGION_ARG=\"\"" >> $RUFUS_SLURM_SCRIPT
else
  # Add a chunk for post-processing
  if [ "$SLURM_ARRAY_JOB_LIMIT_RUFUS_ARG" -lt $((NUM_CHUNKS + 1)) ]; then
    # Always have to subtract one to allow post-processing script queue and shift to 0-index
    ADJ_SLURM_ARRAY_LIMIT=$((SLURM_ARRAY_JOB_LIMIT_RUFUS_ARG - 1))
    # 999

    # Get number of jobs most of the scripts will run (1-based count)
    BASE_COUNT_PER_SCRIPT=$((NUM_CHUNKS / ADJ_SLURM_ARRAY_LIMIT))
    # 3

    # Get remainder that needs to be distributed amongst the last N scripts (0-based count)
    NUM_JOBS_PLUS_ONE=$((NUM_CHUNKS % ADJ_SLURM_ARRAY_LIMIT))
    echo "NUM_JOBS_PLUS_ONE: $NUM_JOBS_PLUS_ONE"
    # 3102 % 999 = 105

    # Get the switch point (i.e. the 0-based array index number where we need to have +1 on the base count)
    NUM_JOBS_BASE_COUNT=$((ADJ_SLURM_ARRAY_LIMIT - NUM_JOBS_PLUS_ONE))
    SWITCH_INDEX=$NUM_JOBS_BASE_COUNT
    # 999 - 105 = 894

    # Sanity check
    RUFUS_CALLS_BASE_COUNT=$((NUM_JOBS_BASE_COUNT * BASE_COUNT_PER_SCRIPT))
    RUFUS_CALLS_PLUS_ONE=$((NUM_JOBS_PLUS_ONE * (BASE_COUNT_PER_SCRIPT + 1)))

    TOTAL_RUFUS_CALLS=$((RUFUS_CALLS_BASE_COUNT + RUFUS_CALLS_PLUS_ONE))
    TOTAL_JOBS=$((NUM_JOBS_BASE_COUNT + NUM_JOBS_PLUS_ONE))

    if [ "$TOTAL_RUFUS_CALLS" != "$NUM_CHUNKS" ] || [ "$TOTAL_JOBS" != "$ADJ_SLURM_ARRAY_LIMIT" ]; then
      echo -e "INFO: $NUM_JOBS_BASE_COUNT slurm array jobs will be run with $BASE_COUNT_PER_SCRIPT rufus calls per script"
      echo -e "INFO: $NUM_JOBS_PLUS_ONE slurm array jobs will be run with $((BASE_COUNT_PER_SCRIPT + 1)) rufus calls per script"
      echo "ERROR: Calculation error in determining number of jobs per script; could not create SLURM scripts"
      exit 1
    else
      echo -e "INFO: $NUM_JOBS_BASE_COUNT slurm array jobs will be run with $BASE_COUNT_PER_SCRIPT rufus calls per script"
      echo -e "INFO: $NUM_JOBS_PLUS_ONE slurm array jobs will be run with $((BASE_COUNT_PER_SCRIPT + 1)) rufus calls per script"
      echo -e "INFO: 1 slurm array job will be run to combine results"
      echo -e "INFO: to fit into the allotted $SLURM_ARRAY_JOB_LIMIT_RUFUS_ARG jobs"
    fi

    # Write out the slurm header
    ADJ_SLURM_ARRAY_END=$((ADJ_SLURM_ARRAY_LIMIT - 1))
    echo -e "#SBATCH -a 0-${ADJ_SLURM_ARRAY_END}%${SLURM_JOB_LIMIT_RUFUS_ARG}" >> $RUFUS_SLURM_SCRIPT
    echo "" >> $RUFUS_SLURM_SCRIPT

    # Write out the region argument and srun command
    echo -e "job_count=$BASE_COUNT_PER_SCRIPT" >> $RUFUS_SLURM_SCRIPT
    echo -e "starting_index=\$((SLURM_ARRAY_TASK_ID * $BASE_COUNT_PER_SCRIPT))" >> $RUFUS_SLURM_SCRIPT
    echo -e "if [ \$SLURM_ARRAY_TASK_ID -eq $SWITCH_INDEX ]; then" >> $RUFUS_SLURM_SCRIPT
    echo -e "   job_count=\$((\$job_count + 1))" >> $RUFUS_SLURM_SCRIPT
    echo -e "elif [ \$SLURM_ARRAY_TASK_ID -gt $SWITCH_INDEX ]; then" >> $RUFUS_SLURM_SCRIPT
    echo -e "   job_count=\$((\$job_count + 1))" >> $RUFUS_SLURM_SCRIPT
    echo -e "   num_jobs_plus_one=\$((\$SLURM_ARRAY_TASK_ID - $SWITCH_INDEX))" >> $RUFUS_SLURM_SCRIPT
    echo -e "   starting_index=\$(($SWITCH_INDEX * $BASE_COUNT_PER_SCRIPT + (\$num_jobs_plus_one * ($BASE_COUNT_PER_SCRIPT + 1))))" >> $RUFUS_SLURM_SCRIPT
    echo -e "fi" >> $RUFUS_SLURM_SCRIPT
    echo -e "for i in \$(seq 0 \$((\$job_count - 1))); do" >> $RUFUS_SLURM_SCRIPT
    echo -e "    curr_job=\$((\$starting_index + \$i))" >> $RUFUS_SLURM_SCRIPT
    echo -e "    region_arg=\$(singularity exec ${CONTAINER_PATH_RUFUS_ARG} bash /opt/RUFUS/singularity/launch_utilities/get_region.sh \"\$curr_job\" \"$WINDOW_SIZE_RUFUS_ARG\" \"$GENOME_BUILD_RUFUS_ARG\")" >> $RUFUS_SLURM_SCRIPT
    echo -e "    REGION_ARG=\"-R \$region_arg\"" >> $RUFUS_SLURM_SCRIPT
    echo -en "   srun --mem=0 singularity exec --bind ${HOST_DATA_DIR_RUFUS_ARG}:/mnt ${CONTAINER_PATH_RUFUS_ARG} bash /opt/RUFUS/runRufus.sh -s /mnt/$SUBJECT_RUFUS_ARG " >> $RUFUS_SLURM_SCRIPT
    echo -en "srun --mem=0 singularity exec --bind ${HOST_DATA_DIR_RUFUS_ARG}:/mnt ${CONTAINER_PATH_RUFUS_ARG} bash /opt/RUFUS/runRufus.sh -s /mnt/$SUBJECT_RUFUS_ARG " >> rufus.cmd
	  echo -en "-pa \$SLURM_ARRAY_TASK_ID " >> $RUFUS_SLURM_SCRIPT
	  echo -en "-cn \$curr_job " >> $RUFUS_SLURM_SCRIPT
    write_out_rest_of_rufus_args
    echo -e "done" >> $RUFUS_SLURM_SCRIPT
  else
      # Write out the slurm header
      ADJ_CHUNK_END=$((NUM_CHUNKS - 1))
      echo -e "#SBATCH -a 0-${ADJ_CHUNK_END}%${SLURM_JOB_LIMIT_RUFUS_ARG}" >> $RUFUS_SLURM_SCRIPT
      echo "" >> $RUFUS_SLURM_SCRIPT

      # Write out the region argument and srun command
      echo -e "region_arg=\$(singularity exec ${CONTAINER_PATH_RUFUS_ARG} bash /opt/RUFUS/singularity/launch_utilities/get_region.sh \"\$SLURM_ARRAY_TASK_ID\" \"$WINDOW_SIZE_RUFUS_ARG\" \"$GENOME_BUILD_RUFUS_ARG\")" >> $RUFUS_SLURM_SCRIPT
      echo -e "REGION_ARG=\"-R \$region_arg\"" >> $RUFUS_SLURM_SCRIPT
      echo -en "srun --mem=0 singularity exec --bind ${HOST_DATA_DIR_RUFUS_ARG}:/mnt ${CONTAINER_PATH_RUFUS_ARG} bash /opt/RUFUS/runRufus.sh -s /mnt/$SUBJECT_RUFUS_ARG " >> $RUFUS_SLURM_SCRIPT
      echo -en "srun --mem=0 singularity exec --bind ${HOST_DATA_DIR_RUFUS_ARG}:/mnt ${CONTAINER_PATH_RUFUS_ARG} bash /opt/RUFUS/runRufus.sh -s /mnt/$SUBJECT_RUFUS_ARG " >> rufus.cmd
      write_out_rest_of_rufus_args
  fi
fi

# Compose post-process slurm wrapper
PP_SLURM_SCRIPT="rufus_post_process.slurm" # Slurm wrapper for post process script
PP_HEADER_LINES=("#!/bin/bash"
"#SBATCH --job-name=rufus_post_process"   
"#SBATCH --account=${SLURM_ACCOUNT_RUFUS_ARG}" 
"#SBATCH --partition=${SLURM_PARTITION_RUFUS_ARG}"
"#SBATCH --output=${WORKING_DIR}/slurm_out/rufus_post_process_%j.out"   
"#SBATCH --error=${WORKING_DIR}/slurm_err/rufus_post_process_%j.err" 
"#SBATCH --nodes=1"
)

for line in "${PP_HEADER_LINES[@]}"
do
    echo -e "$line" >> $PP_SLURM_SCRIPT
done

if [ -n "$EMAIL_RUFUS_ARG" ]; then
    echo -e "#SBATCH --mail-type=ALL" >> $PP_SLURM_SCRIPT
    echo -e "#SBATCH --mail-user=${EMAIL_RUFUS_ARG}" >> $PP_SLURM_SCRIPT
fi
echo "" >> $PP_SLURM_SCRIPT

IFS=$','
CONTROL_STRING="${CONTROLS_RUFUS_ARG[*]}"

BOUND_DATA_DIR="/mnt"

echo -e "srun singularity exec --bind ${HOST_DATA_DIR_RUFUS_ARG}:/mnt ${CONTAINER_PATH_RUFUS_ARG} bash /opt/RUFUS/post_process/post_process.sh -w $WINDOW_SIZE_RUFUS_ARG -r $REFERENCE_RUFUS_ARG -c $CONTROL_STRING -s $SUBJECT_RUFUS_ARG -d ${BOUND_DATA_DIR}" >> $PP_SLURM_SCRIPT

echo -en "##RUFUS_postProcessCommand=" >> rufus.cmd
echo -e "srun singularity exec --bind ${HOST_DATA_DIR_RUFUS_ARG}:/mnt ${CONTAINER_PATH_RUFUS_ARG} bash /opt/RUFUS/post_process/post_process.sh -w $WINDOW_SIZE_RUFUS_ARG -r $REFERENCE_RUFUS_ARG -c $CONTROL_STRING -s $SUBJECT_RUFUS_ARG -d ${BOUND_DATA_DIR}" >> rufus.cmd
mv rufus.cmd "${HOST_DATA_DIR_RUFUS_ARG}"

# Compose invocation script to be executed outside of container
EXE_SCRIPT=launch_rufus.sh
echo -e "#!/bin/bash" > $EXE_SCRIPT
echo -e "" >> $EXE_SCRIPT
echo -e "# This script should be executed after calling the container setup_slurm.sh helper. It requires $PP_SLURM_SCRIPT and $RUFUS_SLURM_SCRIPT to be present in the same directory." >> $EXE_SCRIPT 
echo -e "# Insert command for your system to load singularity here if needed (e.g. module load singularity)" >> $EXE_SCRIPT
echo "" >> $EXE_SCRIPT
echo -e "# Launch calling job" >> $EXE_SCRIPT
echo -e "ARRAY_JOB_ID=\$(sbatch --parsable $RUFUS_SLURM_SCRIPT)" >> $EXE_SCRIPT
echo -e "" >> $EXE_SCRIPT
echo -e "# Launch post-process job - will wait on calling phase to complete" >> $EXE_SCRIPT
echo -e "sbatch --depend=afterany:\$ARRAY_JOB_ID $PP_SLURM_SCRIPT" >> $EXE_SCRIPT

echo -e "Slurm scripts ready to execute with $EXE_SCRIPT. Please make sure singularity is available in your environment, and then run... "
echo -e "bash $EXE_SCRIPT"
