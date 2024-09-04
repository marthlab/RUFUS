#!/bin/bash

#LOCAL_TESTING_UTIL_PATH=/home/ubuntu/RUFUS/singularity/launch_utilities/
#UTIL_PATH=$LOCAL_TESTING_UTIL_PATH
UTIL_PATH=/opt/RUFUS/singularity/launch_utilities/

PARSER=${UTIL_PATH}arg_parser.sh
. $PARSER "$@"

GENOME_HELPERS=${UTIL_PATH}genome_helpers.sh
. $GENOME_HELPERS

CHUNK_UTILITIES=${UTIL_PATH}chunk_utilities.sh
. $CHUNK_UTILITIES

NUM_CHUNKS=$(get_num_chunks $WINDOW_SIZE_RUFUS_ARG $GENOME_BUILD_RUFUS_ARG)

WORKING_DIR=$(pwd)
mkdir -p slurm_out

if [ ! "$WINDOW_SIZE_RUFUS_ARG" = "0" ]; then
	mkdir -p slurm_out/individual_jobs
fi

# TODO: if ARRAY_SIZE_LIMIT < NUM_CHUNKS need to split these into multiple batch scripts

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

# Don't overwrite a run if already exists
if [ -f "$RUFUS_SLURM_SCRIPT" ]; then
	echo "ERROR: run_rufus.slurm already exists - are you overwriting an existing output? Please delete run_rufus.slurm and retry"
	exit 1
fi

for line in "${HEADER_LINES[@]}"
do
    echo -e "$line" >> $RUFUS_SLURM_SCRIPT
done

if [ ! -z $EMAIL_RUFUS_ARG ]; then
	echo -e "#SBATCH --mail-type=ALL" >> $RUFUS_SLURM_SCRIPT
	echo -e "#SBATCH --mail-user=${EMAIL_RUFUS_ARG}" >> $RUFUS_SLURM_SCRIPT
fi

if [ "$WINDOW_SIZE_RUFUS_ARG" = "0" ]; then
	echo "" >> $RUFUS_SLURM_SCRIPT
	echo -e "REGION_ARG=\"\"" >> $RUFUS_SLURM_SCRIPT
else
	echo -e "#SBATCH -a 0-${NUM_CHUNKS}%${SLURM_JOB_LIMIT_RUFUS_ARG}" >> $RUFUS_SLURM_SCRIPT
	echo "" >> $RUFUS_SLURM_SCRIPT
	echo -e "region_arg=\$(singularity exec ${CONTAINER_PATH_RUFUS_ARG}rufus.sif bash /opt/RUFUS/singularity/launch_utilities/get_region.sh \"\$SLURM_ARRAY_TASK_ID\" \"$WINDOW_SIZE_RUFUS_ARG\" \"$GENOME_BUILD_RUFUS_ARG\")" >> $RUFUS_SLURM_SCRIPT
	echo -e "REGION_ARG=\"-R \$region_arg\"" >> $RUFUS_SLURM_SCRIPT
fi

echo -en "srun --mem=0 singularity exec --bind ${HOST_DATA_DIR_RUFUS_ARG}:/mnt ${CONTAINER_PATH_RUFUS_ARG}rufus.sif bash /opt/RUFUS/runRufus.sh -s /mnt/$SUBJECT_RUFUS_ARG " >> $RUFUS_SLURM_SCRIPT
for control in "${CONTROLS_RUFUS_ARG[@]}"; do
	echo -en "-c $control " >> $RUFUS_SLURM_SCRIPT
done

if [ ! -z "$REF_HASH_RUFUS_ARG" ]; then
	echo -en "-f $REF_HASH_RUFUS_ARG " >> $RUFUS_SLURM_SCRIPT
fi
echo -e "-r $REFERENCE_RUFUS_ARG -m $KMER_DEPTH_CUTOFF_RUFUS_ARG -k 25 -t $THREAD_LIMIT_RUFUS_ARG -L -vs \$REGION_ARG" >> $RUFUS_SLURM_SCRIPT

# Compose post-process slurm script
PP_SLURM_SCRIPT="postProcessRufus.slurm"
PP_HEADER_LINES=("#!/bin/bash"
"#SBATCH --job-name=rufus_post_process"   
"#SBATCH --output=rufus_post_process_%j.out"   
"#SBATCH --error=rufus_post_process_%j.err" 
"#SBATCH --nodes=1"
)

IFS=$'\t'
TAB_CONTROLS="${CONTROLS_RUFUS_ARG[*]}"

PP_VAR_LINES=(""
"SUBJECT_FILE=$SUBJECT_RUFUS_ARG"
"WINDOW_SIZE=$WINDOW_SIZE_RUFUS_ARG"
"CONTROL_STRING=$TAB_CONTROLS"
"REFERENCE=$REFERENCE_RUFUS_ARG"
)

for line in "${PP_HEADER_LINES[@]}"
do
    echo -e "$line" >> $PP_SLURM_SCRIPT
done

if [ ! -z $EMAIL_RUFUS_ARG ]; then
    echo -e "#SBATCH --mail-type=ALL" >> $PP_SLURM_SCRIPT
    echo -e "#SBATCH --mail-user=${EMAIL_RUFUS_ARG}" >> $PP_SLURM_SCRIPT
fi

cat "${FILE_STUB_DIR}postProcess.body" >> $PP_SLURM_SCRIPT

# Get rufus run(s) going
# TODO: need to iterate through all batch script args, collect JOBIDs and wait on all
#ARRAY_JOB_ID=$(sbatch --parsable $RUFUS_SLURM_SCRIPT)

#sbatch --depend=afterany:$ARRAY_JOB_ID singularity exec bash ${PP_SLURM_SCRIPT}
echo "All RUFUS runs completed. Beginning post-processing..."
