#!/bin/bash
# This is run within container

# ENV override
: "${RUFUS_ROOT:=/opt/RUFUS}"

# Constants
CMD_OUT="./rufus_temp/rufus.cmd"
GLOBALS_FILE="$RUFUS_ROOT/resources/globals.txt"
USER_SPEC="$(id -u):$(id -g)"

# Required args
CONTAINER_ID="$1"
CONT_ENV_FILE="./rufus_temp/temp_rufus.env"

# Make env variables available
if [ -f "$CONT_ENV_FILE" ]; then
    set -a
    source <(grep -v '^#' $CONT_ENV_FILE | grep -v '^[[:space:]]*$' | sed 's/\r$//')
    set +a
else
    echo "Error: $CONT_ENV_FILE file not found - please provide valid path to rufus.env file"
    exit 1
fi

# Make globals available
if [ -f "$GLOBALS_FILE" ]; then
    set -a
    source <(grep -v '^#' $GLOBALS_FILE | grep -v '^[[:space:]]*$' | sed 's/\r$//')
    set +a
else
    echo "Error: $GLOBALS_FILE file not found - please provide valid path to rufus.env file"
    exit 1
fi

# get control argument or internal version
ctrl_arg=""
if [ "${#CONTROL_FILE_ARRAY[@]}" -eq 0 ]; then
    ctrl_arg="-e internal_$CONTROL_HASH_VERSION"
else
    # Concatenate controls into a single -c delimited string
    for control in "${CONTROL_FILE_ARRAY[@]}"; do
        ctrl_arg+="-c $control"
    done
fi

run_cmd="docker -u ${USER_SPEC} exec ${CONTAINER_ID} bash $RUFUS_ROOT/runRufus.sh -s $SUBJECT_FILE $ctrl_arg -r $REFERENCE_FASTA -k $KMER_LENGTH -m $KMER_DEPTH_CUTOFF -t $THREAD_LIMIT $OTHER_FLAGS -e kg1_$KG1_HASH_VERSION"
post_cmd="docker exec -u ${USER_SPEC} ${CONTAINER_ID} bash $RUFUS_ROOT/post_process/post_process.sh -s $SUBJECT_FILE -r $REFERENCE_FASTA -w $WINDOW_SIZE -d . $ctrl_arg"

mkdir -p ./rufus_supplementals

echo "##RUFUSCommandLine=<ID=rufus, Branch=\"$RUFUS_BRANCH\", Version=\"$RUFUS_VERSION\", Command=\"$run_cmd\">" > "$CMD_OUT"
echo "##RUFUSCommandLine=<ID=rufus, Branch=\"$RUFUS_BRANCH\", Version=\"$RUFUS_VERSION\", Command=\"$post_cmd\">" >> "$CMD_OUT"