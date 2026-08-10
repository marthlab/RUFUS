#!/bin/bash
# This script creates three files in the directory it's run in (pwd):
# 1. A rufus call slurm script
# 2. A rufus post-process slurm script
# 3. A bash script to batch submit the two above slurm scripts
# d1.1.9

#LOCAL_TESTING_UTIL_PATH=/home/ubuntu/RUFUS/singularity/launch_utilities/
#UTIL_PATH=$LOCAL_TESTING_UTIL_PATH

: "${RUFUS_ROOT:=/opt/RUFUS}"
UTIL_PATH=${RUFUS_ROOT}/singularity/launch_utilities/

# TODO: detect the container runtime instead of hard-coding `singularity`.
#
# The SLURM scripts generated below emit `singularity exec ...` (10 sites in this file, plus
# get_region.sh invocations). That works on Apptainer hosts only because Apptainer installs a
# `singularity` compatibility shim -- it is not a guarantee, and a site that ships Apptainer
# without the shim cannot run a generated script. Flipping the literal to `apptainer` just moves
# the breakage to sites still on Singularity CE, so neither hard-coded name is right. Detect once
# here and substitute the result into the generated scripts:
#
#   CONTAINER_CMD="$(command -v apptainer || command -v singularity)" \
#       || { echo "ERROR: neither apptainer nor singularity found on PATH" >&2; exit 1; }
#
# Deliberately deferred: this changes every generated SLURM script, which is exactly the machinery
# the sharded whole-genome release gate exercises. Land it after that run, not into it. The docs
# and the newer functional cases (tests/functional/cases/f*.sh) already say/use `apptainer`; this
# file and singularity/tests/ are the remaining holdouts, and those tests are separately stale
# (hard-coded /home/ubuntu paths, references to Child/Mother/Father.bam that do not exist).

PARSER=${UTIL_PATH}arg_parser.sh
. $PARSER "$@"

GENOME_HELPERS=${UTIL_PATH}genome_helpers.sh
. $GENOME_HELPERS

CHUNK_UTILITIES=${UTIL_PATH}chunk_utilities.sh
. $CHUNK_UTILITIES

NUM_CHUNKS=$(get_num_chunks "$WINDOW_SIZE_RUFUS_ARG" "$GENOME_BUILD_RUFUS_ARG")

# Resolve per-region hashes (download from S3 if needed, validate all exist)
RESOLVE_HASHES=${RUFUS_ROOT}/resource_helpers/resolve_hashes.sh
. $RESOLVE_HASHES

if [ -n "$KG1_HASH_VERSION" ]; then
    KG1_HASH_DIR="$(pwd)/rufus_hashes/kg1"
    download_hashes "kg1" "$KG1_HASH_VERSION" "$WINDOW_SIZE_RUFUS_ARG" "$GENOME_BUILD_RUFUS_ARG" "$KG1_HASH_DIR" \
        || { echo "ERROR: Failed to download KG1 hashes"; exit 1; }
fi

if [ -n "$CONTROL_HASH_VERSION" ]; then
    CONTROL_HASH_DIR="$(pwd)/rufus_hashes/control"
    download_hashes "control" "$CONTROL_HASH_VERSION" "$WINDOW_SIZE_RUFUS_ARG" "$GENOME_BUILD_RUFUS_ARG" "$CONTROL_HASH_DIR" \
        || { echo "ERROR: Failed to download control hashes"; exit 1; }
fi

if [ -n "$KG1_HASH_DIR" ]; then
    validate_all_region_hashes "$KG1_HASH_DIR" "$WINDOW_SIZE_RUFUS_ARG" "$GENOME_BUILD_RUFUS_ARG" \
        || { echo "ERROR: KG1 hash validation failed"; exit 1; }
fi

if [ -n "$CONTROL_HASH_DIR" ]; then
    validate_all_region_hashes "$CONTROL_HASH_DIR" "$WINDOW_SIZE_RUFUS_ARG" "$GENOME_BUILD_RUFUS_ARG" \
        || { echo "ERROR: Control hash validation failed"; exit 1; }
fi

# Re-compute bind mounts now that hash dirs may have been set by S3 downloads
BIND_MOUNTS="$(collect_bind_dirs)"

# Preflight every generator input by running it inside the container.
#
# A generator is an arbitrary shell script RUFUS executes to obtain reads, so not all of its
# dependencies are discoverable by reading it: paths assembled from environment variables, the
# REF_PATH/REF_CACHE lookup samtools uses to decode CRAM without an explicit -T, and any binary it
# shells out to. collect_bind_dirs() binds what it can find statically; this catches the rest.
#
# Without it those failures surface inside a queued job -- in windowed mode, across every task in
# the array -- hours after submission. Running the generator here costs seconds and surfaces the
# generator's own error text.
#
# Uses a detected runtime rather than the hard-coded `singularity` of the generated scripts (see
# the TODO at the top of this file): this is new setup-time code, so it can do the right thing
# without touching the generated-script machinery that deferral is about. If no runtime is on PATH
# the check is skipped with a warning -- setup has never required one, and refusing to generate
# scripts on a submit host without a container runtime would be a regression.
preflight_generators() {
    local generators=()
    local f
    for f in "${SUBJECTS_RUFUS_ARG[@]}" "${CONTROLS_RUFUS_ARG[@]}"; do
        [ "$(get_input_type "$f")" == "generator" ] && generators+=("$f")
    done
    [ ${#generators[@]} -eq 0 ] && return 0

    # setup_slurm.sh is normally invoked from *inside* the container:
    #   apptainer exec rufus.sif bash /opt/RUFUS/singularity/setup_slurm.sh ...
    # as every README example shows. No container runtime exists inside the image, and none is
    # needed -- this is already the environment the generator will run in, samtools included -- so
    # run the generator directly. Only enter the container when setup is running on the host.
    #
    # Note the binds differ between the two cases: in-container we inherit whatever the caller's
    # exec bound (at CHPC the singularity/apptainer module sets *_BINDPATH=/scratch,/uufs), whereas
    # the generated job scripts use BIND_MOUNTS. A generator that reads from a path bound only in
    # the job environment can therefore fail here; the error text below says so.
    local -a run_prefix=()
    if [ -n "${APPTAINER_CONTAINER:-}${SINGULARITY_CONTAINER:-}" ] || [ -d /.singularity.d ]; then
        : # already inside the container -- run_prefix stays empty
    else
        local runtime sif
        runtime="$(command -v apptainer || command -v singularity)" || runtime=""
        if [ -z "$runtime" ]; then
            echo "WARNING: setup is not running inside a container and neither apptainer nor" >&2
            echo "         singularity is on PATH; skipping the generator preflight. Generator" >&2
            echo "         errors will not surface until the jobs run." >&2
            return 0
        fi
        sif="${CONTAINER_PATH_RUFUS_ARG:-rufus.sif}"
        if [ ! -f "$sif" ]; then
            echo "WARNING: container $sif not found; skipping the generator preflight." >&2
            return 0
        fi
        run_prefix=("$runtime" exec --bind "${BIND_MOUNTS}${DEV_BIND_ARGS}" "$sif")
    fi

    local gen out rc checked=0
    for gen in "${generators[@]}"; do
        # A generator with a pre-built hash beside it is never executed: runRufus.sh keeps the
        # <generator><region_postfix> naming and RunJellyForRUFUS.sh skips jellyfish when
        # <generator><region_postfix>.Jhash exists. That is how pre-built DSA/control hashes are
        # supplied, and such a generator is legitimately empty, so there is nothing to preflight.
        # Generators are rejected in windowed mode, so .wg is the only postfix reachable here.
        if [ -e "${gen}.wg.Jhash" ]; then
            echo "Skipping preflight for $(basename "$gen"): pre-built hash ${gen}.wg.Jhash will be used instead."
            continue
        fi
        checked=$((checked + 1))

        # head closes the pipe once it has enough to judge, so the generator is not run to
        # completion; its exit status is therefore not meaningful and the output is what we check.
        out="$(timeout 120 "${run_prefix[@]}" bash -c "bash '$gen' 2>&1 | head -20" 2>&1)"
        rc=$?

        if [ $rc -eq 124 ]; then
            echo "ERROR: generator $gen produced no output within 120s inside the container." >&2
            echo "       A generator that blocks this long is usually waiting on a reference or" >&2
            echo "       index that is not bound into the container." >&2
            exit 1
        fi

        if echo "$out" | grep -qE '^@(HD|SQ|RG|PG|CO)[[:space:]]'; then
            continue
        fi
        if echo "$out" | awk -F'\t' 'NF>=11 { found=1; exit } END { exit !found }'; then
            continue
        fi

        echo "ERROR: generator $gen did not produce SAM inside the container." >&2
        echo "       RUFUS runs generators with 'bash <generator>' and expects SAM on stdout." >&2
        echo "       Output was:" >&2
        if [ -n "$out" ]; then
            echo "$out" | sed 's/^/         /' >&2
        else
            echo "         (no output)" >&2
        fi
        echo "       An empty generator, or one whose data paths are not visible here, fails this" >&2
        echo "       way. If the paths are only bound in the job environment, bind them for setup" >&2
        echo "       too (SINGULARITY_BINDPATH/--bind) or add them with -d." >&2
        exit 1
    done
    if [ "$checked" -gt 0 ]; then
        echo "Generator preflight passed ($checked of ${#generators[@]} generator(s) produced SAM in the container)."
    fi
}
preflight_generators

WORKING_DIR=$(pwd)

echo -en "##RUFUS_callCommand=" > rufus.cmd

# Build subject args string for runRufus.sh (each subject gets its own -s flag)
SUBJECT_ARGS_STRING=""
for subj in "${SUBJECTS_RUFUS_ARG[@]}"; do
    SUBJECT_ARGS_STRING+="-s $subj "
done

# Compose run script(s)
RUFUS_SLURM_SCRIPT="rufus_call.slurm"
HEADER_LINES=("#!/bin/bash"
"#SBATCH --job-name=rufus_call"
"#SBATCH --time=${SLURM_TIME_LIMIT_RUFUS_ARG}" 
"#SBATCH --account=${SLURM_ACCOUNT_RUFUS_ARG}" 
"#SBATCH --partition=${SLURM_PARTITION_RUFUS_ARG}"
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

    # Static exclude hashes (same for all regions)
    if [ -n "$EXCLUDE_HASH_LIST_RUFUS_ARG" ]; then
      for exclude in "${EXCLUDE_HASH_LIST_RUFUS_ARG[@]}"; do
        echo -en "-e $exclude " >> $RUFUS_SLURM_SCRIPT
        echo -en "-e $exclude " >> rufus.cmd
      done
    fi

    # Per-region hash args (resolved at SLURM job runtime via glob)
    if [ -n "$KG1_HASH_DIR" ] || [ -n "$CONTROL_HASH_DIR" ]; then
      echo -en "\$HASH_ARGS " >> $RUFUS_SLURM_SCRIPT
      echo -en "\$HASH_ARGS " >> rufus.cmd
    fi

    # -cr is needed if ANY input is a CRAM, subject or control: runRufus.sh kills the run the moment
    # it decodes a .cram with _arg_cramref unset. Keying this off the first subject alone broke every
    # mixed set (e.g. -s x.generator -c y.cram, or a BAM subject with a CRAM control).
    #
    # Both flags are emitted rather than swapping one for the other. runRufus.sh only assigns
    # _arg_ref from _arg_cramref inside its per-file CRAM branches, which run after it computes
    # _arg_ref_cat="${_arg_ref%.*}". With -cr alone and a non-CRAM subject, _arg_ref_cat would be
    # empty at that point and the BWA prefix would fall back to the reference path instead of the
    # extension-stripped prefix. Passing -r as well sets _arg_ref up front; both take the same path.
    local ref_flags="-r $REFERENCE_RUFUS_ARG"
    local _input
    for _input in "${SUBJECTS_RUFUS_ARG[@]}" "${CONTROLS_RUFUS_ARG[@]}"; do
        if [[ "$_input" == *.cram ]]; then
            ref_flags="-r $REFERENCE_RUFUS_ARG -cr $REFERENCE_RUFUS_ARG"
            break
        fi
    done
    echo -en "$ref_flags -m $KMER_DEPTH_CUTOFF_RUFUS_ARG -k 25 -t $THREAD_LIMIT_RUFUS_ARG -L -vs " >> $RUFUS_SLURM_SCRIPT
    echo -en "$ref_flags -m $KMER_DEPTH_CUTOFF_RUFUS_ARG -k 25 -t $THREAD_LIMIT_RUFUS_ARG -L -vs " >> rufus.cmd

    # Hash size (-hs) for the k-mer count step. Must match the -s of any pre-built control/DSA/exclude
    # hash or runRufus.sh's merge aborts; left unset, RUFUS applies its own default (16G wg / 1G window).
    if [ -n "$HASH_SIZE_RUFUS_ARG" ]; then
      echo -en "-hs $HASH_SIZE_RUFUS_ARG " >> $RUFUS_SLURM_SCRIPT
      echo -en "-hs $HASH_SIZE_RUFUS_ARG " >> rufus.cmd
    fi

    if [ "${PAR_LOW_COV_THRESHOLD_RUFUS_ARG}" != "7" ]; then
      echo -en "-plct $PAR_LOW_COV_THRESHOLD_RUFUS_ARG " >> $RUFUS_SLURM_SCRIPT
      echo -en "-plct $PAR_LOW_COV_THRESHOLD_RUFUS_ARG " >> rufus.cmd
    fi

    if [ "$WINDOW_SIZE_RUFUS_ARG" -ne 0 ]; then
      echo -en "\$REGION_ARG " >> $RUFUS_SLURM_SCRIPT
      echo -en "\$REGION_ARG " >> rufus.cmd
    fi
    
    printf '\n' >> "$RUFUS_SLURM_SCRIPT"
    printf '\n' >> rufus.cmd
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

if [ "$WINDOW_SIZE_RUFUS_ARG" -eq 0 ]; then
  echo -e "#SBATCH --mem=${MEM_PER_JOB}" >> $RUFUS_SLURM_SCRIPT
  echo -e "#SBATCH --cpus-per-task=${CPUS_PER_JOB}" >> $RUFUS_SLURM_SCRIPT
  echo -e "#SBATCH -o ${WORKING_DIR}/slurm_out/rufus_call_%j.out" >> $RUFUS_SLURM_SCRIPT
  echo -e "#SBATCH -e ${WORKING_DIR}/slurm_out/rufus_call_%j.err" >> $RUFUS_SLURM_SCRIPT
  printf '\n' >> $RUFUS_SLURM_SCRIPT

  # Resolve whole-genome hash files at setup time (static paths)
  if [ -n "$KG1_HASH_DIR" ] || [ -n "$CONTROL_HASH_DIR" ]; then
    WG_HASH_ARGS=""
    if [ -n "$KG1_HASH_DIR" ]; then
      wg_kg1_hash=$(resolve_hash_for_region "$KG1_HASH_DIR" "wg") \
          || { echo "ERROR: Could not resolve whole-genome KG1 hash"; exit 1; }
      WG_HASH_ARGS="$WG_HASH_ARGS -e $wg_kg1_hash"
    fi
    if [ -n "$CONTROL_HASH_DIR" ]; then
      wg_ctrl_hash=$(resolve_hash_for_region "$CONTROL_HASH_DIR" "wg") \
          || { echo "ERROR: Could not resolve whole-genome control hash"; exit 1; }
      WG_HASH_ARGS="$WG_HASH_ARGS -e $wg_ctrl_hash"
    fi
    echo -e "HASH_ARGS=\"${WG_HASH_ARGS}\"" >> $RUFUS_SLURM_SCRIPT
  fi

	echo -en "srun --mem=${MEM_PER_JOB} singularity exec --bind ${BIND_MOUNTS}${DEV_BIND_ARGS} ${CONTAINER_PATH_RUFUS_ARG} bash /opt/RUFUS/runRufus.sh $SUBJECT_ARGS_STRING" >> $RUFUS_SLURM_SCRIPT
  echo -en "srun --mem=${MEM_PER_JOB} singularity exec --bind ${BIND_MOUNTS}${DEV_BIND_ARGS} ${CONTAINER_PATH_RUFUS_ARG} bash /opt/RUFUS/runRufus.sh $SUBJECT_ARGS_STRING" >> rufus.cmd
  write_out_rest_of_rufus_args
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
    # echo "NUM_JOBS_PLUS_ONE: $NUM_JOBS_PLUS_ONE"
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
    fi

    # Keeping until final job distribution schema settled on
    #else
    #  echo -e "INFO: $NUM_JOBS_BASE_COUNT slurm array jobs will be run with $BASE_COUNT_PER_SCRIPT rufus calls per script"
    #  echo -e "INFO: $NUM_JOBS_PLUS_ONE slurm array jobs will be run with $((BASE_COUNT_PER_SCRIPT + 1)) rufus calls per script"
    #  echo -e "INFO: 1 slurm array job will be run to combine results"
    #  echo -e "INFO: to fit into the allotted $SLURM_ARRAY_JOB_LIMIT_RUFUS_ARG jobs"
    #fi

    # Write out the slurm header
    ADJ_SLURM_ARRAY_END=$((ADJ_SLURM_ARRAY_LIMIT - 1))
    echo -e "#SBATCH --mem=${MEM_PER_JOB}" >> $RUFUS_SLURM_SCRIPT
    echo -e "#SBATCH --cpus-per-task=${CPUS_PER_JOB}" >> $RUFUS_SLURM_SCRIPT
    echo -e "#SBATCH -o ${WORKING_DIR}/slurm_out/rufus_call_%A_%a.out" >> $RUFUS_SLURM_SCRIPT
    echo -e "#SBATCH -e ${WORKING_DIR}/slurm_err/rufus_call_%A_%a.err" >> $RUFUS_SLURM_SCRIPT
    echo -e "#SBATCH -a 0-${ADJ_SLURM_ARRAY_END}%${SLURM_JOB_LIMIT_RUFUS_ARG}" >> $RUFUS_SLURM_SCRIPT
    printf '\n' >> $RUFUS_SLURM_SCRIPT

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
    echo -e "    region_arg=\$(singularity exec ${CONTAINER_PATH_RUFUS_ARG} bash ${RUFUS_ROOT}/singularity/launch_utilities/get_region.sh \"\$curr_job\" \"$WINDOW_SIZE_RUFUS_ARG\" \"$GENOME_BUILD_RUFUS_ARG\")" >> $RUFUS_SLURM_SCRIPT
    echo -e "    REGION_ARG=\"-R \$region_arg\"" >> $RUFUS_SLURM_SCRIPT

    # Per-region hash resolution (validated at setup time, glob guaranteed to match exactly one file)
    if [ -n "$KG1_HASH_DIR" ] || [ -n "$CONTROL_HASH_DIR" ]; then
      echo -e "    fmtd_region=\$(echo \"\$region_arg\" | tr ':-' '_')" >> $RUFUS_SLURM_SCRIPT
      echo -e "    HASH_ARGS=\"\"" >> $RUFUS_SLURM_SCRIPT
      if [ -n "$KG1_HASH_DIR" ]; then
        echo -e "    kg1_hash=\$(ls ${KG1_HASH_DIR}/*\${fmtd_region}*.Jhash)" >> $RUFUS_SLURM_SCRIPT
        echo -e "    HASH_ARGS=\"\$HASH_ARGS -e \$kg1_hash\"" >> $RUFUS_SLURM_SCRIPT
      fi
      if [ -n "$CONTROL_HASH_DIR" ]; then
        echo -e "    ctrl_hash=\$(ls ${CONTROL_HASH_DIR}/*\${fmtd_region}*.Jhash)" >> $RUFUS_SLURM_SCRIPT
        echo -e "    HASH_ARGS=\"\$HASH_ARGS -e \$ctrl_hash\"" >> $RUFUS_SLURM_SCRIPT
      fi
    fi
    echo -en "   srun --mem=${MEM_PER_JOB} singularity exec --bind ${BIND_MOUNTS}${DEV_BIND_ARGS} ${CONTAINER_PATH_RUFUS_ARG} bash /opt/RUFUS/runRufus.sh $SUBJECT_ARGS_STRING" >> $RUFUS_SLURM_SCRIPT
    echo -en "srun --mem=${MEM_PER_JOB} singularity exec --bind ${BIND_MOUNTS}${DEV_BIND_ARGS} ${CONTAINER_PATH_RUFUS_ARG} bash /opt/RUFUS/runRufus.sh $SUBJECT_ARGS_STRING" >> rufus.cmd
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
      echo -e "region_arg=\$(singularity exec ${CONTAINER_PATH_RUFUS_ARG} bash ${RUFUS_ROOT}/singularity/launch_utilities/get_region.sh \"\$SLURM_ARRAY_TASK_ID\" \"$WINDOW_SIZE_RUFUS_ARG\" \"$GENOME_BUILD_RUFUS_ARG\")" >> $RUFUS_SLURM_SCRIPT
      echo -e "REGION_ARG=\"-R \$region_arg\"" >> $RUFUS_SLURM_SCRIPT

      # Per-region hash resolution (validated at setup time, glob guaranteed to match exactly one file)
      if [ -n "$KG1_HASH_DIR" ] || [ -n "$CONTROL_HASH_DIR" ]; then
        echo -e "fmtd_region=\$(echo \"\$region_arg\" | tr ':-' '_')" >> $RUFUS_SLURM_SCRIPT
        echo -e "HASH_ARGS=\"\"" >> $RUFUS_SLURM_SCRIPT
        if [ -n "$KG1_HASH_DIR" ]; then
          echo -e "kg1_hash=\$(ls ${KG1_HASH_DIR}/*\${fmtd_region}*.Jhash)" >> $RUFUS_SLURM_SCRIPT
          echo -e "HASH_ARGS=\"\$HASH_ARGS -e \$kg1_hash\"" >> $RUFUS_SLURM_SCRIPT
        fi
        if [ -n "$CONTROL_HASH_DIR" ]; then
          echo -e "ctrl_hash=\$(ls ${CONTROL_HASH_DIR}/*\${fmtd_region}*.Jhash)" >> $RUFUS_SLURM_SCRIPT
          echo -e "HASH_ARGS=\"\$HASH_ARGS -e \$ctrl_hash\"" >> $RUFUS_SLURM_SCRIPT
        fi
      fi

      echo -en "srun --mem=${MEM_PER_JOB} singularity exec --bind ${BIND_MOUNTS}${DEV_BIND_ARGS} ${CONTAINER_PATH_RUFUS_ARG} bash ${RUFUS_ROOT}/runRufus.sh $SUBJECT_ARGS_STRING" >> $RUFUS_SLURM_SCRIPT
      echo -en "srun --mem=${MEM_PER_JOB} singularity exec --bind ${BIND_MOUNTS}${DEV_BIND_ARGS} ${CONTAINER_PATH_RUFUS_ARG} bash ${RUFUS_ROOT}/runRufus.sh $SUBJECT_ARGS_STRING" >> rufus.cmd
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
"#SBATCH --cpus-per-task=10"
"#SBATCH --mem=1G"
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

echo -e "srun singularity exec --bind ${BIND_MOUNTS}${DEV_BIND_ARGS} ${CONTAINER_PATH_RUFUS_ARG} bash ${RUFUS_ROOT}/post_process/post_process.sh -s ${SUBJECTS_RUFUS_ARG[0]} -w $WINDOW_SIZE_RUFUS_ARG" >> $PP_SLURM_SCRIPT

echo -en "##RUFUS_postProcessCommand=" >> rufus.cmd
echo -e "srun singularity exec --bind ${BIND_MOUNTS}${DEV_BIND_ARGS} ${CONTAINER_PATH_RUFUS_ARG} bash ${RUFUS_ROOT}/post_process/post_process.sh -s ${SUBJECTS_RUFUS_ARG[0]} -w $WINDOW_SIZE_RUFUS_ARG" >> rufus.cmd

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