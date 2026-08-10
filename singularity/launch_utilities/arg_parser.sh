#!/bin/bash

# Statics
DEFAULT_1MB_CPUS_PER_JOB="12"
DEFAULT_1MB_MEM_PER_JOB="20G"
DEFAULT_WG_CPUS_PER_JOB="40"
DEFAULT_WG_MEM_PER_JOB="150G"


usage() {
  echo "Usage: $0 [-s subject1,subject2,...] [-c control1,control2,control3...] [-b genome_build] [-a slurm_account] [-p slurm_partition] ...options"
  echo "Required Arguments:"
  echo "-s subject(s) A single subject or comma-delimited array of multiple subject BAM/CRAM/FASTQ/generator files (full paths)"
  echo "-c control(s) A single control or comma-delimited array of multiple controls (full paths, same formats as -s)"
  echo "              BAM, CRAM and generator inputs may be mixed; FASTQ inputs must be used alone."
  echo "-b genome_build  The desired genome build; currently only supports GRCh38"
  echo "-r reference  Full path to the reference file matching the genome build"
  echo "-a slurm_account  The account for the slurm job"
  echo "-p slurm_partition    The partition for the slurm job"
  echo "-l slurm_job_array_limit    The maximum amount of jobs slurm allows in an array"
  echo "Optional Arguments:"
  echo "-m kmer_depth_cutoff	The amount of kMers that must overlap the variant to be included in the final call set"
  echo "-w window_size	The size of the windows to run RUFUS on, in units of kilabases (KB); allowed range between 500-5000; defaults to single run of entire genome if not provided"
  echo "-f reference_hash   Full path to Jhash file containing reference kMer hash list"
  echo "-x exclude_hash     Single or comma-delimited list of full paths to Jhash file(s) containing kMers to exclude (static, same for all regions)"
  echo "-K kg1_hash_dir     Full path to directory of per-region KG1 Jhash files (files named *{region}*.Jhash)"
  echo "-G kg1_version      KG1 hash version to download from S3 (e.g., v3.0)"
  echo "-D ctrl_hash_dir    Full path to directory of per-region control Jhash files (files named *{region}*.Jhash)"
  echo "-V ctrl_version     Control hash version to download from S3 (e.g., v1.0)"
  echo "-y path_to_rufus_container   If not provided, will look in current directory for rufus.sif"
  echo "-z rufus_threads	Number of threads provided to RUFUS; defaults to 36 for entire genome; 10 for 1MB windows (NOTE: must be less than cpus_per_call)"
  echo "-e email  The email address to notify with slurm updates"
  echo "-q slurm_job_queue_limit    The maximum amount of jobs able to be ran at once; defaults to 20"
  echo "-t slurm_time_limit   The maximum amount of time to let the slurm job run; defaults to 7 days for full run, or one hour per window (DD-HH:MM:SS)"
  echo "-M memory_per_call    How much memory to allot to the rufus calling stage job; default 150G for entire genome; 20G for 1MB windows (e.g. 150G or 20G)"
  echo "-C cpus_per_call      How many cpus to allot to each rufus calling stage job; default 40 for entire genome; 12 for 1MB windows"
  echo "-d dev_binds     Comma-delimited list of host:container bind mounts for dev testing (e.g., /local/runRufus.sh:/opt/RUFUS/runRufus.sh)"
  echo "-P par_low_cov_threshold  Control k-mer count ceiling below which a variant is flagged as low-coverage-parent/inherited (default 7; set to 0 to disable, e.g. when using an assembly as the control)"
  echo "-H hash_size  jellyfish hash size (-s) for the k-mer count step, e.g. 64G. MUST match the -s any pre-built control/DSA/exclude hash was built with, or the merge fails. Maps to runRufus -hs; RUFUS defaults to 16G whole-genome / 1G windowed if unset."
  echo "-h help	Print usage"
  echo ""
  echo "Output files are written to the current working directory."
	exit 1
}

# Initialize variables
SUBJECTS_RUFUS_ARG=()
CONTROL_STRING_RUFUS_ARG=""
CONTROLS_RUFUS_ARG=()
GENOME_BUILD_RUFUS_ARG="GRCh38"
SLURM_ACCOUNT_RUFUS_ARG=""
SLURM_PARTITION_RUFUS_ARG=""
REFERENCE_RUFUS_ARG=""
KMER_DEPTH_CUTOFF_RUFUS_ARG="5"
WINDOW_SIZE_RUFUS_ARG="0"
EMAIL_RUFUS_ARG=""
SLURM_JOB_LIMIT_RUFUS_ARG="20"
SLURM_ARRAY_JOB_LIMIT_RUFUS_ARG="1000"
SLURM_TIME_LIMIT_RUFUS_ARG=""
CONTAINER_PATH_RUFUS_ARG=""
THREAD_LIMIT_RUFUS_ARG=""
EXCLUDE_HASH_LIST_RUFUS_ARG=()
REFERENCE_HASH_RUFUS_ARG=""
KG1_HASH_DIR=""
KG1_HASH_VERSION=""
CONTROL_HASH_DIR=""
CONTROL_HASH_VERSION=""
MEM_PER_JOB=""
CPUS_PER_JOB=""
PAR_LOW_COV_THRESHOLD_RUFUS_ARG="7"
DEV_BIND_MOUNTS_ARG=()
HASH_SIZE_RUFUS_ARG=""

# Parse command line options using getopts
#
# -C previously lacked its trailing colon ("M:CK:"), so it took no argument while its handler
# still assigned CPUS_PER_JOB=$OPTARG. `-C 36` therefore left OPTARG unset, never set
# CPUS_PER_JOB, and dropped "36" as an unread positional -- the documented cpus_per_call option
# silently did nothing and every job got the default (40 whole-genome / 12 windowed). Now `C:`.
#
# STILL MISWIRED: -h carries a colon ("h:") so it demands an argument. Bare `-h` never reaches
# the h) case; it falls to the missing-argument branch, printing "Option -h requires an argument"
# before the usage text. Usage still prints, so this is cosmetic. Fix is to drop the colon.
while getopts ":s:c:b:a:p:r:m:w:e:l:q:t:f:x:y:z:M:C:K:G:D:V:d:P:H:h" opt; do
    case ${opt} in
        s)
            IFS=',' read -r -a SUBJECTS_RUFUS_ARG <<< "$OPTARG"
            ;;
        c)
            IFS=',' read -r -a CONTROLS_RUFUS_ARG <<< "$OPTARG"
            ;;
        b)
            GENOME_BUILD_RUFUS_ARG=$OPTARG
            ;;
        a)
            SLURM_ACCOUNT_RUFUS_ARG=$OPTARG
            ;;
        p)
            SLURM_PARTITION_RUFUS_ARG=$OPTARG
            ;;
        r)
            REFERENCE_RUFUS_ARG=$OPTARG
            ;;
        m)
            KMER_DEPTH_CUTOFF_RUFUS_ARG=$OPTARG
            ;;
		w)
			WINDOW_SIZE_RUFUS_ARG=$OPTARG
			;;
        e)
            EMAIL_RUFUS_ARG=$OPTARG
            ;;
        q)
            SLURM_JOB_LIMIT_RUFUS_ARG=$OPTARG
            ;;
        l)
            SLURM_ARRAY_JOB_LIMIT_RUFUS_ARG=$OPTARG
            ;;
        t)
            SLURM_TIME_LIMIT_RUFUS_ARG=$OPTARG
            ;;
        y)
            CONTAINER_PATH_RUFUS_ARG=$OPTARG
            ;;
		x)
            IFS=',' read -r -a EXCLUDE_HASH_LIST_RUFUS_ARG <<< "$OPTARG"
			;;
        K)
            KG1_HASH_DIR=$OPTARG
            ;;
        G)
            KG1_HASH_VERSION=$OPTARG
            ;;
        D)
            CONTROL_HASH_DIR=$OPTARG
            ;;
        V)
            CONTROL_HASH_VERSION=$OPTARG
            ;;
		f)
			REFERENCE_HASH_RUFUS_ARG=$OPTARG
		    ;;
		z)
			THREAD_LIMIT_RUFUS_ARG=$OPTARG
			;;
        M)
			MEM_PER_JOB=$OPTARG
			;;
        C)
			CPUS_PER_JOB=$OPTARG
			;;
        d)
            IFS=',' read -r -a DEV_BIND_MOUNTS_ARG <<< "$OPTARG"
            ;;
        P)
            if ! [[ "$OPTARG" =~ ^[0-9]+$ ]]; then
                echo "ERROR: -P par_low_cov_threshold must be a non-negative integer." >&2
                exit 1
            fi
            PAR_LOW_COV_THRESHOLD_RUFUS_ARG=$OPTARG
            ;;
        H)
            if ! [[ "$OPTARG" =~ ^[0-9]+[GMKgmk]?$ ]]; then
                echo "ERROR: -H hash_size must be a jellyfish hash size, e.g. 64G, 500M, or a plain integer." >&2
                exit 1
            fi
            HASH_SIZE_RUFUS_ARG=$OPTARG
            ;;
        h)
            usage
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            usage
            ;;
        :)
            echo "Option -$OPTARG requires an argument." >&2
            usage
            ;;
    esac
done
shift $((OPTIND - 1))

# Check for required strings
if [[ ${#SUBJECTS_RUFUS_ARG[@]} -eq 0 || -z "$GENOME_BUILD_RUFUS_ARG" || -z "$REFERENCE_RUFUS_ARG" || -z "$SLURM_ACCOUNT_RUFUS_ARG" || -z "$SLURM_PARTITION_RUFUS_ARG" || -z "$SLURM_ARRAY_JOB_LIMIT_RUFUS_ARG" ]]; then
    echo "ERROR: Missing required argument(s); please see usage instructions with -h." >&2
    exit 1
fi

# Check that all subject files exist
for subject in "${SUBJECTS_RUFUS_ARG[@]}"; do
	if [ ! -f "$subject" ]; then
		echo "ERROR: subject file $subject does not exist or cannot be read." >&2
		exit 1
	fi
done

# Check that all control files exist
for control in "${CONTROLS_RUFUS_ARG[@]}"; do
	if [ ! -f "$control" ]; then
		echo "ERROR: control file $control does not exist or cannot be read." >&2
		exit 1
	else
		if [ "$CONTROL_STRING_RUFUS_ARG" == "" ]; then
			CONTROL_STRING_RUFUS_ARG="$control"
		else
			CONTROL_STRING_RUFUS_ARG="${CONTROL_STRING_RUFUS_ARG}, $control"
		fi
	fi
done

# Validate input file types.
#
# This used to require every subject and control to be the exact same type. That was broader than
# anything downstream actually needs: bam, cram and generator inputs all funnel into the same
# concatenated generator in runRufus.sh and the Filter stage streams that generator, so those three
# are indistinguishable by the time reads are pulled. Only FASTQ has to stand alone -- see below.
get_input_type() {
	case "$1" in
		*.cram) echo "cram" ;;
		*.bam) echo "bam" ;;
		*.fastq.gz|*.fq.gz) echo "fastq" ;;
		*.fastq|*.fq) echo "fastq" ;;
		*.generator) echo "generator" ;;
		*) echo "unknown" ;;
	esac
}

HAS_FASTQ_INPUT="false"
HAS_NON_FASTQ_INPUT="false"
HAS_GENERATOR_INPUT="false"
for input in "${SUBJECTS_RUFUS_ARG[@]}" "${CONTROLS_RUFUS_ARG[@]}"; do
	case "$(get_input_type "$input")" in
		fastq)
			HAS_FASTQ_INPUT="true"
			;;
		generator)
			HAS_GENERATOR_INPUT="true"
			HAS_NON_FASTQ_INPUT="true"
			;;
		bam|cram)
			HAS_NON_FASTQ_INPUT="true"
			;;
		*)
			echo "ERROR: input file $input has an unrecognized file type. Supported types: .bam, .cram, .generator, .fastq, .fq, .fastq.gz, .fq.gz" >&2
			exit 1
			;;
	esac
done

# FASTQ must be exclusive. A FASTQ subject populates _arg_fastqA/_arg_fastqB in runRufus.sh, and the
# Filter stage then reads ONLY those two mate files -- the generator holding the bam/cram/generator
# reads is never filtered, even though k-mer counting did span it. That yields a HashList whose
# k-mers have no reads to assemble from: fewer calls, no error message. Output naming desyncs as
# well, since post_process is handed subject[0] while runRufus names off the first non-FASTQ subject.
if [ "$HAS_FASTQ_INPUT" == "true" ] && [ "$HAS_NON_FASTQ_INPUT" == "true" ]; then
	echo "ERROR: FASTQ inputs cannot be combined with BAM/CRAM/generator inputs. RUFUS filters reads" >&2
	echo "       from the FASTQ mate files alone in that case, so reads from the other inputs would be" >&2
	echo "       counted but never filtered, silently costing calls. Pass all inputs as FASTQ, or" >&2
	echo "       convert the FASTQ to BAM/generator first. BAM, CRAM and generator inputs may be mixed." >&2
	exit 1
fi

# Check that reference file exists
if [ ! -f "$REFERENCE_RUFUS_ARG" ]; then
	echo "ERROR: reference file $REFERENCE_RUFUS_ARG does not exist or cannot be read." >&2
    exit 1
fi

# Check that BWA indexes exist alongside the reference.
# runRufus.sh prefers the extension-stripped prefix when <prefix>.sa is present and
# otherwise falls back to the reference path itself, so validate whichever it will pick.
# Without this, the missing index only surfaces at the bwa mem step, which in region
# mode is hours into every queued array task.
REFERENCE_BWA_PREFIX="$REFERENCE_RUFUS_ARG"
if [ -e "${REFERENCE_RUFUS_ARG%.*}.sa" ]; then
	REFERENCE_BWA_PREFIX="${REFERENCE_RUFUS_ARG%.*}"
fi

MISSING_BWA_INDEXES=()
for bwa_suffix in amb ann bwt pac sa; do
	if [ ! -e "${REFERENCE_BWA_PREFIX}.${bwa_suffix}" ]; then
		MISSING_BWA_INDEXES+=("${REFERENCE_BWA_PREFIX}.${bwa_suffix}")
	fi
done
if [ ! -e "${REFERENCE_RUFUS_ARG}.fai" ]; then
	MISSING_BWA_INDEXES+=("${REFERENCE_RUFUS_ARG}.fai")
fi

if [ ${#MISSING_BWA_INDEXES[@]} -ne 0 ]; then
	echo "ERROR: reference $REFERENCE_RUFUS_ARG is missing required index files:" >&2
	for missing in "${MISSING_BWA_INDEXES[@]}"; do
		echo "       $missing" >&2
	done
	echo "RUFUS aligns candidate reads with BWA and cannot run without these." >&2
	echo "Build them once with:" >&2
	echo "       bash \${RUFUS_ROOT}/resource_helpers/build_bwa_indexes.sh $REFERENCE_RUFUS_ARG" >&2
	echo "(indexing a human-sized reference takes roughly an hour)" >&2
	exit 1
fi

# Check that window size is in valid range
# Check if time limit has been assigned, if not - use defaults for full mode or windowed mode
if [ "$WINDOW_SIZE_RUFUS_ARG" -eq 0 ]; then
	if [ -z "$SLURM_TIME_LIMIT_RUFUS_ARG" ]; then
		SLURM_TIME_LIMIT_RUFUS_ARG="7-00:00:00"
	fi
    CPUS_PER_JOB=${CPUS_PER_JOB:-$DEFAULT_WG_CPUS_PER_JOB}
    MEM_PER_JOB=${MEM_PER_JOB:-$DEFAULT_WG_MEM_PER_JOB}
    THREAD_LIMIT_RUFUS_ARG=${THREAD_LIMIT_RUFUS_ARG:-36}
elif [ "$WINDOW_SIZE_RUFUS_ARG" -ne 1000 ]; then
	echo "ERROR: only windows of 1000 (1MB) supported currently" >&2
    exit 1
else
	if [ -z "$SLURM_TIME_LIMIT_RUFUS_ARG" ]; then
		SLURM_TIME_LIMIT_RUFUS_ARG="01:00:00"
	fi
    CPUS_PER_JOB=${CPUS_PER_JOB:-$DEFAULT_1MB_CPUS_PER_JOB}
    MEM_PER_JOB=${MEM_PER_JOB:-$DEFAULT_1MB_MEM_PER_JOB}
    THREAD_LIMIT_RUFUS_ARG=${THREAD_LIMIT_RUFUS_ARG:-10}
fi

# Guard the per-job memory against the jellyfish hash floor. jellyfish pre-faults the ENTIRE hash
# array at startup, so the -s size sets a hard RAM floor (not a peak that ramps): if -M is below it
# the count is OOM-killed in the first minute -- a fast, confusing failure that otherwise only shows
# up after submission. Catch it here instead.
#
# _hash_size_to_pow2_gb: parse a jellyfish -s token to the GiB of the power-of-two array jellyfish
# actually rounds up to (so -H 48G is judged as the 64G it allocates, not 48).
_hash_size_to_pow2_gb() {
    local tok="$1" num unit bytes p
    [[ "$tok" =~ ^([0-9]+)([GgMmKk]?)$ ]] || { echo 0; return; }
    num="${BASH_REMATCH[1]}"; unit="${BASH_REMATCH[2]}"
    case "$unit" in
        G|g) bytes=$(( num * 1024 * 1024 * 1024 ));;
        M|m) bytes=$(( num * 1024 * 1024 ));;
        K|k) bytes=$(( num * 1024 ));;
        *)   bytes=$num;;
    esac
    p=1
    while [ "$p" -lt "$bytes" ]; do p=$(( p * 2 )); done
    echo $(( p / (1024 * 1024 * 1024) ))
}
# Approx RAM (GiB) jellyfish -m 25 pre-faults for a hash of the given -s, from measured `jellyfish
# mem` points (16G->64, 32G->123, 64G->237, 128G->~460); 4 GiB per GiB-of-hash elsewhere, which is
# conservative (never lands under the true floor) for the large whole-genome sizes this guards.
_hash_mem_floor_gb() {
    local g; g=$(_hash_size_to_pow2_gb "$1")
    case "$g" in
        16) echo 64;;  32) echo 123;;  64) echo 237;;  128) echo 460;;
        *)  echo $(( g * 4 ));;
    esac
}
# Parse a SLURM memory string (e.g. 300G, 300000M, 1T) to GiB; 0 if it has no unit we recognize.
_slurm_mem_to_gb() {
    local tok="$1" num unit
    [[ "$tok" =~ ^([0-9]+)([GgMmTt])$ ]] || { echo 0; return; }
    num="${BASH_REMATCH[1]}"; unit="${BASH_REMATCH[2]}"
    case "$unit" in T|t) echo $(( num * 1024 ));; G|g) echo "$num";; M|m) echo $(( num / 1024 ));; esac
}

# Effective count-step hash size: -H override, else RUFUS's own default (16G whole-genome, 1G window).
if [ -n "$HASH_SIZE_RUFUS_ARG" ]; then
    _effective_hash_size="$HASH_SIZE_RUFUS_ARG"
elif [ "$WINDOW_SIZE_RUFUS_ARG" -eq 0 ]; then
    _effective_hash_size="16G"
else
    _effective_hash_size="1G"
fi
_hash_floor_gb=$(_hash_mem_floor_gb "$_effective_hash_size")
_mem_gb=$(_slurm_mem_to_gb "$MEM_PER_JOB")
if [ "$_mem_gb" -gt 0 ] && [ "$_hash_floor_gb" -gt "$_mem_gb" ]; then
    echo "ERROR: per-job memory (-M ${MEM_PER_JOB}) is below the jellyfish hash-size floor for -s ${_effective_hash_size}." >&2
    echo "       jellyfish pre-faults the whole ${_effective_hash_size} array (~${_hash_floor_gb} GiB) at startup and would be" >&2
    echo "       OOM-killed within the first minute of the count step." >&2
    echo "       Raise -M to at least $(( _hash_floor_gb + 40 ))G (headroom for reads + spill), or lower -H." >&2
    echo "       Whole-genome floors: -s 16G ~64, 32G ~123, 64G ~237 GiB." >&2
    exit 1
fi

# Generator inputs are whole-genome only, for the same reason FASTQ is: runRufus.sh cannot scope one
# to a region. A generator is an arbitrary shell command producing SAM, so `-R` is simply not applied
# to it -- the contents are used verbatim. In windowed mode every array task would therefore count
# and call the entire genome, burning thousands of node-hours to produce identical per-window VCFs.
if [ "$HAS_GENERATOR_INPUT" == "true" ] && [ "$WINDOW_SIZE_RUFUS_ARG" -ne 0 ]; then
	echo "ERROR: generator inputs cannot be used with windowed mode (-w); they are whole-genome only." >&2
	echo "       A generator cannot be region-scoped, so every window would re-run the whole genome." >&2
	echo "       Drop -w to run whole-genome, or supply the sample as an indexed BAM/CRAM." >&2
	exit 1
fi

if [ "$THREAD_LIMIT_RUFUS_ARG" -ge "$CPUS_PER_JOB" ]; then
	echo "ERROR: thread limit ($THREAD_LIMIT_RUFUS_ARG) must be less than cpus per job ($CPUS_PER_JOB)." >&2
	exit 1
fi
# Validate per-region hash flags: cannot specify both local dir and S3 version for same type
if [ -n "$KG1_HASH_DIR" ] && [ -n "$KG1_HASH_VERSION" ]; then
    echo "ERROR: Cannot specify both -K (local KG1 hash dir) and -G (KG1 S3 version). Use one or the other." >&2
    exit 1
fi
if [ -n "$CONTROL_HASH_DIR" ] && [ -n "$CONTROL_HASH_VERSION" ]; then
    echo "ERROR: Cannot specify both -D (local control hash dir) and -V (control S3 version). Use one or the other." >&2
    exit 1
fi

# Validate local hash directories exist if provided
if [ -n "$KG1_HASH_DIR" ] && [ ! -d "$KG1_HASH_DIR" ]; then
    echo "ERROR: KG1 hash directory does not exist: $KG1_HASH_DIR" >&2
    exit 1
fi
if [ -n "$CONTROL_HASH_DIR" ] && [ ! -d "$CONTROL_HASH_DIR" ]; then
    echo "ERROR: Control hash directory does not exist: $CONTROL_HASH_DIR" >&2
    exit 1
fi

# Check that if path to image not provided, it's in the current dir
if [ -z $CONTAINER_PATH_RUFUS_ARG ]; then
	if [ ! -f "rufus.sif" ]; then
		echo "ERROR: rufus.sif not in current directory - please provide path to container or put it in this one under rufus.sif"
	fi
fi

# Validate dev bind mounts and build singularity --bind args string
DEV_BIND_ARGS=""
if [ ${#DEV_BIND_MOUNTS_ARG[@]} -gt 0 ]; then
    for bind_spec in "${DEV_BIND_MOUNTS_ARG[@]}"; do
        host_path="${bind_spec%%:*}"
        container_path="${bind_spec#*:}"
        if [ "$host_path" == "$bind_spec" ]; then
            echo "ERROR: dev bind mount '$bind_spec' must be in host:container format (e.g., /local/runRufus.sh:/opt/RUFUS/runRufus.sh)" >&2
            exit 1
        fi
        if [ ! -e "$host_path" ]; then
            echo "ERROR: dev bind mount host path does not exist: $host_path" >&2
            exit 1
        fi
        DEV_BIND_ARGS+=" --bind ${bind_spec}"
    done
    echo "DEV MODE: additional bind mounts:${DEV_BIND_ARGS}"
fi

# Directories that must never be bind-mounted from the host, including anything beneath them.
# A bind shadows whatever the container has at that path, so binding /opt would hide the entire
# /opt/RUFUS install and binding /usr would swap the container's toolchain for the host's.
# Generators legitimately reference paths under these (e.g. /opt/RUFUS/scripts/FastqToSam.pl),
# so they have to be filtered out silently rather than treated as an error.
BIND_DENYLIST=(/ /bin /boot /dev /etc /lib /lib64 /opt /proc /root /run /sbin /srv /sys /usr /var)

is_denied_bind_dir() {
	local dir="$1" denied
	for denied in "${BIND_DENYLIST[@]}"; do
		if [ "$denied" == "/" ]; then
			[ "$dir" == "/" ] && return 0
		else
			[ "$dir" == "$denied" ] && return 0
			[[ "$dir" == "$denied"/* ]] && return 0
		fi
	done
	return 1
}

# Echo the directories referenced *inside* a generator file, one per line.
#
# A generator is a shell script RUFUS executes (`bash <generator>`) to produce SAM, so its data
# dependencies live in the file body rather than on the command line -- the launcher would
# otherwise bind the generator itself and none of the data it reads.
#
# Deliberately best-effort and strictly additive: a token is used only if it resolves to something
# that exists on the host. Paths assembled at runtime ("$DATA/sample.bam" yields the non-existent
# "/sample.bam"), globs, and process substitutions are skipped silently, leaving the bind set
# exactly as it was before. Nothing here can turn a working bind set into a broken one, and the
# container preflight in setup_slurm.sh is what catches whatever this misses.
generator_referenced_dirs() {
	local gen="$1"
	local token resolved dir
	while IFS= read -r token; do
		[ -n "$token" ] || continue
		[ -e "$token" ] || continue
		resolved="$(realpath "$token" 2>/dev/null)" || continue
		if [ -d "$resolved" ]; then
			dir="$resolved"
		else
			dir="$(dirname "$resolved")"
		fi
		is_denied_bind_dir "$dir" && continue
		echo "$dir"
	done < <(grep -o "/[^[:space:]'\";|&<>()\`]*" "$gen" 2>/dev/null | sort -u)
}

# Collect unique parent directories for all input files to use as bind mounts.
# Singularity --bind preserves host paths inside the container (no remapping needed).
collect_bind_dirs() {
    local -A seen_dirs
    local dirs=()

    # Always include pwd for output
    seen_dirs["$(pwd)"]=1
    dirs+=("$(pwd)")

    local files=("${SUBJECTS_RUFUS_ARG[@]}" "$REFERENCE_RUFUS_ARG")
    for control in "${CONTROLS_RUFUS_ARG[@]}"; do
        files+=("$control")
    done
    if [ -n "$REFERENCE_HASH_RUFUS_ARG" ]; then
        files+=("$REFERENCE_HASH_RUFUS_ARG")
    fi
    for exclude in "${EXCLUDE_HASH_LIST_RUFUS_ARG[@]}"; do
        files+=("$exclude")
    done

    # Add hash directories directly (not individual files)
    local hash_dirs=()
    if [ -n "$KG1_HASH_DIR" ]; then
        hash_dirs+=("$KG1_HASH_DIR")
    fi
    if [ -n "$CONTROL_HASH_DIR" ]; then
        hash_dirs+=("$CONTROL_HASH_DIR")
    fi

    for hd in "${hash_dirs[@]}"; do
        local resolved_hd
        resolved_hd="$(realpath "$hd")"
        if [ -z "${seen_dirs[$resolved_hd]+x}" ]; then
            seen_dirs["$resolved_hd"]=1
            dirs+=("$resolved_hd")
        fi
    done

    for f in "${files[@]}"; do
        local d
        d="$(dirname "$(realpath "$f")")"
        if [ -z "${seen_dirs[$d]+x}" ]; then
            seen_dirs["$d"]=1
            dirs+=("$d")
        fi
    done

    # Generator inputs carry their data references inside the file, so bind those dirs too.
    local gen_added=()
    for f in "${SUBJECTS_RUFUS_ARG[@]}" "${CONTROLS_RUFUS_ARG[@]}"; do
        [ "$(get_input_type "$f")" == "generator" ] || continue
        local gd
        while IFS= read -r gd; do
            [ -n "$gd" ] || continue
            if [ -z "${seen_dirs[$gd]+x}" ]; then
                seen_dirs["$gd"]=1
                dirs+=("$gd")
                gen_added+=("$gd")
            fi
        done < <(generator_referenced_dirs "$f")
    done
    # collect_bind_dirs runs inside a command substitution, so this cannot set a flag the caller
    # would see; the caller silences the repeat by setting _REPORTED_GEN_BINDS after the first call.
    if [ ${#gen_added[@]} -gt 0 ] && [ -z "${_REPORTED_GEN_BINDS:-}" ]; then
        echo "Binding directories referenced inside generator input(s): ${gen_added[*]}" >&2
        echo "  (best-effort scan; add any it missed with -d host:container)" >&2
    fi

    # Join with commas
    local IFS=','
    echo "${dirs[*]}"
}

BIND_MOUNTS="$(collect_bind_dirs)"
# setup_slurm.sh recomputes BIND_MOUNTS once the S3 hash dirs are known; the generator-bind notice
# above has already been shown, so suppress it on that second pass.
_REPORTED_GEN_BINDS=1

# Export variables for use in the main script
export BIND_MOUNTS
export SUBJECTS_RUFUS_ARG
export CONTROL_STRING_RUFUS_ARG
export CONTROLS_RUFUS_ARG
export GENOME_BUILD_RUFUS_ARG
export SLURM_ACCOUNT_RUFUS_ARG
export SLURM_PARTITION_RUFUS_ARG
export REFERENCE_RUFUS_ARG
export KMER_DEPTH_CUTOFF_RUFUS_ARG
export WINDOW_SIZE_RUFUS_ARG
export EMAIL_RUFUS_ARG
export SLURM_JOB_LIMIT_RUFUS_ARG
export SLURM_ARRAY_JOB_LIMIT_RUFUS_ARG
export SLURM_TIME_LIMIT_RUFUS_ARG
export CONTAINER_PATH_RUFUS_ARG
export THREAD_LIMIT_RUFUS_ARG
export EXCLUDE_HASH_LIST_RUFUS_ARG
export REFERENCE_HASH_RUFUS_ARG
export KG1_HASH_DIR
export KG1_HASH_VERSION
export CONTROL_HASH_DIR
export CONTROL_HASH_VERSION
export MEM_PER_JOB
export CPUS_PER_JOB
export PAR_LOW_COV_THRESHOLD_RUFUS_ARG
export DEV_BIND_ARGS
export HASH_SIZE_RUFUS_ARG