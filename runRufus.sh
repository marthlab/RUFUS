#!/bin/bash

# Set root directory to env path if exists
: "${RUFUS_ROOT:=/opt/RUFUS}"
echo "RUFUS_ROOT is $RUFUS_ROOT"

WORK_ROOT="${PWD}"
cd "$WORK_ROOT"

# Import globals
GLOBALS_FILE="$RUFUS_ROOT/resources/globals.txt" # Path to globals file inside container
set -a
source <(grep -v '^#' $GLOBALS_FILE | grep -v '^[[:space:]]*$' | sed 's/\r$//')
set +a
echo -n "You are running the $RUFUS_BRANCH"
echo " version of RUFUS: $RUFUS_VERSION"
echo " root dir is $WORK_ROOT"
# Provenance: the exact commit this image was built from (see /opt/RUFUS/BUILD_INFO, written at
# image build time). Absent on pre-provenance images, so guard for it. Read with sed rather than
# sourcing so it can't clobber RUFUS_VERSION/RUFUS_BRANCH from globals above.
if [ -f "$RUFUS_ROOT/BUILD_INFO" ]; then
	_gsha=$(sed -n 's/^GIT_SHA=//p'    "$RUFUS_ROOT/BUILD_INFO")
	_gbr=$(sed -n 's/^GIT_BRANCH=//p'  "$RUFUS_ROOT/BUILD_INFO")
	_gbt=$(sed -n 's/^BUILD_TIME=//p'  "$RUFUS_ROOT/BUILD_INFO")
	echo " image built from git ${_gsha:-unknown} (${_gbr:-unknown}) at ${_gbt:-unknown}"
fi

set -e 

start_time=$(date +"%s")
echo -e "RUFUS command was: $0 $@"
date

MaxHashDepth=1200; #need to make this a passed option
RDIR="$RUFUS_ROOT"
##########################__SET_EXECUTABLE_PATHS__##############################
# Internal scripts
RunJelly=$RDIR/scripts/RunJellyForRUFUS.sh
PullSampleHashes=$RDIR/scripts/CheckJellyHashList.sh
RUFUSOverlap=$RDIR/scripts/Overlap.shorter.sh
RemoveCoInheritedVars=$RDIR/scripts/remove_coinherited.sh

# RUFUS Specific Modules
RufAlu=$RDIR/bin/externals/rufalu/src/rufalu_project/src/aluDetect
modifiedJelly=$RDIR/bin/externals/modified_jellyfish/src/modified_jellyfish_project/bin/jellyfish

# TODO: already in PATH, shouldn't have to rename here, should be directly accessible
RUFUSmodel=$RDIR/bin/ModelDist
RUFUSfilterFASTQ=$RDIR/bin/RUFUS.Filter
RUFUSfilterFASTQse=$RDIR/bin/RUFUS.Filter.single

# Standard tools - TODO: move to image install
fastp=$RDIR/bin/externals/fastp/src/fastp_project/fastp
samblaster=$RDIR/bin/externals/samblaster/src/samblaster_project/samblaster
bwa=$RDIR/bin/externals/bwa/src/bwa_project/bwa
############################################################################################

die()
{
    local _ret=$2
    test -n "$_ret" || _ret=1
    test "$_PRINT_HELP" = yes && print_help >&2
    echo "$1" >&2
    exit ${_ret}
}

begins_with_short_option()
{
    local first_option all_short_options
    all_short_options='srtkmh'
    first_option="${1:0:1}"
    test "$all_short_options" = "${all_short_options/$first_option/}" && return 1 || return 0
}


# THE DEFAULTS INITIALIZATION - POSITIONALS
_positionals=()
_arg_exclude=()
# THE DEFAULTS INITIALIZATION - OPTIONALS
_arg_controls=()
_arg_subjects=()
_arg_subject_fastqs=()
_arg_control_fastqs=()
_arg_ref=
_arg_threads=3
_arg_kmersize=25
# NOTE: this default is load-bearing beyond its face value. The model phase (see MODEL PHASE
# below) only runs when _arg_min is EMPTY, so initialising it here means that branch is never
# taken and ModelDist never runs. Clearing this default would silently enable the model and
# change genotyping behaviour across the board -- do not "tidy" it without reading that block.
_arg_min=5
_arg_refhash=
_arg_saliva="FALSE"
_arg_exome="FALSE"
_MaxAlleleSize="1000"
_arg_mosaic="FALSE"
_assemblySpeed="full"
_parallel_jelly="no"
_pairedEnd="true"
_arg_region=
_use_region_hash="FALSE"
_arg_filterK=1
_arg_ParLowK=2
_arg_ParLowCovThreshold=7
_arg_hash_size=
_filterMinQ=15
_arg_stop="nope"
_arg_dev_reporting="FALSE"
_arg_dev_file_output="FALSE"
_arg_slurm_array_index=0
_arg_abs_coord_index=0
print_help ()
{
	printf "%s\n" "The general script's help msg"
	printf 'Usage: %s [-s|--subject <arg>] ... [-r|--ref <arg>] [-t|--threads <arg>] [-k|--kmersize <arg>] [-m|--min <arg>] [-h|--help] [<controls-1>] ... [<controls-n>] ...\n' "$0"
	printf "\t%s\n" "-s,--subject: bam/cram/fastq(or pair of fastq files)/generator file containing the subject of interest (can be used multiple times for split files from the same sample)"
	printf "\t%s\n" "-c, --controls: bam/cram/fastq(or pair of fastq files)/generator file for the sequence data of the control sample (can be used multipe times)"
	printf "\t%s\n" "-e,--exclude: Jhash file of kmers to exclude from mutation list, k must be  (no default, can be used multiple times)"
	printf "\t%s\n" "-eR, --exclude-region-hash: Region-specific Jhash file of kmers to exclude from mutation list; only support 1mb region sizes"
	printf "\t%s\n" "-se, --single_end_reads: subject bam file is single end reads, not paired (default is to assume paired end data)"
	printf "\t%s\n" "-r,--ref: file path to the desired reference file (no default)"
	printf "\t%s\n" "-cr,--cramref: file path to the desired reference file to decompress input cram files (no default)"
	printf "\t%s\n" "-t,--threads: number of threads to use (no default) (min 3)"
	printf "\t%s\n" "-k,--kersize: size of k-mer to use (no default)"
	printf "\t%s\n" "-m,--min: minimum k-mer count to call a variant (default 5 -- see NOTE below)"
	printf "\t%s\n" "-i, --saliva: flag to indicate that the subject sample is a buccal swab and likely contains a significant fraction of contaminant DNA"
	printf "\t%s\n" "-mx, --MaxAllele: Max size for insert/deletion events to put the entire alt sequence in. (default 1000)"
	printf "\t%s\n" "-L, --Report_Low_Freq: Reprot Mosaic/Low Frequency/Somatic variants (default FALSE)"
	printf "\t%s\n" "-CLEAN: Does not do a rufus run but cleans up intermediate files created by RUFUS"
	printf "\t%s\n" "-o, --devOutput: Prints very verbose run information to stdout"
	printf "\t%s\n" "-h,--help: Print help"
	printf "\t%s\n" "-d: Print dev help"
	printf "\n"
	printf "\t%s\n" "NOTE on -m/--min and the coverage model: -m always carries a value (default 5), and the"
	printf "\t%s\n" "coverage-model phase only runs when -m is unset, so in practice the model is never built."
	printf "\t%s\n" "RUFUS uses -m directly as the minimum k-mer count. A consequence is that the GT and"
	printf "\t%s\n" "FILTER columns of the VCF are not model-derived; see MODEL PHASE in runRufus.sh."
}

print_devhelp ()
{
	printf "%s\n" "The general script's help msg"
	printf 'Usage: %s [-s|--subject <arg>] ... [-r|--ref <arg>] [-t|--threads <arg>] [-k|--kmersize <arg>] [-m|--min <arg>] [-h|--help] [<controls-1>] ... [<control\
s-n>] ...\n' "$0"
	printf "\t%s\n" "-s,--subject: bam/cram/fastq(or pair of fastq files)/generator file containing the subject of interest (can be used multiple times for split files from the same sample)"
	printf "\t%s\n" "-c, --controls: bam/cram/fastq(or pair of fastq files)/generator file for the sequence data of the control sample (can be used multiple times)"
	printf "\t%s\n" "-e,--exclude: Jhash file of kmers to exclude from mutation list, k must be  (no default, can be used multiple times)"
	printf "\t%s\n" "-se, --single_end_reads: subject bam file is single end reads, not paired (default is to assume paired end data)"
	printf "\t%s\n" "-r,--ref: file path to the desired reference file (no default)"
	printf "\t%s\n" "-cr,--cramref: file path to the desired reference file to decompress input cram files (no default)"
	printf "\t%s\n" "-t,--threads: number of threads to use (no default) (min 3)"
	printf "\t%s\n" "-k,--kmersize: size of k-mer to use (no default)"
	printf "\t%s\n" "-m,--min: minimum k-mer count to call a variant (default 5 -- see NOTE below)"
	printf "\t%s\n" "-i, --saliva: flag to indicate that the subject sample is a buccal swab and likely contains a significant fraction of contaminant DNA"
	printf "\t%s\n" "-mx, --MaxAllele: Max size for insert/deletion events to put the entire alt sequence in. (default 1000)"
	printf "\t%s\n" "-L, --Report_Low_Freq: Report Mosaic/Low Frequency/Somatic variants (default FALSE)"
	printf "\t%s\n" "-z, --Dev output: Keep all intermediate files produced by RUFUS (default FALSE)"
	printf "\t%s\n" "-CLEAN: Does not do a rufus run but cleans up intermediate files created by RUFUS" 

	printf "\t%s\n" "################################################################################################"	
	printf "\t%s\n" "Extra options, mosstly experimental or algorithm parameters that you normaly dont need to adjust"
	printf "\t%s\n" "################################################################################################"

	printf "\t%s\n" "-f,--refhash: Jhash file containing reference hashList (no default)"
	printf "\t%s\n" "-mx, --MaxAllele: Max size for insert/deletion events to put the entire alt sequence in. (default 1000)"
	printf "\t%s\n" "-ex, --exome: flag to set if your input data is exome sequencing. Sets saliva fix and max kmer depth of 1 million (EXPERIMENTAL values used here have not been exhaustivly tested). The distribution model is not used -- though see the NOTE below, it is not used in any mode. The intended -m = 20 override is also inactive: -m already has a default, so an exome run uses -m 5 unless you pass -m explicitly."
	printf "\t%s\n" "-q1,--fastq1: If starting from fastq files, a list of the mate1 fastq files to improve RUFUS.filter"
	printf "\t%s\n" "-q2,--fastq2: If starting from fastq files, a list of the mate2 fastq files to improve RUFUS.filter"
	printf "\t%s\n" "-vs, --Very_Short_Assembly: use very short assembly methods, recommended when you are expecting over 10,000 variants "
	printf "\t%s\n" "-pj, --Parallelize_Jelly: parallelize jellyfish step, only use if you have more than 96G of ram"
	printf "\t%s\n" "-R, --Region: Run RUFUS only on a samtools style region"
	printf "\t%s\n" "-fk, --filterK: kmer threshold for number of kmers required to keep a read during filtering (default = 1)"
	printf "\t%s\n" "-fq, --filterMinQ: Minimum base quality for filter step, any kmer with any bases lower than this quality will be ignored (default = 15)"
	printf "\t%s\n" "-pl, --ParLowK: Lowest kmer count to be kept when counting parent jellyfish tables (default = 2, using 1 will SIGNIFICANTLY increase run time and is not advised)"
	printf "\t%s\n" "-plct, --ParLowCovThreshold: k-mer count ceiling in controls below which a variant is flagged as low-coverage-parent/inherited (default = 7, set to 0 to disable, e.g. when using an assembly as the control)"
	printf "\t%s\n" "-hs, --hash_size: jellyfish initial hash size (-s) for the subject/control count step, e.g. 32G (default: 1G in region mode, 8G in whole-genome mode). Must match the -s used to build any pre-made control/DSA hash being merged against."
	printf "\t%s\n" "-StJ: Stop run after jellyfish steps" #TODO: dont require reference and other non needed options if this is set
	printf "\t%s\n" "-StH: Stop run after hash compare steps" #TODO: dont require reference and other non needed options if this is set
	printf "\t%s\n" "-StF: Stop run after filter steps" 
	printf "\t\t%s\n" "This can be useful when you know you have low level contamination and want to remove kmers up to a certain count"
	printf "\t%s\n" "-h,--help: HELP!!!!!!!!!!!!!!!"
	printf "\t%s\n" "-d,--devhelp: HELP!!! for developers"
	printf "\t%s\n" "-pa,--passArray: pass the slurm array index to the script for debugging purposes"
	printf "\t%s\n" "-cn,--currAbsNum: pass the calculated absolute coordinate value to the script for debugging purposes"
	printf "\n"
	printf "\t%s\n" "################################################################################################"
	printf "\t%s\n" "NOTE: the coverage-distribution model (ModelDist) is never built, in any run mode."
	printf "\t%s\n" "################################################################################################"
	printf "\t%s\n" "The model phase is gated on '[ -z \$_arg_min ] && [ \$_arg_exome == FALSE ]', but _arg_min is"
	printf "\t%s\n" "initialised to 5 in the defaults block and never cleared, so the first test is never true and"
	printf "\t%s\n" "the else branch always runs. That branch writes a 4-line placeholder .7.7.model and never"
	printf "\t%s\n" "produces the .7.7.dist that RUFUS.interpret is passed via -mod, so ProcessDist fails to open"
	printf "\t%s\n" "it (non-fatally) and the Bayesian genotyper is inert. Observable effects: GT is '.' and FILTER"
	printf "\t%s\n" "is '.' on essentially every SNV/indel record, FILTER=PASS is unreachable on that path, and"
	printf "\t%s\n" "Dist1XCutoff falls back to 100000 which disables the repeat filter in PickDepthSomatic."
	printf "\t%s\n" "This is documented, not fixed: enabling the model activates several latent defects in"
	printf "\t%s\n" "RUFUS.interpret at the same time. See docs/RUFUS.interpret.audit.md section 2.1."

}
re='^[0-9]+$'; 
parse_commandline ()
{
	while test $# -gt 0
	do
	_key="$1"
	case "$_key" in
	-s|--subject)
		test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
		# Consume values up to the next -flag. FASTQ (unaligned) -> separate list, assembled into a
		# whole-genome generator below; bam/cram/generator -> _arg_subjects directly.
		while [[ $# -gt 1 && $2 != -* ]]; do
			_ext="${2##*.}"
			if [[ "$_ext" == "fastq" || "$_ext" == "fq" || "$_ext" == "gz" ]]; then
				_arg_subject_fastqs+=("$2")
			else
				_arg_subjects+=("$2")
			fi
			shift
		done
		;;
	-r|--ref)
		test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
		_arg_ref="$2"
		shift
		;;
			
	-cr|--cramref)
		test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
		_arg_cramref="$2"
		shift
		;;
		
	-q1|--fastq1)
		test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
		_arg_fastqA="$2"
		shift
		;;
	-q2|--fastq2)
		test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
		_arg_fastqB="$2"
		shift
		;;
	-t|--threads)
		test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
		_arg_threads="$2"
		if ! [[ $_arg_threads =~ $re ]] ; then
			echo "Threads must be a number "
			exit 100
		fi
		shift
		;;
			
	-f|--refhash)
		test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
		_arg_refhash="$2"
		shift
		;;
	-k|--kmersize)
		test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
		_arg_kmersize="$2"
		shift
		;;
	-c|--controls)
		test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
		while [[ $# -gt 1 && $2 != -* ]]; do
			_ext="${2##*.}"
			if [[ "$_ext" == "fastq" || "$_ext" == "fq" || "$_ext" == "gz" ]]; then
				_arg_control_fastqs+=("$2")
			else
				_arg_controls+=("$2")
			fi
			shift
		done
		;;
	-e|--exclude)
		test $# -lt 2 && die "Missing value for the optional argument '$key'." 1
		_arg_exclude+=("$2")
		shift
		;;
	-m|--min)
		test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
		_arg_min="$2"
		if ! [[ $_arg_min =~ $re ]] ; then
			echo "arg -m must be a number "
			exit 100
		fi
		shift
		;;
	-fk|--filterK)
		test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
		_arg_filterK=$2
		if ! [[ $_arg_filterK =~ $re ]] ; then
			echo "arg -fk or --filterK must be a number "
			exit 100
		fi
		shift
		;;
	-fq|--filterMinQ)
		test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
		_filterMinQ=$2
		if ! [[ $_filterMinQ =~ $re ]] ; then
			echo "arg -fq or --filterMinQ must be a number "
			exit 100
		fi
		shift
		;;	
	-pl|--ParLowK)
		test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
		_arg_ParLowK=$2
		if ! [[ $_arg_ParLowK =~ $re ]] ; then
			echo "arg -pl or --ParLowK must be a number "
			exit 100
		fi
		shift
		;;
	-plct|--ParLowCovThreshold)
		test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
		_arg_ParLowCovThreshold=$2
		if ! [[ $_arg_ParLowCovThreshold =~ $re ]] ; then
			echo "arg -plct or --ParLowCovThreshold must be a number "
			exit 100
		fi
		shift
		;;
	-hs|--hash_size)
		test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
		_arg_hash_size=$2
		if ! [[ $_arg_hash_size =~ ^[0-9]+[GMKgmk]?$ ]] ; then
			echo "arg -hs or --hash_size must be a jellyfish hash size, e.g. 32G, 500M, or a plain integer"
			exit 100
		fi
		shift
		;;
	-R|--region)
			test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
			_arg_region="$2"
			if  [[ -z $_arg_region  ]] ; then
				echo "arg region must not be empyty"
				exit 100
			fi
			shift
			;;
	-pa|--passArray)
		test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
			_arg_slurm_array_index=$2
			if ! [[ $_arg_slurm_array_index =~ $re ]] ; then
		echo "arg -pa or --passArray must be a number "
		exit 100
		fi
		shift
			;;
	-cn|--currAbsNum)
		test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
		_arg_abs_coord_index=$2
		if ! [[ $_arg_abs_coord_index =~ $re ]] ; then
		echo "arg -cn or --currAbsNum must be a number "
		exit 100
		fi
		shift
		;;
	-i|--saliva)
		_arg_saliva="TRUE"
		echo "INFO: Saliva subject sample provided"
		;;
	-vs|--Very_Short_Assembly|--vs)
		_assemblySpeed="veryfast"
		echo "INFO: Very fast assembly being used"
	;;
	-se|--single_end_reads|--se)
		_pairedEnd="false"
		echo "INFO: Sample Bam file is single end data"
		;;
	-pj|--Parallelize_Jelly|--pj)
		_parallel_jelly="yes"
		echo "INFO: Paralellizing jellyfish, assuming 3 samples"
		;; 
	-ex|--exome)
		_arg_exome="TRUE"
		echo "INFO: Exome run"
		;;
	-A|--MaxAllele)
		test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
		_MaxAlleleSize="$2"
		if ! [[ $_MaxAlleleSize =~ $re ]] ; then
			echo "MaxAlleleSize must be a number "
			exit 100
		fi
		echo "INFO: MaxAlleleSize set to $_MaxAlleleSize"
		shift
		;;
	-L|--Report_Low_Freq)
		_arg_mosaic="TRUE"
		echo "INFO: Reporting mosaic/low frequency variants"
		;;
	-h|--help)
		print_help
		exit 0
		;;
		-h*)
		print_help
		exit 0
		;;
	-d|--devhelp)
		print_devhelp
		exit 0
		;;
	-StJ)
		_arg_stop="jelly"; 
		echo "Stopping run after jellyfish steps"; 
		;;
	-StH)
		_arg_stop="hash"; 
		echo "Stopping run after hash compare steps";
		;;
	
	-StF)
		_arg_stop="filter"
		echo "Stopping run after filter steps"; 
		;;
	-o)
		_arg_dev_reporting="TRUE"
		echo "Verbose developer reporting to stdout"
		;;
	-z)
		_arg_dev_file_output="TRUE"
		echo "Retain all intermediate files created by RUFUS run"
		;;
	-CLEAN)
		echo "Cleaning up intermediate files";
		rm *generator.Jhash *generator.Jhash.histo *generator.Jhash.histo.7.7.dist *generator.Jhash.histo.7.7.out *generator.Jhash.histo.7.7.prob *generator.k25_c4.HashList *generator.Mutations.Mate1.fastq *generator.Mutations.Mate2.fastq *.generator.temp *.generator.temp.mate1.fastq *.temp.mate2.fastq *.generator.V2.overlap.fastq *.generator.V2.overlap.fastqd *.generator.V2.overlap.hashcount.fastq *.generator.V2.overlap.hashcount.fastq.bam.vcf *.generator.V2.overlap.hashcount.fastq.bam.vcf.bed;  
		echo "Cleanup done";
		exit 1; 
		;;
	*)
		echo "ERROR: Unkown argument $1"; 
		exit 100
		_positionals+=("$1")
		;;

	esac
	shift
	done
}

assign_positional_args ()
{
    _positional_names=()
    _our_args=$((${#_positionals[@]} - ${#_positional_names[@]}))
    for ((ii = 0; ii < _our_args; ii++))
    do
	_positional_names+=("_arg_exclude[$((ii + 0))]")
	done

    for (( ii = 0; ii < ${#_positionals[@]}; ii++))
    do
	eval "${_positional_names[ii]}=\${_positionals[ii]}" || die "Error during argument parsing, possibly an Argbash bug." 1
	done
}


# Combined EXIT trap: logs region status, then cleans up intermediary files
# NOTE: the actual trap is set later, after region logging variables are initialized

# This function wraps the Jellyfish hash table creation script in order to keep track of exit statuses.
# It writes all exit statuses for controls to a single file, jelly_exit_code_controls.log and the exit status for the subject to jelly_exit_code_subject.log
# Reports whether a hash table is empty if looking in a specific region to stdout.
make_jelly_hash ()
{
  local generator="$1"
  local k="$2"
  local threads="$3"
  local lowK="$4" # The minimum number of kmers to keep in count step
  local formRegionArg="$5"
  local isControl="$6"
  local exit_file="$7"

  # If we're in windowed mode, make hash smaller to make intersections with 1kg possible
  hash_size="16G"
  if [ "$formRegionArg" != "wg" ]; then
	echo "We're in windowed mode, use a smaller hash size to allow for 1kg comparison"
	hash_size="1G"
  fi

  # Allow the region-based default to be overridden by the -hs/--hash_size argument.
  # NOTE: this must match the -s used to build any pre-made control/DSA hash being merged against.
  if [ -n "$_arg_hash_size" ]; then
	echo "Overriding default hash size ($hash_size) with provided -hs/--hash_size value: $_arg_hash_size"
	hash_size="$_arg_hash_size"
  fi

  set +e
  bash $RunJelly "$generator" "$k" "$threads" "$lowK" "$hash_size"
  local exitCode=$?
  set -e

  echo "$exitCode" > "$exit_file"

  if [[ $exitCode -ne 0 && -n "$formRegionArg" ]]; then
      echo "RUFUS could not find any kmers in the provided region $formRegionArg in the file $generator; this usually means there is no coverage"
  fi
}

# Normalize a jellyfish -s size token (e.g. 64G, 500M, 1000000) to the power-of-two slot
# count jellyfish actually allocates, so two sizes compare the way "jellyfish merge" compares
# them. Echoes the normalized integer, or nothing if the token is unparseable.
normalize_hash_size ()
{
  local tok="$1" num unit bytes p
  [[ "$tok" =~ ^([0-9]+)([GgMmKk]?)$ ]] || { echo ""; return; }
  num="${BASH_REMATCH[1]}"; unit="${BASH_REMATCH[2]}"
  case "$unit" in
    G|g) bytes=$(( num * 1024 * 1024 * 1024 ));;
    M|m) bytes=$(( num * 1024 * 1024 ));;
    K|k) bytes=$(( num * 1024 ));;
    *)   bytes=$num;;
  esac
  p=1
  while [ "$p" -lt "$bytes" ]; do p=$(( p * 2 )); done
  echo "$p"
}

# Read the -s (hash size) a pre-built Jhash was created with, straight from its header via
# "jellyfish info", normalized to the allocated power-of-two. Empty if it can't be determined.
read_built_hash_size ()
{
  local hash="$1" tok
  tok=$($modifiedJelly info "$hash" 2>/dev/null | sed -n 's/.* -s \([0-9]\+[GMKgmk]\?\) .*/\1/p' | head -1)
  [ -n "$tok" ] && normalize_hash_size "$tok"
}

check_empty_hashes ()
{
	local region_arg="$1"
	shift
	local subject_file="${@: -1}"
	local control_files=("${@:1:$#-1}")

	# RunJellyForRUFUS.sh exit-code contract: 0 = counted OK, 1 = ran but the region
	# genuinely has no k-mers, 2 = the counting tool itself failed (OOM, disk, crash).
	# A tool failure must never be reported as an empty region -- in a sharded run that
	# records a lost shard as legitimately variant-free, which is a silent wrong answer
	# rather than a visible one. Treat anything that is not a clean 0 or 1 as failure,
	# so an unexpected code (e.g. a signal-derived 137) also fails loudly.
	local ef ef_rc
	for ef in "${control_files[@]}" "$subject_file"; do
		[ -f "$ef" ] || continue
		ef_rc="$(cat "$ef" 2>/dev/null)"
		case "$ef_rc" in
			0|1) ;;
			*)
				echo "ERROR: k-mer counting failed (exit '${ef_rc:-<empty>}', from $ef) in region $region_arg;" \
				     "this is a tool failure, NOT an absence of coverage" >&2
				_region_exit_reason="jellyfish_failed"
				exit 1
				;;
		esac
	done

	# Check that at least one control has hashes (i.e. has a zero exit code)
	found_zero=false

	# Controls
	if [ "${#control_files[@]}" -gt 0 ]; then
		found_zero=false
		for f in "${control_files[@]}"; do
		if [ -f "$f" ] && [ "$(cat "$f")" -eq 0 ]; then
			found_zero=true
			break
		fi
		done

		if [ "$found_zero" = false ]; then
		echo "No control kmers found in region $region_arg" >&2
		_region_exit_reason="no_control_kmers"
		exit 0
		fi
	fi

	# Subject
	if [ ! -f "$subject_file" ] || [ "$(cat "$subject_file")" -ne 0 ]; then
		echo "No subject kmers found in region $region_arg" >&2
		_region_exit_reason="no_subject_kmers"
		exit 0
	fi

    # Cleanup
	rm -f "$subject_file"
	for f in "${control_files[@]}"; do
		rm -f "$f"
	done
}


parse_commandline "$@"
region_postfix=""
if [ ! -z "${_arg_region}" ]; then
	formatted_region=$(echo "${_arg_region}" | tr : _ | tr - _)
	region_postfix=".${formatted_region}"
else
	formatted_region="wg"
	region_postfix=".${formatted_region}"
fi

# Region status logging: each invocation appends a line to a shared log in WORK_ROOT.
# The EXIT trap handles all logging so we only need to set _region_exit_reason before exiting.
REGION_LOG="${WORK_ROOT}/region_status.log"
_region_exit_reason=""

on_exit() {
  local exit_code=$?

  # --- Region status logging (only in region mode) ---
  if [ "$formatted_region" != "wg" ] && [ -n "$formatted_region" ]; then
    local status reason
    if [ "$exit_code" -eq 0 ]; then
      if [ "$_region_exit_reason" = "success" ]; then
        status="VARIANTS_CALLED"
      else
        status="NO_VARIANTS"
      fi
      reason="${_region_exit_reason:-unknown}"
    else
      status="ERROR"
      reason="exit_code=${exit_code};${_region_exit_reason:-unknown}"
    fi
    printf '%s\t%s\t%s\n' "$formatted_region" "$status" "$reason" >> "$REGION_LOG"
  fi

  # --- File cleanup ---
  echo "Cleaning up..." >&2
  if [ "$exit_code" -ne 0 ]; then
    echo "Script exited with error (exit code $exit_code). Preserving WORK_DIR for debugging: $WORK_DIR" >&2
    return
  fi
  if [ "$_arg_dev_file_output" == "FALSE" ]; then
    # Must cd out of WORK_DIR before removing it (rm -rf silently fails on CWD in some filesystems/containers)
    cd "$WORK_ROOT"
    if [[ "$WORK_DIR" == "$WORK_ROOT"/rufus_* ]]; then
      rm -rf "$WORK_DIR"
    else
      echo "Refusing to remove unsafe WORK_DIR: $WORK_DIR" >&2
    fi
  fi
}
trap 'on_exit' EXIT

# Make region specific directory so no file overlaps
WORK_DIR="${WORK_ROOT}/rufus_${formatted_region}"
export WORK_DIR="$WORK_DIR"

# Make all sub-dirs w/ abs path
mkdir -p "$WORK_DIR"
RUFUS_TMP="$WORK_DIR/rufus_temp"
mkdir -p "$RUFUS_TMP"
mkdir -p "$WORK_DIR/Intermediates"
mkdir -p "$WORK_DIR/TempOverlap"

cd "$WORK_DIR"

rufus_invoc_file="$RUFUS_TMP/rufus_command_$formatted_region.txt"
echo "$RUFUS_BRANCH" > $rufus_invoc_file
echo "$RUFUS_VERSION" >> $rufus_invoc_file
printf "%q " "$@"     >> "$rufus_invoc_file"

# [ <-- needed because of Argbash

##############################__Print out all parameters__#################################################
#echo "~~~~~~~~~~~~ printing out paramater values used in script ~~~~~~~~~~~~~~~~"
#echo "  _arg_exclude=:"
#	for each in "${_arg_exclude[@]}"                           
#	do                                                         
#	  	echo "		$each"
#	done           
#echo "  _arg_controls=:"
#    arg_control_string=""
#	for each in "${_arg_controls[@]}"                           #
#	do                                                          #
#		echo "		$each"
#        arg_control_string="$arg_control_string $each"                                              #
#	done           

if [ "$_arg_dev_reporting" = "TRUE" ]; then
	echo "Verbose developer reporting on..."
	echo "  _arg_subjects=${_arg_subjects[*]}"
	echo "  _arg_ref=$_arg_ref" 
	echo "  _arg_threads=$_arg_threads" 
	echo "  _arg_kmersize=$_arg_kmersize" 
	echo "  _arg_min=$_arg_min" 
	echo "  _arg_refhash=$_arg_refhash" 
	echo "  _arg_saliva=$_arg_saliva" 
	echo "  _arg_exome=$_arg_exome" 
	echo "  _MaxAlleleSize=$_MaxAlleleSize" 
	echo "  _arg_mosaic=$_arg_mosaic" 
	echo "  _assemblySpeed=$_assemblySpeed" 
	echo "  _parallel_jelly=$_parallel_jelly" 
	echo "  _pairedEnd=$_pairedEnd" 
	echo "  _arg_region=$_arg_region" 
	echo "  _arg_filterK=$_arg_filterK" 
	echo "  _arg_ParLowK=$_arg_ParLowK"
	echo "  _arg_hash_size=$_arg_hash_size"
	echo "  _filterMinQ=$_filterMinQ"
	echo "  _arg_dev_file_output=$_arg_dev_file_output"
fi

##############################__CHECK_FOR_MANDATORY_PARAMS__#################################################
if [ -z $_arg_kmersize ]
then
 echo "INFO: You did not profide a kmer size, running with the default 25"
 _arg_kmersize="25"
 #echo "You must provide a kmer size [--kmersize|-k] (we recommend a kmer size of 25)"
 #echo "Killing run with non-zero exit status"
 #kill -9 $$
fi

if [ -z $_arg_threads ]
then
    echo "INFO: no thread argument given, running with all avaialble threads = $(nproc)"; 
    _arg_threads=$(nproc); 
fi

if [ ${#_arg_subjects[@]} -eq 0 ] && [ ${#_arg_subject_fastqs[@]} -eq 0 ]
then
	echo "ERROR: you must provide at least one subject sample (sample you want to call variants in)"
	kill -9 $$
fi

if [ ${#_arg_exclude[@]} -eq "0" ] && [ ${#_arg_controls[@]} -eq "0" ] && [ ${#_arg_control_fastqs[@]} -eq "0" ]
then
    echo "You must provide RUFUS with at least one control or exclude sample"
    echo "Killing run with non-zero exit status"
    kill -9 $$
fi


#########__remove -e and --exclude from _arg_exclude array__################
new_array=()
for value in "${_arg_exclude[@]}"
do
    [[ $value != --exclude ]] && new_array+=($value)
done
ExcludeTemp=("${new_array[@]}")
unset new_array

new_array=()
for value in "${ExcludeTemp[@]}"
do 
    [[ $value != -e ]] && new_array+=($value)
done
_arg_exclude=("${new_array[@]}")
unset new_array
unset ExcludeTemp
########################Setting up Exome Run EXPERIMENTAL ##################################

if [ "$_arg_exome" == "TRUE" ]; then 
	echo "Exome run set.  Setting max kmer to 1M and saliva = true and making sure a lower cutoff was set "
	MaxHashDepth=100000000
	_arg_saliva="TRUE"

	# NOTE: this block is dead for the same reason the model phase is (see MODEL PHASE below):
	# _arg_min is initialised to 5 in the defaults block and never cleared, so it is never
	# empty and the 20 is never applied. An exome run therefore uses -m 5 like everything
	# else unless the user passes -m explicitly. Documented, not changed -- raising the
	# effective exome minimum from 5 to 20 is a real behaviour change, not a cleanup.
	if [ -z $_arg_min ]; then
		echo "Minimum not provided, picking a min of 20 for the alt count"
		_arg_min="20"
	fi
fi


#############################################################################################################

###############__PRINTING_OUT_ARG_BASH_VALUES__##############
#echo "Value of --subject: $_arg_subject"                    #
#echo "Value of --controls:"                                 #
#for each in "${_arg_controls[@]}"                           #
#do                                                          #
#  echo "$each"                                              #
#done                                                        #
#echo "Value of --threads: $_arg_threads"                    #
#echo "Value of --kmersize: $_arg_kmersize"                  #
#echo "Value of --ref: $_arg_ref"                            #
#echo "Value of --minCount $_arg_min"                        #
#############################################################

Parents=("${_arg_controls[@]}")

#########__CREATE_ALL_GENERATOR_FILES_AND_VARIABLES__#############
# Use first subject file for naming conventions
if [ ${#_arg_subjects[@]} -gt 0 ]; then
	ProbandFileName=$(basename "${_arg_subjects[0]}")
else
	ProbandFileName=$(basename "${_arg_subject_fastqs[0]}")
fi
ProbandExtension="${ProbandFileName##*.}"
ProbandGenerator="${ProbandFileName}${region_postfix}.generator"

# Build concatenated generator from all subject files.
# Exactly one command in the generator may emit a SAM header: the body is run as a single
# stream (`bash "$ProbandGenerator" | samtools ...`) and samtools aborts on a second @HD
# mid-stream -- collate discards the whole stream, so Filter sees zero reads. _subj_hdr
# carries the header flag for the first emitting command and is cleared thereafter; the
# FASTQ block below continues the same flag so a bam+fastq mix stays single-headered.
_subj_hdr="-h "
> "$ProbandGenerator"
for subject in "${_arg_subjects[@]}"
do
    subjectFileName=$(basename "$subject")
    subjectExtension="${subjectFileName##*.}"

    if [[ "$subjectExtension" != "cram" ]] && [[ "$subjectExtension" != "bam" ]] && [[ "$subjectExtension" != "generator" ]] || [[ ! -e "$subject" ]]
    then
        echo "The proband bam/cram/generator file $subject was not provided or does not exist; killing run with non-zero exit status"
        kill -9 $$
    elif [[ "$subjectExtension" == "bam" ]]
    then
        if [[ ! -e "$subject".bai ]]
        then
            echo "Index file for subject bam file $subject not found. Please place in data directory and rerun."
            _region_exit_reason="missing_subject_bam_index"
            exit 1
        fi
        echo "samtools view ${_subj_hdr}-@ 8 -F 3328 $subject $_arg_region" >> "$ProbandGenerator"
        _subj_hdr=""
    elif [[ "$subjectExtension" == "cram" ]]
    then
        if [[ ! -e "$subject".crai ]]
        then
            echo "Index file for subject cram file $subject not found. Please place in data directory and rerun."
            _region_exit_reason="missing_subject_cram_index"
            exit 1
        fi
        if [ "$_arg_cramref" == "" ]
        then
            echo "ERROR cram reference not provided for cram input"
            kill -9 $$
        fi
        echo "samtools view ${_subj_hdr}-@ 8 -F 3328 -T $_arg_cramref $subject $_arg_region" >> "$ProbandGenerator"
        _subj_hdr=""
        _arg_ref="$_arg_cramref"
    elif [[ "$subjectExtension" == "generator" ]]
    then
        # A pre-built generator carries its own header-emitting command. Keep it only if it
        # lands first; otherwise strip the header flags as it is appended.
        if [ -n "$_subj_hdr" ]
        then
            cat "$subject" >> "$ProbandGenerator"
        else
            sed -e 's/^\(samtools view\) -h /\1 /' -e 's/ header$//' "$subject" >> "$ProbandGenerator"
        fi
        _subj_hdr=""
    else
        echo "unknown error during generator generation, killing run with non-zero exit status"
        kill -9 $$
    fi
done

# FASTQ subject(s): whole-genome only -- unaligned reads cannot be region-scoped. The loop above
# skipped them; build the generator here as unmapped SAM records, headed by a single @HD if no
# bam/cram subject above has already claimed it (_subj_hdr).
if [ ${#_arg_subject_fastqs[@]} -gt 0 ]; then
	if [ -n "$_arg_region" ]; then
		echo "ERROR: FASTQ input is whole-genome only and cannot be region-scoped; remove -R/--region (or supply an aligned bam/cram)."
		_region_exit_reason="fastq_with_region"
		exit 1
	fi
	if [ -z "$_arg_ref" ]; then
		echo "ERROR: FASTQ subject input requires a reference via -r/--ref."
		exit 1
	fi
	for fq in "${_arg_subject_fastqs[@]}"; do
		[ -e "$fq" ] || { echo "FASTQ subject file $fq does not exist; killing run"; kill -9 $$; }
		_hdr=""; [ -n "$_subj_hdr" ] && _hdr=" header"; _subj_hdr=""
		if [[ "$fq" == *.gz ]]; then
			echo "perl $RDIR/scripts/FastqToSam.pl <(zcat $fq)$_hdr" >> "$ProbandGenerator"
		else
			echo "perl $RDIR/scripts/FastqToSam.pl <(cat $fq)$_hdr" >> "$ProbandGenerator"
		fi
	done
	# Paired FASTQ: filter directly from the raw mate files via the -q1/-q2 path (avoids the
	# collate/mate-split, which needs pairing flags the counting generator does not carry). Single-end
	# (-se) skips this and uses the RUFUS.Filter.single path on the generator. An explicit -q1/-q2 wins.
	if [ "$_pairedEnd" = "true" ] && [ -z "${_arg_fastqA:-}" ]; then
		if [ ${#_arg_subject_fastqs[@]} -lt 2 ]; then
			echo "ERROR: paired FASTQ subject needs mate1 and mate2 files (-s R1 R2); for single-end reads add -se."
			exit 1
		fi
		_arg_fastqA="${_arg_subject_fastqs[0]}"
		_arg_fastqB="${_arg_subject_fastqs[1]}"
	fi
fi

# Have to do this after proband check in case cram reference is used
_arg_ref_cat="${_arg_ref%.*}"

ParentGenerators=()
ParentJhash=()
ParentFileNames=""
space=" "

for parent in "${Parents[@]}"
do 
    parentFileName=$(basename "$parent")
    ParentFileNames=$ParentFileNames$space$parent
    parentExtension="${parentFileName##*.}"

    if  [[ "$parentExtension" != "cram" ]] && [[ "$parentExtension" != "bam" ]]  && [[ "$parentExtension" != "generator" ]] 
    then
		echo "The control bam/generator file" "$parent" " was not provided, or does not exist; killing run with non-zero exit status"
		kill -9 $$
    elif [[ "$parentExtension" == "bam" ]]
    then
		# check for index file (needed for mpileup in post processing)
		if [[ ! -e "$parent".bai ]]
		then
			echo "Index file for control bam file $parentFileName not found. Please place in data directory and rerun."
			_region_exit_reason="missing_control_bam_index"
			exit 1
		fi
	    	parentGenerator="${parentFileName}${region_postfix}.generator"
	    	ParentGenerators+=("$parentGenerator")
	    	echo "samtools view -h -@ 8 -F 3328 $parent $_arg_region" > "$parentGenerator"
    elif [[ "$parentExtension" == "cram" ]] 
    then
		# check for index file (needed for mpileup in post processing)
		if [[ ! -e "$parent".crai ]]
		then
			echo "Index file for control cram file $parentFileName not found. Please place in data directory and rerun."
			_region_exit_reason="missing_control_cram_index"
			exit 1
		fi
		parentGenerator="${parentFileName}${region_postfix}.generator"
		ParentGenerators+=("$parentGenerator")
	    if [ "$_arg_cramref" == "" ]
	    then
			echo "ERROR cram reference not provided for cram input"; 
			kill -9 $$ 
	    fi
		echo "samtools view -h -@ 8 -F 3328 -T $_arg_cramref $parent  $_arg_region" > "$parentGenerator"
		_arg_ref="$_arg_cramref"
    elif [[ "$parentExtension" = "generator" ]]
    then
		# The historical name is <control><region_postfix> as a FULL path, and it is load-bearing:
		# RunJellyForRUFUS.sh early-returns when $GEN.Jhash exists, so a pre-built control hash
		# placed next to the input as <control><region_postfix>.Jhash (e.g. a DSA hash symlinked to
		# DSA_SMHT004.1.generator.wg.Jhash) is picked up and jellyfish is skipped entirely. The
		# generator body is never executed in that case, which is the point -- the hash already
		# exists and the reads it came from may not even be on this filesystem.
		#
		# A caller may equally supply a pre-scoped generator at that path. Only when neither is
		# present does the generator actually have to run, and only then is a runnable copy
		# materialised -- in the working directory, not next to the user's input.
		#
		# The generator is used verbatim: -R/--region is NOT applied to it, exactly as for generator
		# subjects above. Scoping a generator to a region is the caller's responsibility.
		parentGenerator="${parent}${region_postfix}"
		if [ ! -e "$parentGenerator" ] && [ ! -e "${parentGenerator}.Jhash" ]
		then
			if [[ ! -e "$parent" ]]
			then
				echo "The control generator file $parent does not exist; killing run with non-zero exit status"
				kill -9 $$
			fi
			parentGenerator="${parentFileName}${region_postfix}.generator"
			cat "$parent" > "$parentGenerator"
		fi
		ParentGenerators+=("$parentGenerator")
    fi
done
#################################################################

# FASTQ control: whole-genome only; a single control built from its mate fastq(s).
if [ ${#_arg_control_fastqs[@]} -gt 0 ]; then
	if [ -n "$_arg_region" ]; then
		echo "ERROR: FASTQ control input is whole-genome only; remove -R/--region."
		exit 1
	fi
	ctrlFqGen="$(basename "${_arg_control_fastqs[0]}")${region_postfix}.generator"
	ParentGenerators+=("$ctrlFqGen")
	> "$ctrlFqGen"
	_fq_first=1
	for fq in "${_arg_control_fastqs[@]}"; do
		[ -e "$fq" ] || { echo "FASTQ control file $fq does not exist; killing run"; kill -9 $$; }
		_hdr=""; [ "$_fq_first" -eq 1 ] && _hdr=" header"; _fq_first=0
		if [[ "$fq" == *.gz ]]; then
			echo "perl $RDIR/scripts/FastqToSam.pl <(zcat $fq)$_hdr" >> "$ctrlFqGen"
		else
			echo "perl $RDIR/scripts/FastqToSam.pl <(cat $fq)$_hdr" >> "$ctrlFqGen"
		fi
	done
fi

# Note: BWA index presence is validated below, once _arg_ref_bwa is resolved.

###### when we add PB need to check its reference stuff here 
###########################################################################

refFileName=$(basename "$parent")
refExtension="${refFileName##*.}"

if [[ "$refExtension" != fa ]] && [[ -e "$_arg_ref".fa ]] 
then 
    _arg_ref="$_arg_ref".fa
elif [[ "$refExtension" != fasta ]] && [[ -e "$_arg_ref".fasta ]]
then
    _arg_ref="$_arg_ref".fasta
fi

if [[ ! -e "$_arg_ref" ]]
then
    echo "cannot find reference file " "$_arg_ref"
    echo "tried extensions: "
    echo "no extension "
    echo ".fa"
    echo ".fasta"
    echo "killing run with non-zero exit status"
    kill -9 $$
fi
if [[ -e "$_arg_ref_cat".sa ]]
then
    _arg_ref_bwa=$_arg_ref_cat
else
    _arg_ref_bwa=$_arg_ref
fi

# Validate BWA/samtools indexes before any expensive stage. The launch scripts
# check this too, but runRufus.sh can be invoked directly, and without this the
# missing index only surfaces at the bwa mem step -- far into the run, after the
# jellyfish, model, hashlist and filter stages. _arg_ref_bwa above resolved the
# exact prefix bwa will load; the .fai hangs off the reference samtools reads.
missing_bwa_indexes=()
for bwa_suffix in amb ann bwt pac sa
do
    if [[ ! -e "$_arg_ref_bwa".$bwa_suffix ]]
    then
        missing_bwa_indexes+=("$_arg_ref_bwa.$bwa_suffix")
    fi
done
if [[ ! -e "$_arg_ref".fai ]]
then
    missing_bwa_indexes+=("$_arg_ref.fai")
fi
if [[ ${#missing_bwa_indexes[@]} -ne 0 ]]
then
    echo "ERROR: reference $_arg_ref is missing required index files:"
    for missing in "${missing_bwa_indexes[@]}"
    do
        echo "       $missing"
    done
    echo "RUFUS aligns candidate reads with BWA and cannot run without these."
    echo "Build them with:"
    echo "       bash $RDIR/resource_helpers/build_bwa_indexes.sh $_arg_ref"
    _region_exit_reason="missing_bwa_indexes"
    exit 1
fi


#echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
#echo "Reference supplied: " "$_arg_ref"
#echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"

# echo -e "slurm array index is: $_arg_slurm_array_index"
# echo -e "absolute coordinate index is: $_arg_abs_coord_index"
# echo -e "slurm array index is: $_arg_slurm_array_index" >&2
# echo -e "absolute coordinate index is: $_arg_abs_coord_index" >&2

################__COPY_ARG_BASH_VARIABLES_TO_SCRIPT_VARIABLES__##################

K=$_arg_kmersize
Threads=$_arg_threads
ref=$_arg_ref
#################################################################################

if [[ -z "$K" ]]
then
    echo "@@@@@@@@@@@__WARNING__@@@@@@@@@@@@@"
    echo "kmer size ([-k|kmersize]) was not provided, killing run with non-zero exit status"
    echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
    kill -9 $$

fi

if [[ -z "$Threads" ]]
then
    echo "@@@@@@@@@@@__WARNING__@@@@@@@@@@@@@"
    echo "number of threads ([-t|--threads]) was not provided, killing run with non-zero exit status"
    echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
    kill -9 $$
fi

if [[ -z "$ref" ]]
then
    echo "@@@@@@@@@@@__WARNING__@@@@@@@@@@@@@"
    echo "reference file ([-r|--ref]) was not provided, killing run with non-zero exit status"
    echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
    kill -9 $$

fi


###################__PRINT_VARIABLES_USED__######################################
#echo "~~~~~~~~~~~~ printing out parameter values used in script ~~~~~~~~~~~~~~~~"
#echo "value of ProbandGenerator $ProbandGenerator"
#echo "Value of ParentGenerators:"
#for parent  in "${ParentGenerators[@]}"
#do
#  echo " $parent"
#done
#echo "Value of K is: $K"
#echo "value of ref is: $ref"
#echo "value of min is: $_arg_min" 
#echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
#################################################################################


if [ "$_arg_refhash" ]
then
    echo "Reference hash provided: " "$_arg_refhash"
fi

if [ "$_arg_min" ]
then
      MutantMinCov=$_arg_min
fi
######################################


############__BUILD_JHASH_STRING__################
parentsString=""
parentsExcludeString=""
space=" "
jhash=".Jhash"

for parent in "${ParentGenerators[@]}"
do
  #echo "parent is  $parent "
  parentsString=$parentsString$space$parent$jhash
done
for exclude in "${_arg_exclude[@]}"
do
    parentsExcludeString=$parentsExcludeString$space$exclude
done

##################################################


############__PREFLIGHT: HASH SIZE MATCH__################
# Jellyfish can only merge/diff hashes built at the same -s. The subject hash is built fresh here
# (hours for a whole-genome sample), but pre-built control/DSA/exclude hashes carry a fixed size
# from when they were made. If those disagree with the subject size, "$modifiedJelly merge" aborts
# with "Can't merge hash with different size", leaving an empty HashList that only surfaces much
# later as the misleading "No mutant hashes pulled from fastqs". Catch it here, in seconds, before
# paying for the subject build.
#
# Only pre-built hashes that already exist on disk can mismatch; control generators counted in this
# run are built at the subject size and match by construction, so those (not yet on disk) are skipped.

# Intended subject hash size -- mirror make_jelly_hash: -hs override, else 16G whole-genome / 1G region.
if [ -n "$_arg_hash_size" ]; then
	_subject_hash_size="$_arg_hash_size"
elif [ -z "$_arg_region" ]; then
	_subject_hash_size="16G"
else
	_subject_hash_size="1G"
fi
_subject_slots=$(normalize_hash_size "$_subject_hash_size")

_hash_size_mismatch=0
for _control_hash in $(echo $parentsString) $(echo $parentsExcludeString); do
	[ -f "$_control_hash" ] || continue          # not-yet-built generator control -> will match by construction
	_control_slots=$(read_built_hash_size "$_control_hash")
	[ -n "$_control_slots" ] || continue         # size unreadable -> don't block the run
	if [ "$_control_slots" != "$_subject_slots" ]; then
		echo "ERROR: hash size mismatch. The subject hash will be built at -s $_subject_hash_size, but a pre-built control/exclude hash was built at a different -s:" >&2
		echo "         $_control_hash" >&2
		_hash_size_mismatch=1
	fi
done
if [ "$_hash_size_mismatch" -ne 0 ]; then
	echo "Jellyfish cannot merge hashes of different sizes, so the k-mer subtraction would silently yield zero mutant k-mers." >&2
	echo "Fix: re-run with -hs/--hash_size set to the control's size (e.g. -hs 64G), or rebuild the control(s) at -s $_subject_hash_size." >&2
	_region_exit_reason="hash_size_mismatch"
	exit 100
fi
########################################################


####################__GENERATE_JHASH_FILES_FROM_JELLYFISH__#####################
_region_exit_reason="jellyfish_stage"
CONTROL_EXIT_FILES=()
SUBJECT_EXIT_FILE="$RUFUS_TMP/jelly_exit_subject.${formatted_region}.txt"

if [ $_parallel_jelly == "yes" ]
then 
	######## TODO instead of assuming 3 samples
	JThreads=$(( Threads / 3 ))
	if [ "$JThreads" -lt 3 ]
	then
    JThreads=3
	fi

	for parent in "${ParentGenerators[@]}"
	do
		ef="$RUFUS_TMP/jelly_exit_control.${formatted_region}.$(basename "$parent").txt"
		CONTROL_EXIT_FILES+=("$ef")
    	make_jelly_hash $parent $K $(echo $JThreads -2 | bc) $_arg_ParLowK $formatted_region true "$ef"  &
	done

    make_jelly_hash $ProbandGenerator $K $(echo $JThreads -2 | bc) 2 $formatted_region false "$SUBJECT_EXIT_FILE" &
    
	# wait here for all hashes to finish being made
	wait

    check_empty_hashes "$_arg_region" "${CONTROL_EXIT_FILES[@]}" "$SUBJECT_EXIT_FILE"

else
  JThreads=$Threads
  if [ "$JThreads" -lt 3 ]
    then
      JThreads=3
    fi

    for parent in "${ParentGenerators[@]}"
    do
	  ef="$RUFUS_TMP/jelly_exit_control.${formatted_region}.$(basename "$parent").txt"
	  CONTROL_EXIT_FILES+=("$ef")
      make_jelly_hash $parent $K $(echo $JThreads -2 | bc) $_arg_ParLowK $formatted_region true "$ef"
    done
    
	make_jelly_hash $ProbandGenerator $K $(echo $JThreads -2 | bc) 2 $formatted_region false "$SUBJECT_EXIT_FILE"

    check_empty_hashes "$_arg_region" "${CONTROL_EXIT_FILES[@]}" "$SUBJECT_EXIT_FILE"
fi
##############################################################################


###########################_EMPTY_JHASH_CHECK##############################
########TODO just checking file size isn't a great idea, when jellyfish fails the fields arent zero size
for parent in "${ParentGenerators[@]}"
do
    ## Check Jhash files are not empty
     if [ ! -s "$parent".Jhash ]
     then
		if [ -z $_arg_region ]
		then
			echo "ERROR: No hashes identified for $parent. This is very unlikely for a whole genome run and quality of input files should be examined. Exiting with non-zero status..."
			_region_exit_reason="no_hashes_for_parent"
			exit 100
		else
			echo "WARNING:$parent.Jhash is empty - this can happen if $parent has zero coverage for the region provided to RUFUS. Stopping run."
			_region_exit_reason="parent_empty_hash"
			exit 0
		fi
     fi
done

if [ ! -s "$ProbandGenerator".Jhash ]
then
    echo "@@@@@@@@@@@__WARNING__@@@@@@@@@@@@@"
    echo "$ProbandGenerator.Jhash  is empty"
    echo "Killing run with exit status 1"
    echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
    kill -9 $$
fi
##############################################################################



##################__GENERATE_JHASH_HISTOGRAMS__#################################
######TODO I can probably get rid of this if I just make model read either tab or space
# .Jhash.histo files are created in RunJellyForRUFUS.sh call
perl -ni -e 's/ /\t/;print' "$ProbandGenerator".Jhash.histo
for parent in "${ParentGenerators[@]}"
do
  perl -ni -e 's/ /\t/;print' "$parent".Jhash.histo
done
##############################################################################



_region_exit_reason="model_stage"
#######################__RUFUS_Model__############################################
# MODEL PHASE -- READ THIS BEFORE CHANGING THE CONDITION BELOW.
#
# As written, the `then` branch is UNREACHABLE and ModelDist never runs, in any mode.
# `_arg_min` is initialised to 5 in the defaults block near the top of this file and is
# never cleared, so `[ -z "$_arg_min" ]` is never true and the `else` branch always runs.
# The commented-out line directly below is the previous condition: originally every
# non-exome run built the model, and adding the `-z "$_arg_min"` conjunct silently
# disabled it for everyone, because of that default.
#
# Confirmed empirically across every run in resources/reg_test_files/runs/: zero
# *.7.7.dist files, fifteen *.7.7.model placeholders (written by the else branch), zero
# logs containing "Starting model phase", six containing "min was provided".
#
# Downstream effect: the else branch writes a 4-line placeholder .7.7.model and never
# produces the .7.7.dist that Overlap.shorter.sh passes to RUFUS.interpret via -mod.
# ProcessDist treats the failed open as non-fatal, so DistGlobal stays empty and the
# Bayesian genotyper is inert -- GT and FILTER come out '.' on essentially every
# SNV/indel record, FILTER=PASS is unreachable on that path, and Dist1XCutoff falls back
# to 100000, which disables the repeat filter in PickDepthSomatic.
#
# This is DOCUMENTED, NOT FIXED, and deliberately so: restoring the condition would
# enable the model for the first time and simultaneously activate several latent defects
# in RUFUS.interpret (uninitialised read in BayseanGenotyper, out-of-bounds depth-clamp
# reads, unconstrained GT ploidy in ParseGenotype). Those have to be fixed in the same
# change or output gets worse, not better. See docs/RUFUS.interpret.audit.md section 2.1.
#if [ $_arg_exome == "FALSE" ] #[	-z "$_arg_min" ]
if [ -z "$_arg_min" ]  && [ $_arg_exome == "FALSE" ]
then
	echo "No exome supplied, running whole genome mode" #echo "min not provided, building model"
	if [ -e "$ProbandGenerator.Jhash.histo.7.7.model" ]
	then
	 	echo "Skipping model phase..."
	else
		echo "Starting model phase..."
		"$RUFUSmodel" "$ProbandGenerator".Jhash.histo $K 150 $Threads > "$ProbandGenerator".Jhash.histo.7.7.out 
		for parent in "${ParentGenerators[@]}"
		do
			"$RUFUSmodel" "$parent".Jhash.histo $K 150 $Threads > "$parent".Jhash.histo.7.7.out & 
		done
		echo "Model phase complete"
	fi 

	if [ -z "$_arg_min" ]
	then
		if [ -e "$ProbandGenerator".Jhash.histo.7.7.model ]
		then
			if [ "$_arg_dev_reporting" == "TRUE" ]; then
				echo "$(grep Best\ Model "$ProbandGenerator".Jhash.histo.7.7.out)"
			fi

			MutantMinCov=$(head -2 "$ProbandGenerator".Jhash.histo.7.7.model | tail -1 )
			
			if [ "$_arg_dev_reporting" == "TRUE" ]; then
				echo "INFO: mutant min coverage from generated model is $MutantMinCov"
	 		fi
			MutantSC=$(head -4 "$ProbandGenerator".Jhash.histo.7.7.model | tail -1 )

			
			if [ "$_arg_dev_reporting" == "TRUE" ]; then
				echo "INFO: mutant SC coverage from generated model is $MutantSC"
			fi

			MaxHashDepth=$(echo "$MutantSC * 5" | bc)

			if [ "$_arg_dev_reporting" == "TRUE" ]; then
				echo "INFO: MaxHashDepth = $MaxHashDepth"
			fi
		else
			echo "ERROR Model didnt run correctly, exiting"
			return -1
		fi
	else
		echo "min coverage provided of $_arg_min, setting min kmer to that"
		MutantMinCov="$_arg_min"
	fi 
		
else 
	if [	-z "$_arg_min" ]
	then 
		echo "min coverage must be provided with an exome run"
		return -1; 
	else
# NOTE: this is a placeholder because filter step needs the arg_min value in a certain order
# The actual other values here are not(?) important or used
		echo "3" > "$ProbandGenerator".Jhash.histo.7.7.model; 
		echo "$_arg_min" >> "$ProbandGenerator".Jhash.histo.7.7.model;
		echo "3.1392e+09" >> "$ProbandGenerator".Jhash.histo.7.7.model;
		echo "1000000" >> "$ProbandGenerator".Jhash.histo.7.7.model; 
		echo "min was provided, min is $_arg_min" 
		MutantMinCov="$_arg_min" 
		#touch "$ProbandGenerator".Jhash.histo.7.7.model
	fi
fi
########################################################################################

if [ "$_arg_stop" = "jelly" ];
then
        echo "-StJ used, stopping run";
        _region_exit_reason="user_stop_after_jelly"
        exit 1;
fi
#######################################################################################

if [ -z $MutantMinCov ]; then
	echo "ERROR: No min coverage set, possible error in Model"
	_region_exit_reason="no_min_coverage_from_model"
	exit 100
fi
if [ "$MutantMinCov" -lt "2" ]
then
	echo "ERROR, model couldn't pick a sensible lower cutoff, check your subject bam file"
        _region_exit_reason="model_bad_cutoff"
        exit 1
fi
#################################__HASH_LIST_FILTER__#####################################

_region_exit_reason="hash_list_stage"
echo "Identifying unique subject kMers..."
FIFO_PREFIX="$WORK_DIR/${ProbandGenerator}.$$"
FIFO_MAIN="${FIFO_PREFIX}.temp"
FIFO_M1="${FIFO_PREFIX}.temp.mate1.fastq"
FIFO_M2="${FIFO_PREFIX}.temp.mate2.fastq"

if [ -s "$ProbandGenerator".k"$K"_c"$MutantMinCov".HashList ]
then 
    echo "$ProbandGenerator.HashList exists, skipping creation..."
else
    if [ -e "${FIFO_MAIN}" ]
    then 
    	rm  "${FIFO_MAIN}"
    fi

	# echo "About to make MAIN fifo"
    mkfifo "${FIFO_MAIN}"

	# NOTE: the modifiedJelly merge is actually the opposite of a merge
	# It intersects all of the provided Jhash files, and keeps only complements, or unique kmers
	# These unique kmers are then queried from the subject Jhash and kept only if they are from the subject (and pass min/max thresholds)
	# This was done because any attempt to simply filter the subject Jhash was prohibitively slow
    $modifiedJelly merge -o "${ProbandGenerator}.mer_counts_merged.jf" "$ProbandGenerator".Jhash $(echo $parentsString) $(echo $parentsExcludeString)  > "${FIFO_MAIN}" & 
    bash $PullSampleHashes $ProbandGenerator.Jhash "${FIFO_MAIN}" $MutantMinCov $MaxHashDepth > "$ProbandGenerator".k"$K"_c"$MutantMinCov".HashList &
	exec 5>&-
    wait
fi

########################################################################################

if [ $(head  "$ProbandGenerator".k"$K"_c"$MutantMinCov".HashList | wc -l | awk '{print $1}') -eq "0" ]; then
  if [ -z $_arg_region ]
  then
    echo "ERROR: No mutant hashes pulled from fastqs. This is very unlikely for a whole genome run and quality of input files should be examined. Exiting with non-zero status..."
    _region_exit_reason="no_mutant_hashes_wg"
    exit 100
  else
	# debug: where amy run exited
    echo "WARNING: No unique hashes identified in region $_arg_region. Stopping run."
    _region_exit_reason="no_unique_hashes"
    exit 0
  fi
fi
########################################################################################
if [ "$_arg_stop" = "hash" ];
then
        echo "-StH used, stopping run";
        _region_exit_reason="user_stop_after_hash"
        exit 1;
fi
######################__RUFUS_FILTER__##################################################
_region_exit_reason="filter_stage"
echo "Filtering unique reads containing unique kMers..."

if [ "$_pairedEnd" = "true" ]
then 
	if [ -e "$ProbandGenerator".Mutations.Mate1.fastq ]
	then
		echo "skipping filter"
	else
		if [ -z $_arg_fastqA ]
		then
			# NOTE: this is the one usually run when starting with a bam file - extracts reads from PassThroughSamCheck.stranded and puts them into mate1 and 2 respectively, which go into filter
		    if [ -e "${FIFO_M1}" ]; then 
		    	rm  "${FIFO_M1}"
		    fi
		    if [ -e "${FIFO_M2}" ]; then
	            rm  "${FIFO_M2}"
	        fi

		    mkfifo "${FIFO_M1}" "${FIFO_M2}"
		    # bash "$ProbandGenerator" | "$RDIR"/bin/PassThroughSamCheck.stranded "$WORK_DIR/$ProbandGenerator".filter.chr "${FIFO_M1}" "${FIFO_M2}" &
	    # NOTE: collate|fastq validated == PassThroughSamCheck.stranded (pairs+RC orientation, job 16684004). collate temp prefix is on WORK_DIR (lustre) — large at genome scale; redirect to node-local if needed.
	    bash "$ProbandGenerator" | samtools collate -@ "$Threads" -u -O - "$WORK_DIR/$ProbandGenerator".collate.tmp | samtools fastq -@ "$Threads" -n -1 "${FIFO_M1}" -2 "${FIFO_M2}" -s /dev/null -0 /dev/null - &
		      $RUFUSfilterFASTQ "$WORK_DIR/$ProbandGenerator".k"$K"_c"$MutantMinCov".HashList "${FIFO_M1}" "${FIFO_M2}" "$ProbandGenerator" "$K" $_filterMinQ $_arg_filterK "$(echo $Threads -2 | bc)" &
		    wait
			rm -f "$FIFO_M1" "$FIFO_M2"
		else
			echo "Running RUFUS.filter from paired FASTQ files"
			FileName=$(basename $_arg_fastqA)
			Extension="${FileName##*.}"
			if [[ $Extension == 'gz' ]]
			then
				echo "Compressed fastq files found"
				$RUFUSfilterFASTQ "$WORK_DIR/$ProbandGenerator".k"$K"_c"$MutantMinCov".HashList  <(zcat $_arg_fastqA) <(zcat $_arg_fastqB) "$ProbandGenerator" $K $_filterMinQ $_arg_filterK "$(echo $Threads -2 | bc)"

			else
				echo "Uncompressed fastq files found" 
				$RUFUSfilterFASTQ "$WORK_DIR/$ProbandGenerator".k"$K"_c"$MutantMinCov".HashList  $_arg_fastqA $_arg_fastqB "$ProbandGenerator" $K $_filterMinQ $_arg_filterK "$(echo $Threads -2 | bc)"
			fi
			wait
		fi
	fi
	if [ $(head "$ProbandGenerator".Mutations.Mate1.fastq | wc -l | awk '{print $1}') -eq "0" ]; then
		if [ -z $_arg_region ]
		then
			echo "ERROR: No reads passed the filtering step in the entire genome. This is extremely unlikely and input files should be examined."
			_region_exit_reason="no_reads_passed_filter_wg"
			exit 100
		else
			echo "No reads passed the filtering step in region $_arg_region. Stopping run."
			_region_exit_reason="no_reads_passed_filter"
			exit 0
		fi
	fi

	shortinsert="false"
	if [ -e "$ProbandGenerator".Mutations.fastq.bam ]
	then 
		echo "skipping mapping mates" 
	else
	# Sort fastq mates
        sortedMate1Fastq="$ProbandGenerator".sorted.Mutations.Mate1.fastq
        sortedMate2Fastq="$ProbandGenerator".sorted.Mutations.Mate2.fastq

        cat "$ProbandGenerator".Mutations.Mate1.fastq | paste - - - - | sort -k1 -S 8G | tr "\t" "\n" > $sortedMate1Fastq
        cat "$ProbandGenerator".Mutations.Mate2.fastq | paste - - - - | sort -k1 -S 8G | tr "\t" "\n" > $sortedMate2Fastq
		if [ $shortinsert = "false" ]
		then
			echo "skipping fastp fix"
					# NOTE: this is the one run because shortinsert is hard coded - why?
	                $bwa mem -t $Threads $_arg_ref_bwa $sortedMate1Fastq $sortedMate2Fastq | $samblaster | samtools sort -T "$ProbandGenerator".Mutations.fastq -O bam - > "$ProbandGenerator".Mutations.fastq.bam 
	                samtools index "$ProbandGenerator".Mutations.fastq.bam	
		else
			echo "using fastp fix"
			#cat "$ProbandGenerator".Mutations.Mate1.fastq "$ProbandGenerator".Mutations.Mate2.fastq > "$ProbandGenerator".Mutations.fastq
	        #$bwa mem -t $Threads $_arg_ref_bwa <( cat "$ProbandGenerator".Mutations.Mate1.fastq "$ProbandGenerator".Mutations.Mate2.fastq)  | samtools sort -T "$ProbandGenerator".Mutations.fastq -O bam - > "$ProbandGenerator".Mutations.fastq.bam
	        
			$fastp -i $sortedMate1Fastq -I $sortedMate2Fastq -m -o "$ProbandGenerator".Mutations.Mate1.fastq.fastp.fastq -O "$ProbandGenerator".Mutations.Mate2.fastq.fastp.fastq --merged_out "$ProbandGenerator".Mutations.Mate1.fastq.merged.fastq
			$bwa mem -t $Threads $_arg_ref_bwa "$ProbandGenerator".Mutations.Mate1.fastq.fastp.fastq "$ProbandGenerator".Mutations.Mate2.fastq.fastp.fastq  | $samblaster | samtools sort -T "$ProbandGenerator".Mutations.fastq -O bam - > "$ProbandGenerator".Mutations.fastq.pared.bam
			$bwa mem -t $Threads $_arg_ref_bwa "$ProbandGenerator".Mutations.Mate1.fastq.merged.fastq  | samtools sort -T "$ProbandGenerator".Mutations.fastq -O bam - > "$ProbandGenerator".Mutations.fastq.merged.bam
			samtools merge "$ProbandGenerator".Mutations.fastq.bam "$ProbandGenerator".Mutations.fastq.merged.bam "$ProbandGenerator".Mutations.fastq.pared.bam 
			
			#$bwa mem -t $Threads $_arg_ref_bwa "$ProbandGenerator".Mutations.Mate1.fastq "$ProbandGenerator".Mutations.Mate2.fastq  | $samblaster | samtools sort -T "$ProbandGenerator".Mutations.fastq -O bam - > "$ProbandGenerator".Mutations.fastq.bam
			samtools index "$ProbandGenerator".Mutations.fastq.merged.bam
			samtools index "$ProbandGenerator".Mutations.fastq.pared.bam
			samtools index "$ProbandGenerator".Mutations.fastq.bam
		fi
	fi
else
#########put se pipe here
	if [ -e "$ProbandGenerator".Mutations.fastq ]
	then
		echo "skipping filter"
	else
		if [ -z $_arg_fastqA ]
		then

		    #echo "running this one filer SE" 
	            if [ -e "${FIFO_MAIN}" ]; then
	                rm  "${FIFO_MAIN}"
	            fi
	            mkfifo "${FIFO_MAIN}"
				# NB: do NOT `exec 6<>"$FIFO_MAIN"` here. Holding a read+write fd on the FIFO keeps a
				# writer open, so RUFUS.Filter.single never sees EOF after samtools fastq finishes and
				# hangs forever (the `wait` below then never returns). The backgrounded writer and reader
				# rendezvous on the FIFO on their own.
	            # bash "$ProbandGenerator" | "$RDIR"/bin/PassThroughSamCheck.stranded.se "$WORK_DIR/$ProbandGenerator".filter.chr > "${FIFO_MAIN}" &
	            # NOTE: single-end analogue (UNTESTED for equivalence — validate before relying on it)
	            bash "$ProbandGenerator" | samtools fastq -@ "$Threads" - > "${FIFO_MAIN}" &
	              $RUFUSfilterFASTQse  "$WORK_DIR/$ProbandGenerator".k"$K"_c"$MutantMinCov".HashList "${FIFO_MAIN}" "$ProbandGenerator" "$K" $_filterMinQ $_arg_filterK "$(echo $Threads -2 | bc)" &
		    wait
			rm -f "$FIFO_MAIN"
		else
			echo "Running RUFUS.filter from single FASTQ files"
			echo "havent written this yet EXITing"
			_region_exit_reason="se_fastq_not_implemented"
			exit 1
			#########WRITE THIS##########	
			wait
		fi
	fi
	
	if [ $(head "$ProbandGenerator".Mutations.fastq | wc -l  | awk '{print $1}') -eq "0" ]; then
	  	if [ -z $_arg_region ]
		then
			echo "ERROR: No mutant hashes pulled from fastqs. This is extremely unlikely for a whole genome run and quality of input files should be examined."
			_region_exit_reason="no_mutant_fastq_reads_wg"
			exit 100
		else
			echo "WARNING: No mutant fastq reads identified in region $_arg_region. This usually means unique kmers came from reads that did not pass the sam check. Stopping RUFUS run."
			_region_exit_reason="no_mutant_fastq_reads"
			exit 0
		fi
	fi
	
	shortinsert="false"
	if [ -e "$ProbandGenerator".Mutations.fastq.bam ]
	then 
		echo "skipping mapping mates" 
	else
	  $bwa mem -t $Threads $_arg_ref_bwa "$ProbandGenerator".Mutations.fastq | $samblaster | samtools sort -T "$ProbandGenerator".Mutations.fastq -O bam - > "$ProbandGenerator".Mutations.fastq.bam
	  samtools index "$ProbandGenerator".Mutations.fastq.bam
			
	fi

fi 

########################################################################################
if [ "$_arg_saliva" = "TRUE" ]
then 
	echo "saliva sample provided, only using aligned mutant contigs"
	if [ -e  "${ProbandGenerator}".Mutations.fastq.FULL.bam ]
	then 
		echo "skipping saliva filter"
	else
		mv "${ProbandGenerator}".Mutations.fastq.bam "${ProbandGenerator}".Mutations.fastq.FULL.bam 
		samtools index "${ProbandGenerator}".Mutations.fastq.FULL.bam
		rm "${ProbandGenerator}".Mutations.fastq.bam.bai
		samtools view -F 12 -b "${ProbandGenerator}".Mutations.fastq.FULL.bam > "$ProbandGenerator".Mutations.fastq.bam
		samtools index "${ProbandGenerator}".Mutations.fastq.bam
	fi
fi



if [ $( samtools view "${ProbandGenerator}".Mutations.fastq.bam | head | wc -l | awk '{print $1}') -eq "0" ]; then
		if [ -z $_arg_region ]
		then
			echo "ERROR: All reads failed to align to the reference genome. This is extremely unlikely for a whole genome run and something likely went wrong."
			_region_exit_reason="no_reads_aligned_wg"
			exit 100
		else
       		echo "WARNING: No reads aligned to the reference for ${ProbandGenerator}.Mutations.fastq for the region $_arg_region. Stopping RUFUS run."
			_region_exit_reason="no_reads_aligned"
			exit 0
		fi
fi 
#################################################################################
if [ "$_arg_stop" = "filter" ];
then
        echo "-StF used, stopping run";
        _region_exit_reason="user_stop_after_filter"
        exit 1;
fi
###################__RUFUS_OVERLAP__#############################################
if [ -e ${ProbandGenerator}.V2.overlap.hashcount.fastq.bam.FINAL.vcf.gz ]
then
    echo "########### Skipping overlap step ###########"
else
    _region_exit_reason="overlap_stage"
    echo "########### Starting RUFUS overlap ###########"

	# Have to assign something here to maintain argument order
	if [ -z "$_arg_refhash" ]; then
		_arg_refhash="empty"
	fi

    echo " bash  $RUFUSOverlap $_arg_ref ${ProbandGenerator}.Mutations.fastq 5 $ProbandGenerator ${ProbandGenerator}.k${K}_c${MutantMinCov}.HashList $K $Threads $_MaxAlleleSize $_assemblySpeed $_arg_ref_bwa $rufus_invoc_file $_arg_refhash ${ProbandGenerator}.Jhash $parentsString $_arg_ParLowCovThreshold"
    bash  $RUFUSOverlap "$_arg_ref" "$ProbandGenerator".Mutations.fastq 5 $ProbandGenerator "$ProbandGenerator".k"$K"_c"$MutantMinCov".HashList "$K" "$Threads" "$_MaxAlleleSize" "$_assemblySpeed" "$_arg_ref_bwa" "$rufus_invoc_file" "$_arg_refhash" "$ProbandGenerator".Jhash "$parentsString" "$_arg_ParLowCovThreshold"
    #bash  $RUFUSOverlap "$_arg_ref" "$ProbandGenerator".Mutations.fastq 3 $ProbandGenerator "$ProbandGenerator".k"$K"_c"$MutantMinCov".HashList "$K" "$Threads" "$_MaxAlleleSize" "$_assemblySpeed" "$ProbandGenerator".Jhash "$parentsString" "$_arg_ref_bwa" "$_arg_refhash"
    echo "Done with RUFUS overlap"
fi
##############################################################################################


############################__RUFALU__#############################
#aluList=$RDIR/resources/primate_non-LTR_Retrotransposon.fasta
#fastaHackPath=$RDIR/bin/externals/fastahack/src/fastahack_project/bin/tools/fastahack
#jellyfishPath=$RDIR/src/externals/jellyfish-2.2.5/bin/jellyfish
#echo "$RufAlu $ProbandFileName $ProbandGenerator.V2.overlap.hashcount.fastq  $aluList $_arg_ref $jellyfishPath $(echo $ParentFileNames) "
#$RufAlu $_arg_subject $_arg_subject.generator.V2.overlap.hashcount.fastq  $aluList $_arg_ref $fastaHackPath $jellyfishPath  $(echo $ParentFileNames)
########################################################################

_region_exit_reason="vcf_processing_stage"
intermed_vcf="${WORK_DIR}/${ProbandGenerator}.V2.overlap.hashcount.fastq.bam.vcf"

if [[ -s "$intermed_vcf" ]]; then
	count=$(bcftools view -H "$intermed_vcf" | wc -l)
	if [ "$count" -eq 0 ]; then
	  	echo "Intermediate vcf contains no variants, indicating no variants found for this region." >&2
		_region_exit_reason="no_interpret_passing_vars_e"
		exit 0
	fi
  	# safe to proceed (file exists, non-empty, has variants)
else
	echo "Intermediate vcf not present, indicating no variants found for this region." >&2
	_region_exit_reason="no_interpret_passing_vars_dne"
	exit 0
fi

# Sanitize intermediate VCF: remove malformed records from RUFUS.Interpret output
# Keeps header lines as-is. For data lines, requires:
#   - at least 8 tab-delimited fields (CHROM POS ID REF ALT QUAL FILTER INFO)
#   - POS (col 2) is a positive integer
#   - REF (col 4) is non-empty and contains only valid bases (ACGTN)
#   - ALT (col 5) is non-empty, not just ".", and contains only valid VCF ALT characters
#   - If INFO (col 8) contains END=<n>, then END >= POS (prevents tabix "end < begin" error)
sanitized_vcf="${intermed_vcf}.sanitized.vcf"
awk -F'\t' '
/^#/ { print; next }
{
	if (NF < 8) next
	if ($2 !~ /^[0-9]+$/ || $2+0 < 1) next
	if ($4 == "" || $4 == "." || $4 !~ /^[ACGTNacgtn]+$/) next
	if ($5 == "" || $5 == ".") next
	if ($5 !~ /^[ACGTNacgtn.,*<>\[\]0-9:]+$/) next
	# Check END >= POS if END tag is present in INFO field
	info = $8
	if (match(info, /END=[0-9]+/)) {
		end_val = substr(info, RSTART+4, RLENGTH-4) + 0
		if (end_val < $2+0) next
	}
	print
}
' "$intermed_vcf" > "$sanitized_vcf"

sanitized_count=$(grep -vc "^#" "$sanitized_vcf" || true)
original_count=$(grep -vc "^#" "$intermed_vcf" || true)
removed_count=$((original_count - sanitized_count))
if [ "$removed_count" -gt 0 ]; then
	echo "WARNING: Removed $removed_count malformed VCF record(s) from RUFUS.Interpret output ($sanitized_count of $original_count records kept)." >&2
fi
if [ "$sanitized_count" -eq 0 ]; then
	echo "No valid VCF records remain after sanitization for this region." >&2
	_region_exit_reason="no_records_after_sanitization"
	exit 0
fi
mv "$sanitized_vcf" "$intermed_vcf"

# Trim off generator postfix
DEDUPED_VCF="$WORK_DIR/deduped.${formatted_region}.vcf"

# TODO: do I really need this? can I just sort?
grep "^#" "$intermed_vcf" > $WORK_DIR/Intermediates/$ProbandGenerator.V2.overlap.hashcount.fastq.bam.sorted.vcf
grep -v "^#" "$intermed_vcf" | sort -k1,1V -k2,2n >> $WORK_DIR/Intermediates/$ProbandGenerator.V2.overlap.hashcount.fastq.bam.sorted.vcf

echo "arg_mosaic = $_arg_mosaic"
if [ "$_arg_mosaic" == "TRUE" ]
then
	echo "including mosaic"
	bash $RDIR/scripts/VilterAutosomeOnly $WORK_DIR/Intermediates/$ProbandGenerator.V2.overlap.hashcount.fastq.bam.sorted.vcf | perl $RDIR/scripts/ColapsDuplicateCalls.stream.pl > $DEDUPED_VCF
else
	echo "excluding mosaic"
	bash $RDIR/scripts/VilterAutosomeOnly.withoutMosaic $WORK_DIR/Intermediates/$ProbandGenerator.V2.overlap.hashcount.fastq.bam.sorted.vcf | perl $RDIR/scripts/ColapsDuplicateCalls.stream.pl > $DEDUPED_VCF
fi

bgzip -f "$DEDUPED_VCF"
# Index with tabix, iteratively removing records that cause indexing failures
# This catches any malformed records that slip past the awk sanitizer
tabix_max_retries=50
tabix_attempt=0
while true; do
	tabix_stderr=$(tabix -C "$DEDUPED_VCF.gz" 2>&1) && break

	tabix_attempt=$((tabix_attempt + 1))
	if [ "$tabix_attempt" -ge "$tabix_max_retries" ]; then
		echo "ERROR: tabix failed after removing $tabix_attempt malformed record(s). Giving up." >&2
		echo "Last tabix error: $tabix_stderr" >&2
		_region_exit_reason="tabix_max_retries"
		exit 100
	fi

	# Parse the 1-based sequence number from: "Invalid record on sequence #N"
	bad_seq=$(echo "$tabix_stderr" | grep -oP 'sequence #\K[0-9]+' | head -1)
	if [ -z "$bad_seq" ]; then
		echo "ERROR: tabix failed with unexpected error: $tabix_stderr" >&2
		_region_exit_reason="tabix_unexpected_error"
		exit 100
	fi

	echo "WARNING: tabix indexing failed on data record #${bad_seq}, removing it and retrying (attempt $tabix_attempt)." >&2
	echo "  tabix error: $tabix_stderr" >&2

	# Decompress, remove the offending data line, recompress
	tmp_fix_vcf="${DEDUPED_VCF}.tabixfix.vcf"
	zcat "$DEDUPED_VCF.gz" | awk -v bad="$bad_seq" '
		/^#/ { print; next }
		{ data_line++; if (data_line != bad) print }
	' > "$tmp_fix_vcf"
	bgzip -f "$tmp_fix_vcf"
	mv "$tmp_fix_vcf.gz" "$DEDUPED_VCF.gz"
done

# Update reference alleles
REF_VCF="$WORK_DIR/ref.${formatted_region}.vcf"
bcftools +fill-from-fasta "$DEDUPED_VCF.gz" -- -c REF -f "$_arg_ref" > "$REF_VCF"

# Get rid of break-ends
TYPE_VCF="$WORK_DIR/snv_indel.${formatted_region}.vcf"
bcftools view -e "TYPE='bnd'" "$REF_VCF" > "$TYPE_VCF"

# Check for empty gt field
GX_VCF="$WORK_DIR/gx.${formatted_region}.vcf"
bash $RDIR/post_process/remove_no_genotype.sh "$TYPE_VCF" > "$GX_VCF"
bgzip -f "$GX_VCF"
bcftools index -f "$GX_VCF.gz"

# Trim calls to region. In whole-genome mode _arg_region is empty; `bcftools view -r ""` segfaults,
# and there is nothing to trim to, so pass the calls through unchanged.
TRIMMED_VCF="$WORK_DIR/trimed.${formatted_region}.vcf.gz"
if [ -n "$_arg_region" ]; then
	bcftools view -r "$_arg_region" "$GX_VCF.gz" -Oz -o "$TRIMMED_VCF"
else
	cp "$GX_VCF.gz" "$TRIMMED_VCF"
fi
bcftools index "$TRIMMED_VCF"

NO_CO_VCF="$WORK_DIR/no_coinheriteds.vcf.gz"
# remove_coinheriteds pileups each control at the variant sites, so it needs an alignable BAM/CRAM.
# A control given as a pre-built hash (a .generator stub) or fastq has no BAM to pile up (bwa would
# align an empty file -> mpileup fails on the empty bam). Collect only the BAM/CRAM controls and run
# the filter over those; skip entirely if none -- the HashList subtraction has already removed those
# controls' k-mers, so the co-inherited pileup is a secondary check with nothing to pile up.
_bamcram_controls=()
for _ctrl in "${Parents[@]}"; do
	case "$_ctrl" in
		*.bam|*.cram) _bamcram_controls+=("$_ctrl") ;;
	esac
done
if [ ${#_bamcram_controls[@]} -ne "0" ]; then
	bash ${RDIR}/post_process/remove_coinheriteds.sh -t $_arg_threads -r "$formatted_region" -f "$_arg_ref" -i "$TRIMMED_VCF" -o "$NO_CO_VCF" -w "1000" -c "$(IFS=','; echo "${_bamcram_controls[*]}")"
else
	[ ${#_arg_controls[@]} -ne "0" ] && echo "Skipping remove_coinheriteds: no BAM/CRAM control to pile up (controls are hash/generator/fastq); HashList subtraction already handled them." >&2
	mv "$TRIMMED_VCF" "$NO_CO_VCF"
fi

# Left align & atomize
ATOM_VCF="$WORK_DIR/atomed.${formatted_region}.vcf"
bcftools norm -m- -f "$_arg_ref" "$NO_CO_VCF" -Ou | bcftools norm -a -Oz -o "$ATOM_VCF"

# Add HD_AF field
ATOM_VCF_BASENAME=$(basename "$ATOM_VCF")
HDAF_VCF="$WORK_DIR/hd_af.$ATOM_VCF_BASENAME"
SUBJECT_SAMPLE_NAME=$(bcftools view -h "$ATOM_VCF" | tail -n 1 | awk -F'\t' '{ print $10 }')
bash ${RDIR}/post_process/add_hd_med.add_hd_af.sh "$ATOM_VCF" "$SUBJECT_SAMPLE_NAME" "$formatted_region"

# Sort
SORTED_VCF="$WORK_DIR/sorted.${formatted_region}.vcf.gz"
bcftools sort "$HDAF_VCF" -Oz -o "$SORTED_VCF"

# Rename final vcf and zip/index
PREFINAL_VCF="$WORK_DIR/temp.RUFUS.Final.${ProbandFileName}${region_postfix}.vcf.gz"
mv "$SORTED_VCF" "$PREFINAL_VCF"
bcftools index "$PREFINAL_VCF"

FINAL_BASENAME="$(basename "$PREFINAL_VCF")"

cp "$PREFINAL_VCF" "$WORK_ROOT/$FINAL_BASENAME"
cp "$PREFINAL_VCF.csi" "$WORK_ROOT/$FINAL_BASENAME.csi"

end_time=$(date +"%s")
time_delta=$(( $end_time - $start_time ))
hours=$(( time_delta / 3600 ))
minutes=$(( (time_delta % 3600) / 60 ))
seconds=$(( time_delta % 60 ))
printf "RUFUS call stage completed in: %02d:%02d:%02d\n" $hours $minutes $seconds
_region_exit_reason="success"
exit 0
# ] <-- needed because of Argbash
