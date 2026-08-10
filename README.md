RUFUS
=====

K-mer based variant detection. v1.2.0.

Developed by Stephanie Georges, MS\
Based on the thesis project of Andrew Farrell, PhD\
Marth Lab, University of Utah\
*Publication in preparation*

For questions and feature requests, please contact [stephanie.georges@genetics.utah.edu](stephanie.georges@genetics.utah.edu)

## RUFUS Overview

RUFUS is a reference-bias-free, K-mer based variant detection algorithm, for short-read DNA sequence data. RUFUS is intended to run on a high performance computing (HPC) cluster with Apptainer (formerly Singularity) or Docker installed. At a high level, you'll need to download the pre-built container (detailed below) and either use the provided setup script to generate a SLURM script that runs RUFUS or manually create an execution script directly. 

RUFUS calls variants in a single subject (but may be split across multiple files) against one or more control samples. Input files may be FASTQ, CRAM, BAM, or a RUFUS generator file. The reference genome must be in FASTA format, and must be indexed by BWA. If the BWA indexes are not detected in the same directory as the reference genome, RUFUS will create them.

RUFUS has two stages: a variant calling stage, and a post-processing stage. Separation of the stages is necessary because the calling stage may be run in a windowed fashion, requiring multiple parallel RUFUS jobs over all of the windows. The combination stage must wait to proceed until all calling jobs are complete. Algorithmic runtime increases roughly linearly with sample coverage. Generally with whole-genome mode, a 300x sample run will take 1 day. Windowed mode completes significantly faster given sufficient resources (see below for recommendations).


## Running RUFUS

### Obtaining the RUFUS Image

The pre-built RUFUS container is published to two places: Docker Hub, and
[Zenodo](https://doi.org/10.5281/zenodo.13694210) for archival and citation. Either route gives you
the same image.

**From Docker Hub (recommended on HPC).** `apptainer` builds the SIF for you — no `sudo`, no
manual `.def`:
```bash
apptainer pull rufus.sif docker://stefinfection/rufus:latest
apptainer pull rufus.sif docker://stefinfection/rufus:v1.2.0
```
`:latest` always points at the most recent release; pin a specific version instead for a
reproducible analysis. Published versions are listed at
https://hub.docker.com/r/stefinfection/rufus/tags.

**From Zenodo.** Zenodo
does not expose a fixed download path for "latest", so ask its API which file to fetch rather than
building a URL by hand:
```bash
CONCEPT_DOI=13694210   # the concept DOI suffix, 10.5281/zenodo.13694210
URL=$(curl -fsSL --retry 3 --retry-delay 5 "https://zenodo.org/api/records/${CONCEPT_DOI}" \
      | python3 -c "import sys,json;print(next(f['links']['self'] for f in json.load(sys.stdin)['files'] if f['key'].endswith('.sif')))")
[ -n "$URL" ] || { echo "could not resolve the latest RUFUS SIF from Zenodo" >&2; exit 1; }
curl -fL --retry 3 --retry-delay 5 -o rufus.sif "$URL"
```

To browse releases or download via the browsser instead, see https://zenodo.org/records/13694210/latest.

### Input Data

RUFUS requires the following data to run:
1) A subject sample in FASTQ/BAM/CRAM/generator format. BAM/CRAM may be unaligned when using whole genome mode; FASTQ is whole-genome only and cannot be combined with `-R/--region`. If the sample is split across several files, pass `-s` once per file — they are treated as one sample, not as separate subjects.
2) At least one of: one or more control samples (`-c`, same formats as the subject), or an exclude hash (`-e`). Multiple distinct controls are supported, e.g. mother and father for a trio. Supplying only `-e` — typically pre-built 1000G/control hashes — is single-sample mode.

>BAM, CRAM and generator inputs may be freely mixed across `-s` and `-c` — they all feed the same read stream internally. FASTQ must stand alone: if any input is a FASTQ, all of them must be. RUFUS filters reads from the FASTQ mate files directly in that case, so reads from any BAM/CRAM/generator alongside them would be counted but never filtered, quietly costing calls.

>A generator file is a shell script that writes SAM to stdout (e.g. a single line `samtools view -h -F 3328 /path/sample.bam`); RUFUS runs it to obtain reads. Generators are whole-genome only — `-R/--region` is not applied to them, so the SLURM launcher (helper file, details included below) rejects them in windowed (`-w`) mode. This input style is generally used for special cases and not applicable to most users.

3) A reference fasta file (this must be indexed by BWA) - for use in reporting the called variants. *It's recommended to provide the BWA indexes in the same directory as the reference to save time creating them during the RUFUS run.*\
\
To create the BWA indexes, run the following commands:
   ```
   bwa index -a bwtsw {REFERENCE.fa}
   samtools faidx {REFERENCE.fa}
   ```

All input files are specified by their full paths. The necessary host directories are automatically bind-mounted into the Apptainer container.

### Output Data

RUFUS will, by default, output the following files *in the current working directory*:
1) A VCF file containing the called variants named `RUFUS.Final...vcf.gz`
2) A supplemental directory with:
    * A pre-filtered VCF file
    * A BAM file containing the raw reads containing the mutant kmers
    * A BAM file containing the assembled contigs from the raw reads containing the mutant kmers
    * A hash table containing the unique subject kmers and their counts
3) In windowed (region) mode, a `region_status.log` file summarizing the outcome of each region (variants called, no variants found with reason, or error with exit code)

 
### The Two Stages of RUFUS

RUFUS has two execution stages:
1) The calling stage, invoked by the following
```
apptainer exec {PATH_TO_RUFUS_CONTAINER}/rufus.sif bash /opt/RUFUS/runRufus.sh [-s|--subject <arg>] [-r|--ref <arg>] [-t|--threads <arg>] [-k|--kmersize <arg>] [-m|--min <arg>] [-h|--help] [-c|<controls-1>] ... [-c|<controls-n>] ...OPTIONS
```
With the following usage:
```
Required Arguments:
    -s,--subject:   bam/cram/fastq/generator file(s) containing the subject of interest. Use
                    multiple times ONLY for split files of the same sample (e.g.
                    -s subject.part1.bam -s subject.part2.bam); they are combined into one
                    subject, not called separately
    -r,--ref:       file path to the desired reference file
    -t,--threads:   number of threads to use (min 3)

Optional Arguments:
    -c,--controls:  bam/cram/fastq/generator file(s) for the sequence data of a control sample
                    (can be used multiple times for distinct controls, e.g.
                    -c mother.bam -c father.bam). Optional ONLY if -e/--exclude is supplied
                    instead: RUFUS requires at least one control or exclude source and will exit
                    if given neither. Supplying only -e is single-sample mode
    -k,--kmersize:  length of k-mer to use (defaults to 25)
    -m,--min:       overwrites the minimum k-mer depth count to call variant (defaults to 5)
    -e,--exclude:   Jhash file of kmers to exclude from mutation list (can be used multiple
                    times, e.g. -e Jhash1 -e Jhash2)
    -f,--refhash:   Jhash file containing reference hashList
    -R,--region:    genomic region to call variants on (e.g. chr1:1-1000000); used in windowed
                    mode
    -h,--help:      print help
```

2) The post-processing stage, invoked by the following
```
apptainer exec {PATH_TO_RUFUS_CONTAINER}/rufus.sif bash /opt/RUFUS/post_process/post_process.sh -s <subject> -w <window_size>
```
With the following usage:
```
Required Arguments:
    -s subject_file  the name of the subject file; must be the same as that supplied to the
                     RUFUS run
    -w window_size   the size of the window used in the RUFUS run (0 for whole genome mode)

Optional Arguments:
    -h help          print help message
```

In windowed mode, the post-processing stage will print a region status summary showing how many regions called variants, how many had no variants (with breakdown by reason), and how many encountered errors. The final VCF header will include a `##RUFUS_runMode` line indicating region mode was used.


## Using the SLURM Helper Script & Executing the SLURM Batch Scripts

The SLURM helper script automatically creates the two SLURM batch scripts necessary to run RUFUS on a SLURM-managed HPC cluster, as well as a bash script to execute them. To use:
1) Execute the helper script (see full usage options below):
``` 
apptainer exec {PATH_TO_RUFUS_CONTAINER}/rufus.sif bash /opt/RUFUS/singularity/setup_slurm.sh [-s subject] [-c control1,control2,control3...] [-b genome_build] [-a slurm_account] [-p slurm_partition] ...OPTIONS
```

2) Then execute the generated bash script:
```
bash launch_rufus.sh
```

The full usage options for the helper script are as follows:
```
Required Arguments:
    -s subject                  full path to the subject sample BAM/CRAM
    -b genome_build             the desired genome build; currently only supports GRCh38
    -r reference                full path to the reference file matching the genome build
    -a slurm_account            the account for the slurm job
    -p slurm_partition          the partition for the slurm job
    -l slurm_job_array_limit    the maximum amount of jobs slurm allows in an array

Optional Arguments:
    -c control(s)               a single control or comma-delimited array of multiple controls
                                (full paths). Omit for single-sample mode, in which case you must
                                supply a hash source instead -- see "Single-sample mode" below
    -m kmer_depth_cutoff        the amount of kMers that must overlap the variant to be included
                                in the final call set
    -w window_size              the size of the windows to run RUFUS on, in units of kilobases
                                (KB); allowed range between 500-5000; defaults to a single run of
                                the entire genome if not provided
    -f reference_hash           Jhash file containing reference kMer hash list
    -x exclude_hash             single or comma-delimited list of Jhash file(s) containing kMers
                                to exclude (static, same for all regions)
    -K kg1_hash_dir             full path to directory of per-region KG1 Jhash files (files named
                                *{region}*.Jhash)
    -G kg1_version              KG1 hash version to download from S3 (e.g. v3.0)
    -D ctrl_hash_dir            full path to directory of per-region control Jhash files (files
                                named *{region}*.Jhash)
    -V ctrl_version             control hash version to download from S3 (e.g. v1.0)
    -M memory_per_call          how much memory to allot to the rufus calling stage job (e.g.
                                150G or 20G)
    -C cpus_per_call            how many cpus to allot to each rufus calling stage job; defaults
                                to 40 for whole genome, 12 for 1MB windows. Must be greater than
                                the RUFUS thread count (-z), or setup will exit with an error
    -y path_to_rufus_container  if not provided, will look in the current directory for rufus.sif
    -z rufus_threads            number of threads provided to RUFUS; defaults to 36 for the
                                entire genome, 10 for 1MB windows
    -e email                    the email address to notify with slurm updates
    -q slurm_job_queue_limit    the maximum amount of jobs able to be run at once; defaults to 20
    -t slurm_time_limit         the maximum amount of time to let the slurm job run; defaults to
                                7 days for a full run, or one hour per window (DD-HH:MM:SS)
    -h help                     print usage
```
\
*Notes on SLURM arguments*:\
This script utilizes SLURM arrays to batch RUFUS call runs and thus requires the SLURM job array limit to comply with user settings. To find your SLURM job array limit:
```
scontrol show config | grep "MaxArraySize"
```

To maximize parallelism, filling in the slurm job queue limit (-q) is recommended. You can find your limit by typing:
```
scontrol show config | grep "default_queue_depth"
```

### Example Invocations of the helper script

Basic windowed mode:
```
apptainer exec /home/my_container_path/rufus.sif bash /opt/RUFUS/singularity/setup_slurm.sh -s /home/subjects/subject.bam -c /home/controls/control_a.bam,/home/controls/control_b.bam -r /refs/GRCh38_reference.fa -a my-slurm-account -p my-slurm-partition -w 1000 -t "00:30:00" -m 5 -l 20 -z 36 -e "my_email@utah.edu"
```

With local per-region hash directories:
```
apptainer exec /home/my_container_path/rufus.sif bash /opt/RUFUS/singularity/setup_slurm.sh -s /home/subjects/subject.bam -c /home/controls/control_a.bam -r /refs/GRCh38_reference.fa -a my-slurm-account -p my-slurm-partition -w 1000 -l 1000 -K /data/kg1_hashes/v3.0/ -D /data/ctrl_hashes/v1.0/
```

With S3-downloaded hashes (downloaded at run time):
```
apptainer exec /home/my_container_path/rufus.sif bash /opt/RUFUS/singularity/setup_slurm.sh -s /home/subjects/subject.bam -c /home/controls/control_a.bam -r /refs/GRCh38_reference.fa -a my-slurm-account -p my-slurm-partition -w 1000 -l 1000 -G v3.0 -V v1.0
```

### Single-sample mode

RUFUS does not require a matched control. If you have no control sample, omit `-c` and supply a
pre-built k-mer hash source instead — RUFUS subtracts against those hashes rather than against a
control you sequenced. Any of `-x`, `-K`/`-G`, or `-D`/`-V` satisfies this; they are passed through
to the calling stage as `-e/--exclude` arguments.

You must supply at least one of a control or a hash source. Providing neither is rejected at the
start of the calling stage.

Single-sample, windowed, with S3-downloaded 1000 Genomes and control hashes:
```
apptainer exec /home/my_container_path/rufus.sif bash /opt/RUFUS/singularity/setup_slurm.sh -s /home/subjects/subject.bam -r /refs/GRCh38_reference.fa -a my-slurm-account -p my-slurm-partition -b GRCh38 -w 1000 -l 1000 -G v3.0 -V v1.0
```

The same run against hash directories you already hold locally:
```
apptainer exec /home/my_container_path/rufus.sif bash /opt/RUFUS/singularity/setup_slurm.sh -s /home/subjects/subject.bam -r /refs/GRCh38_reference.fa -a my-slurm-account -p my-slurm-partition -b GRCh38 -w 1000 -l 1000 -K /data/kg1_hashes/v3.0/ -D /data/ctrl_hashes/v1.0/
```

*Note*: `-K`/`-G` (KG1 hashes) and `-D`/`-V` (control hashes) are mutually exclusive per type. You may mix local and S3 across types (e.g., `-K /local/kg1/ -V v1.0`). Hash files in local directories must be named with the region string (e.g., `*chr1_1_1000000*.Jhash` for region `chr1:1-1000000`, or `*wg*.Jhash` for whole-genome mode).

=======
