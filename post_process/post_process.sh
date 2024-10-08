#!/bin/bash

usage() {
	echo "Usage: $0 [-w window_size] [-r reference] [-s subject] [-c control1,control2,control3...] [-d source_dir] [-h]"
	echo "Options:"
	echo " -w window_size	Required: The size of the window used in the RUFUS run"
	echo " -r reference	Required: The reference used in the RUFUS run"
	echo " -c controls	Required: The control bam files used in the RUFUS run"
	echo " -s subject_file	Required: The name of the subject file: must be the same as that supplied to the RUFUS run"
	echo " -d source_dir	Required: The source directory where the RUFUS vcf(s) are located"
	echo " -h help	Print help message"
	exit 1
}

# initialize vars
CONTROLS=()
WINDOW_SIZE=0
REFERENCE=""
SUBJECT_FILE=""
SOURCE_DIR="/mnt"

# parse command line arguments
while getopts ":w:r:c:s:d:h" option; do 
	case $option in 
		h) usage;;
		w) WINDOW_SIZE=$OPTARG;;
		r) REFERENCE=$OPTARG;;
		s) SUBJECT_FILE=$OPTARG;;
		d) SOURCE_DIR=$OPTARG;;
		c) IFS=',' read -r -a CONTROLS <<< "$OPTARG";;
		\?) echo "Invalid option: -$OPTARG" >&2
		    usage;;	
		*) echo "Option -$OPTARG requires an argument" >&2
		   usage;;
	esac
done
shift $((OPTIND-1))

# check for mandatory command line arguments
if [[ -z "$WINDOW_SIZE" ]]; then
	    echo "Error: Missing required option -w (window size)" >&2
	        usage
fi

if [[ -z "$REFERENCE" ]]; then
	    echo "Error: Missing required option -r (reference)" >&2
	        usage
fi

if [[ -z "$SOURCE_DIR" ]]; then
	echo "Error: Missing required option -d (source directory for RUFUS vcf(s))" >&2
	        usage
fi

if [ ${#CONTROLS[@]} -eq 0 ]; then
	    echo "Error: Must supply at least one control bam" >&2
	        usage
fi

cd $SOURCE_DIR
echo "RUFUS post-process version C.0.1"
date
start_time=$(date +"%s")

POST_PROCESS_DIR=/opt/RUFUS/post_process/
TEMP_FINAL_VCF="temp.RUFUS.Final.${SUBJECT_FILE}.combined.vcf.gz"
TEMP_PREFILTERED_VCF="temp.RUFUS.Prefiltered.${SUBJECT_FILE}.combined.vcf.gz"

if [ "$WINDOW_SIZE" != "0" ]; then
	IFS=$'\t'
	TAB_DELIM_CONTROL_STRING="${CONTROLS[*]}"
	echo "Windowed run performed, trimming and combining region vcfs..."
	bash ${POST_PROCESS_DIR}trim_and_combine.sh $SUBJECT_FILE $TAB_DELIM_CONTROL_STRING $WINDOW_SIZE
else 
	# Slight name change if not doing a windowed run
	TEMP_FINAL_VCF="temp.RUFUS.Final.${SUBJECT_FILE}.vcf.gz"
	TEMP_PREFILTERED_VCF="${SUPP_DIR}temp.RUFUS.Prefiltered.${SUBJECT_FILE}.vcf.gz"
fi

# TODO: check to see if we had any variants in final, and if not, stop and report
VARS_REPORTED=$(bcftools view -H $TEMP_FINAL_VCF | wc -l)
if [ "$VARS_REPORTED" = "0" ]; then
	echo "RUFUS did not find any variants for the provided parameters. Please adjust and try again."
	rm /mnt/${SUBJECT_FILE}*.generator*
	for control in "${CONTROLS[@]}"; do
    	rm /mnt/${control}*.generator*
	done
	rm -r /mnt/Intermediates
	rm -r /mnt/TempOverlap	
	rm -r /mnt/rufus.cmd
	rm /mnt/temp*.vcf*

	echo "RUFUS did not find any variants for the provided parameters. Please adjust and try again." > fail.out
	exit 0
fi

# Check for empty lines
echo "Checking vcf formatting..."
bash ${POST_PROCESS_DIR}remove_no_genotype.sh $TEMP_FINAL_VCF "final_no_gx.vcf"
bash ${POST_PROCESS_DIR}remove_no_genotype.sh $TEMP_PREFILTERED_VCF "prefiltered_no_gx.vcf"
rm $TEMP_FINAL_VCF
rm $TEMP_PREFILTERED_VCF
mv "final_no_gx.vcf.gz" $TEMP_FINAL_VCF
mv "prefiltered_no_gx.vcf.gz" $TEMP_PREFILTERED_VCF

# Sort
echo "Sorting..."
bcftools sort $TEMP_FINAL_VCF | bgzip > "sorted.${TEMP_FINAL_VCF}"
# TODO: when fix formatting on prefiltered vcf, comment two lines below back in
#bcftools sort $TEMP_PREFILTERED_VCF | bgzip > "sorted.${TEMP_PREFILTERED_VCF}"

rm $TEMP_FINAL_VCF
#rm $TEMP_PREFILTERED_VCF
bcftools index "sorted.$TEMP_FINAL_VCF"

# Remove coinheriteds
echo "Removing coinheriteds..."
IFS=$','
CONTROL_STRING="${CONTROLS[*]}"
COINHERITED_REMOVED_VCF="coinherited_removed.vcf.gz"
bash ${POST_PROCESS_DIR}remove_coinheriteds.sh "$REFERENCE" "sorted.${TEMP_FINAL_VCF}" "$COINHERITED_REMOVED_VCF" "$SOURCE_DIR" "$CONTROL_STRING"

# Add HD_AF field
echo "Adding kmer-based allele frequencies..." 
AF_ADDED_VCF="hd_af.${COINHERITED_REMOVED_VCF}"
SUBJECT_SAMPLE_NAME=$(bcftools view -h $COINHERITED_REMOVED_VCF | tail -n 1 | awk -F'\t' '{ print $10 }')
bash ${POST_PROCESS_DIR}add_hd_med.add_hd_af.sh "$COINHERITED_REMOVED_VCF" "$SUBJECT_SAMPLE_NAME"
bcftools index $AF_ADDED_VCF 

# Compose final vcfs
SUBJECT_STRING=$(basename $SUBJECT_FILE)
FINAL_VCF="RUFUS.Final.${SUBJECT_STRING}.combined.vcf"
PREFILTERED_VCF="RUFUS.Prefiltered.${SUBJECT_STRING}.combined.vcf"

# Inject RUFUS command into header
echo "Composing final vcfs..."
bcftools view -h $AF_ADDED_VCF | head -n -1 > $FINAL_VCF
cat /mnt/rufus.cmd >> $FINAL_VCF
bcftools view -h $AF_ADDED_VCF | tail -n 1 >> $FINAL_VCF
bcftools view -H $AF_ADDED_VCF >> $FINAL_VCF
bgzip $FINAL_VCF
bcftools index "$FINAL_VCF.gz"

#TODO: Comment back in after prefiltered vcf cleaned up
#bcftools view -h $TEMP_PREFILTERED_VCF | head -n -1 > $PREFILTERED_VCF
#cat /mnt/rufus.cmd >> $PREFILTERED_VCF
#bcftools view -h $TEMP_PREFILTERED_VCF | tail -n 1 >> $PREFILTERED_VCF
#bcftools view -H $TEMP_PREFILTERED_VCF >> $PREFILTERED_VCF
#bgzip $PREFILTERED_VCF
#bcftools index "$PREFILTERED_VCF.gz"
#mv "$PREFILTERED_VCF.gz"* rufus_supplementals/

# Only need to move and rename if did a windowed run
if [ "$WINDOW_SIZE" != "0" ]; then
	mv $TEMP_PREFILTERED_VCF prefiltered.vcf.gz
	mv $TEMP_PREFILTERED_VCF.tbi prefiltered.vcf.gz.tbi
	mv prefiltered.vcf.gz* rufus_supplementals/
fi


# TODO: Separate SVs and SNV/Indels
#echo "Separating snvs/indels and SVs..."

# Cleanup
echo "Cleaning up intermediate post-processing files..."
#rm $TEMP_PREFILTERED_VCF*
rm $TEMP_FINAL_VCF*
#rm "sorted.$TEMP_PREFILTERED_VCF"*
rm "sorted.$TEMP_FINAL_VCF"*
rm $COINHERITED_REMOVED_VCF*
rm "normed.sorted.$TEMP_FINAL_VCF"*
rm -r "/mnt/Intermediates"
rm -r "/mnt/TempOverlap"
rm "/mnt/rufus.cmd"
rm "$AF_ADDED_VCF"*

# Combining supplementals
SUPPLEMENTAL_DIR=/mnt/rufus_supplementals/
# TODO: only do this if not reporting in developer mode
ls ${SUPPLEMENTAL_DIR}*generator.V2.overlap.hashcount.fastq.bam | xargs samtools merge ${SUPPLEMENTAL_DIR}unique_contigs.bam
ls ${SUPPLEMENTAL_DIR}*generator.Mutations.fastq.bam | xargs samtools merge ${SUPPLEMENTAL_DIR}unique_reads.bam
samtools sort ${SUPPLEMENTAL_DIR}unique_contigs.bam -o ${SUPPLEMENTAL_DIR}unique_contigs.sorted.bam
samtools sort ${SUPPLEMENTAL_DIR}unique_reads.bam -o ${SUPPLEMENTAL_DIR}unique_reads.sorted.bam
rm ${SUPPLEMENTAL_DIR}unique_contigs.bam
rm ${SUPPLEMENTAL_DIR}unique_reads.bam
rm ${SUPPLEMENTAL_DIR}*generator.V2.overlap.hashcount.fastq.bam*
rm ${SUPPLEMENTAL_DIR}*generator.Mutations.fastq.bam*

cat ${SUPPLEMENTAL_DIR}*.HashList > ${SUPPLEMENTAL_DIR}unique_kmer_counts.txt
rm ${SUPPLEMENTAL_DIR}*.HashList

# TODO: hack for now to get rid of any lingering intermediate files 
# TODO: this actually needs to be fixed with a graceful handling of no variants for certain windows
rm /mnt/${SUBJECT_FILE}*.generator*
for control in "${CONTROLS[@]}"; do
	rm /mnt/${control}*.generator*
done

echo "Post-processing complete."
end_time=$(date +"%s")
time_delta=$(( $end_time - $start_time ))
hours=$(( time_delta / 3600 ))
minutes=$(( (time_delta % 3600) / 60 ))
seconds=$(( time_delta % 60 ))
printf "RUFUS call stage completed in: %02d:%02d:%02d\n" $hours $minutes $seconds
