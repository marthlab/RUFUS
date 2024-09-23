# This script removes inherited variants called by rufus which lie on the same contig as
# a somatic variant. It accomplishes this task by performing a pileup and variant call
# in the control bam for each of the rufus-identified variants, then taking the complement
# of the control set in an intersection.

# parse command line arguments
REFERENCE_FILE=$1  # Reference file (.fa) - must be the same as RUFUS run that generated RUFUS_VCF
RUFUS_VCF=$2       # The 'FINAL' vcf generated by rufus, to be isec'd
OUT_VCF=$3
SRC_DIR=$4
ARG_LIST=("$@")
CONTROL_BAM_LIST=("${ARG_LIST[@]:4}") # Remaining args, all control bams

cd $SRC_DIR

# static vars
CONTROL_ALIGNED="temp_aligned.bam"
CONTROL_VCF="isec_control.vcf.gz"
NORMED_VCF="normed.${RUFUS_VCF}"
BWA=/opt/RUFUS/bin/externals/bwa/src/bwa_project/bwa
PILEUP_SCRIPT="/opt/RUFUS/post_process/single_pileup.sh"

# make intersection directory
ISEC_OUT_DIR="temp_isecs"
mkdir -p $ISEC_OUT_DIR

#format final rufus vcf for intersections
vt normalize -n $RUFUS_VCF -r $REFERENCE_FILE | vt decompose_blocksub - | bgzip > $NORMED_VCF
bcftools index -t $NORMED_VCF

#for loop for each control file provided by user
for CONTROL in "${CONTROL_BAM_LIST[@]}"; do
    MADE_ALIGN_CONTROL=false
	CONTROL_BAM=""   
 
    #check to see if the provided bam file is aligned
    if [ "$(samtools view -H "$CONTROL" | grep -c '^@SQ')" -gt 0 ]; then
        CONTROL_BAM=$CONTROL
    else
        $BWA mem -t 40 $REFERENCE_FILE $CONTROL | samtools view -S -@ 12 -b - > $CONTROL_ALIGNED
        CONTROL_BAM=$CONTROL_ALIGNED
	    MADE_ALIGN_CONTROL=true
    fi
    
    #run pileup and call variants
	MERGED_PILEUP="merged_pileup.vcf"
	echo "Starting parallel mpileup..."
	
	# Split pileup by chromosomes
	bcftools query -f '%CHROM\n' $NORMED_VCF | sort | uniq | \
	awk -v bam="$CONTROL_BAM" -v ref="$REFERENCE_FILE" '{print $1 "\t" bam "\t" ref}' > arguments.txt
	cat arguments.txt | parallel -j +0 --colsep '\t' bash $PILEUP_SCRIPT {1} {2} {3}
	
	# TODO: if we have a small amount of sites, might be more efficient to do this, bricks if too many though
	# Split pileup by sites
	#bcftools query -f '%CHROM\t%POS0\t%POS\n' $NORMED_VCF | sort | uniq | \
	#awk -v bam="$CONTROL_BAM" -v ref="$REFERENCE_FILE" '{print $1 "\t" $2 "\t" $3 "\t" bam "\t" ref}' > arguments.txt
	# NOTE: have to do out of order now, because single pileup start and end are optional
	#cat arguments.txt | parallel -j +0 --colsep '\t' bash $PILEUP_SCRIPT {1} {4} {5} {2} {3}

	# combine pileups
	bcftools concat -o $MERGED_PILEUP -Ov mpileup*.vcf
	bcftools sort -o "sorted.$MERGED_PILEUP" "$MERGED_PILEUP"
	bgzip "sorted.$MERGED_PILEUP"
	bcftools index "sorted.$MERGED_PILEUP.gz"
	rm $MERGED_PILEUP
	
	#TODO: comment back in after testing
	#rm mpileup_*.vcf
	#rm arguments.txt		

	# call variants from merged pileup vcf
	echo "Starting pileup call..."
	bcftools call -cv -Oz -o $CONTROL_VCF "sorted.$MERGED_PILEUP.gz"
    bcftools index -t $CONTROL_VCF
	rm "sorted.$MERGED_PILEUP.gz"*    
	
    #intersect the control vcf with formatted rufus vcf
	echo "Starting intersection..."
    bcftools isec -Oz -w1 -n=1 -p $ISEC_OUT_DIR $NORMED_VCF $CONTROL_VCF    
 
    # save the new vcf as rufus final vcf
    OUTFILE="$ISEC_OUT_DIR/0000.vcf.gz"
    OUT_INDEX="$ISEC_OUT_DIR/0000.vcf.gz.tbi"

    cp "$OUTFILE" "${OUT_VCF}" 
    cp "$OUT_INDEX" "${OUT_VCF}.tbi" 
	rm "$CONTROL_VCF"*

    # clean up aligned control file, if it exists
	if [ "$MADE_ALIGN_CONTROL" = "true" ]; then
    	rm $CONTROL_ALIGNED
	fi

	rm -r $ISEC_OUT_DIR
done
