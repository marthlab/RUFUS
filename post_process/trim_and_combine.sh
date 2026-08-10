#!/bin/bash

# This script concatenates individual vcfs produced by piecemeal RUFUS runs,
# into a single combined vcf. Only variants from within the region targeted by the run
# are kept (at this time, RUFUS may make erroneous calls out of the given range). 
# The final vcf contains the header from the first chr1_1_{CHUNK_SIZE} final 
# individual vcf. Note, the chr and chr length arrays, as well as the chunk size
# must be equal to those run for the piecemeal run.
# Compresses and indexes the final file.

# ENV override
: "${RUFUS_ROOT:=/opt/RUFUS}"

SUBJECT_FILE=$1
WINDOW_SIZE=$2
shift
shift
CONTROLS=("$@")

BCFTOOLS="/opt/bcftools/bcftools"

COMBINED_VCF="temp.RUFUS.Final.${SUBJECT_FILE}.combined.vcf"
#COMBINED_PRE_VCF="temp.RUFUS.Prefiltered.${SUBJECT_FILE}.combined.vcf"

COMBINED_SAMPLE_STRING=""
if [ "$CONTROL_STRING" == "internal" ]; then
    COMBINED_SAMPLE_STRING="${SUBJECT_FILE}"
else
    COMBINED_SAMPLE_STRING=$(printf '%s' "$SUBJECT_FILE"; printf '\t%s' "${CONTROLS[@]}")
fi

#SUPP_DIR="rufus_supplementals/"

# Headers that get written to vcf
COMBINED_HEADER="combined.header"
COMBINED_PRE_HEADER="combined.preheader"

# Start of headers
HEADER_STUB="$RUFUS_ROOT/resources/vcf_header.txt"

# Records that get written to vcf (non-header)
COMBINED_RECORDS="combined.records"
COMBINED_PRE_RECORDS="combined.prerecords"

contig_temp="contig_temp.txt"

CHRS=(
"1"  
"2" 
"3" 
"4" 
"5" 
"6" 
"7" 
"8" 
"9" 
"10" 
"11" 
"12" 
"13" 
"14" 
"15" 
"16" 
"17" 
"18" 
"19" 
"20" 
"21" 
"22" 
"X" 
"Y" 
)

# GRCh38
CHR_LENGTHS=( 
248956422  
242193529 
198295559 
190214555 
181538259 
170805979 
159345973 
145138636 
138394717 
133797422 
135086622 
133275309 
114364328 
107043718 
101991189 
90338345 
83257441 
80373285 
58617616 
64444167 
46709983 
50818468 
156040895 
57227415 
)

# Initialize combined headers
cat "$HEADER_START" > $COMBINED_HEADER
#cat "$PRE_HEADER_START" > $COMBINED_PRE_HEADER
TEMP_TRIMMED="temp.trimmed"

# Adjust chunk size to bp
CHUNK_SIZE=$((WINDOW_SIZE * 1000))

n=${#CHRS[@]}
for (( i = 0; i < n; i++ ))
do
    curr_len="${CHR_LENGTHS[i]}"
    curr_chr="${CHRS[i]}"
    ((chunks=curr_len/CHUNK_SIZE+1))
 
    start_coord=1
    for (( c = 1; c <= chunks; c++ ))
    do
        ((end_coord=start_coord+CHUNK_SIZE-1))
        if [ $end_coord -gt $curr_len ]; then
            echo "passed end coordinate, setting to $curr_len"
            end_coord=$curr_len
        fi
	
        CURR_VCF="temp.RUFUS.Final.${SUBJECT_FILE}.chr${curr_chr}_${start_coord}_${end_coord}.vcf.gz"
        #CURR_PRE_VCF="temp.RUFUS.Prefiltered.${SUBJECT_FILE}.chr${curr_chr}_${start_coord}_${end_coord}.vcf.gz"
        if [[ -f "${CURR_VCF}" ]]; then

            # Write out trimmed region to final vcf
            bcftools view -r "chr${curr_chr}:${start_coord}-${end_coord}" "${CURR_VCF}" > $TEMP_TRIMMED
            bcftools view -H $TEMP_TRIMMED >> $COMBINED_RECORDS
            bcftools view -h $TEMP_TRIMMED | grep "##contig" >> $contig_temp

			# Write out trimmed region to prefiltered vcf 
            #bcftools view -r "chr${curr_chr}:${start_coord}-${end_coord}" "${CURR_PRE_VCF}" > $TEMP_TRIMMED
            #bcftools view -H $TEMP_TRIMMED >> $COMBINED_PRE_RECORDS
            #bcftools view -h $TEMP_TRIMMED | grep "##contig" >> $contig_temp
            #sort -V $contig_temp | uniq >> $COMBINED_PRE_HEADER

	    	# Remove vcf and indexes
	    	rm "$CURR_VCF"*
	    	#rm $CURR_PRE_VCF*
        fi
    
        # Advance start coordinate
        ((start_coord=end_coord+1))
    done
    echo "Done combining chr${curr_chr}"
done

sort -hu $contig_temp >> "$COMBINED_HEADER"
cat $COMBINED_HEADER | uniq > "$COMBINED_VCF"
echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$COMBINED_SAMPLE_STRING" >> "$COMBINED_VCF"
cat $COMBINED_RECORDS >> "$COMBINED_VCF"

# cat $COMBINED_PRE_HEADER | uniq > $COMBINED_PRE_VCF
# echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$COMBINED_SAMPLE_STRING" >> $COMBINED_PRE_VCF
# cat $COMBINED_PRE_RECORDS >> $COMBINED_PRE_VCF

# todo: left off here - unsure why the combined vcf here is out of order and cannot be indexed...
# maybe try to do a full redo of individual vcfs and try again?
bgzip $COMBINED_VCF
bcftools index -t "${COMBINED_VCF}.gz"

# bgzip $COMBINED_PRE_VCF
# bcftools index -t "${COMBINED_PRE_VCF}.gz"

# Clean up temp files
rm $contig_temp
rm $TEMP_TRIMMED
rm $COMBINED_HEADER
# rm $COMBINED_PRE_HEADER
rm $COMBINED_RECORDS
#rm $COMBINED_PRE_RECORDS