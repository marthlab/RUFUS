#!/bin/bash
# NOTE: This script requires bash (process substitution is used)

# Despite it's name, this is the only script currently utilized by RUFUS, despite any mode. When speed is true, veryfast
# mode is used.
#
# veryfast mode consists of the following criteria:
# 1. Reads going into OverlapSam are filtered by length - must be 150bp or shorter
# 2. OverlapSam is run with the following parameters:
#   a. Min Percentage is 99%
#   b. Min Overlap is 25bp
#   c. Min Coverage is 3bp - todo: this needs to be set to command line arg
#
# normal mode consists of the following criteria:
# 1. OverlapSam is run with the following parameters:
#   a. Min Percentage is 95%
#   b. Min Overlap is 20bp
#   c. Min Coverage is 1bp


set -euo pipefail

: "${WORK_DIR:?WORK_DIR must be set}"


humanRef=$1
File=$2 # e.g. WGS_IL_T_1.bwa.dedup.bam.generator.Mutations.fastq
FinalCoverage=$3 #todo: what is the difference between this and minOverlap
NameStub=$4.V2 # e.g. WGS_IL_T_1.bwa.dedup.bam.generator.Mutations.fastq
HashList=$5 # e.g. $ProbandGenerator".k"$K"_c"$MutantMinCov".HashList
HashSize=$6
Threads=$7
MaxAlleleSize=$8
speed=$9
humanRefBwa=${10}
invocFilePath=${11}
refHash=${12} # Will say "empty" if not provided
SampleJhash=${13}
ParentsJhash=${14} # this is optional
ParLowCovThreshold=${15:-7}


MaxCov=100000
#echo " you gave
#File=$2
#FinalCoverage=$3
#NameStub=$4.V2
#HashList=$5
#HashSize=$6
#Threads=$7
#"

echo "final coverage is $FinalCoverage"
#echo "final coveage is $FinalCoverage"


#echo "RUNNING THIS ONE"
#echo "@@@@@@@@@@@@@__IN_OVERLAP__@@@@@@@@@@@@@@@"
#echo "human ref in Overlap is $humanRef"
#echo "bwa human ref in Overlap is $humanRefBwa"
#echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
echo "Overlaping $File..."

: "${RUFUS_ROOT:=/opt/RUFUS}"
RDIR="$RUFUS_ROOT"

AddSA=$RDIR/scripts/AddSAtoReadSame.pl
OverlapHash=$RDIR/bin/Overlap
OverlapRebion2=$RDIR/bin/OverlapRegion
OverlapRegionSmall=$RDIR/bin/OverlapRegion.small
ReplaceQwithDinFASTQD=$RDIR/bin/ReplaceQwithDinFASTQD
ConvertFASTqD=$RDIR/bin/ConvertFASTqD.to.FASTQ
AnnotateOverlap=$RDIR/bin/AnnotateOverlap
bwa=$RDIR/bin/externals/bwa/src/bwa_project/bwa
RUFUSinterpret=$RDIR/bin/RUFUS.interpret
CheckHash=$RDIR/scripts/CheckJellyHashList.sh
OverlapSam=$RDIR/bin/OverlapSam
JellyFish=$RDIR/bin/externals/jellyfish/src/jellyfish_project/bin/jellyfish
MOBList=$RDIR/resources/primate_non-LTR_Retrotransposon.fasta

if [ -s "$WORK_DIR/$File.bam" ] 
then 
	echo "skipping align"
else
    sortedFastq="$WORK_DIR/sorted."$File
    cat "$WORK_DIR/$File" | paste - - - - | sort -k1 -S 8G | tr "\t" "\n" > $sortedFastq
    
	"$bwa" mem -t $Threads "$humanRefBwa" "$sortedFastq" | samtools view -h - | samtools sort -T $File -O bam - > "$File.bam"
	samtools index "$File.bam" 
fi

if [ $( samtools view "$File.bam" | head | wc -l | awk '{print $1}') -eq "0" ]; then
        echo "ERROR: BWA failed on $File .  Either the files are exactly the same of something went wrong in previous step" 
        exit 100
fi

if [ "$speed" = "veryfast" ]
then
	echo "running very fast assembly"; 
	if [ -s "$WORK_DIR/TempOverlap/$NameStub.sam.fastqd" ]
	then
	        echo "skipping sam assemble"
	else
		"$OverlapSam" <( samtools view  -F 3328 "$File.bam" | awk '$9 > 150 || $9 < -150 '  ) .99 25 $FinalCoverage "$WORK_DIR/TempOverlap/$NameStub.sam" $NameStub 1 "$HashList" $Threads
	fi 
	# todo: instead of hash here, first do overlapRegion looking at next 5 reads
	if [ -s "$WORK_DIR/TempOverlap/$NameStub.final.fastqd" ]
	then 
		echo "skipping second assemble"
	else 
		# This will fix the gaps and also do trans-chromosomal alignments
		"$OverlapHash" "$WORK_DIR/TempOverlap/$NameStub.sam.fastqd" .99 75 $FinalCoverage $NameStub 15 1 "$WORK_DIR/TempOverlap/$NameStub.final" 1 $Threads 
	fi
	
	if [ -s "$WORK_DIR/$NameStub.overlap.hashcount.fastq" ]
	then
	        echo "skipping final overlap work"
	else
	        $ReplaceQwithDinFASTQD "$WORK_DIR/TempOverlap/$NameStub.final.fastqd" > "$WORK_DIR/$NameStub.overlap.fastqd"
	        $ConvertFASTqD "$WORK_DIR/$NameStub.overlap.fastqd" > "$WORK_DIR/$NameStub.overlap.fastq"
	
	        #echo "$AnnotateOverlap "$HashList" $WORK_DIR/$NameStub.overlap.fastq $WORK_DIR/TempOverlap/$NameStub.overlap.asembly.hash.fastq > $WORK_DIR/$NameStub.overlap.hashcount.fastq"              
	        $AnnotateOverlap "$HashList" "$WORK_DIR/$NameStub.overlap.fastq" "$WORK_DIR/TempOverlap/$NameStub.overlap.asembly.hash.fastq" > "$WORK_DIR/$NameStub.overlap.hashcount.fastq"
	fi	
else
	echo "Running full assembly"; 
	
	if [ -s "$WORK_DIR/TempOverlap/$NameStub.sam.fastqd" ]  
	then 
		echo "skipping sam assemble"
	else
		"$OverlapSam" <( samtools view  -F 3328 "$File.bam" ) .95 20 1 "$WORK_DIR/TempOverlap/$NameStub.sam" $NameStub 1 "$HashList" $Threads
	fi

	# todo: the problem here is that this is empty
	#if [ $( wc -l $WORK_DIR/TempOverlap/$NameStub.sam.fastqd | awk '{print $1}') -eq "0" ]; then
  if [ $( head "$WORK_DIR/TempOverlap/$NameStub.sam.fastqd" | wc -l | awk '{print $1}') -eq "0" ]; then
		echo "ERROR Assembly produce output for $WORK_DIR/TempOverlap/$NameStub.sam.fastqd"
		exit 100
	fi
	
	if [ -s "$WORK_DIR/TempOverlap/$NameStub.1.fastqd" ]
	then 	
		echo "skipping first overlap"
	else
		echo "$OverlapHash $WORK_DIR/TempOverlap/$NameStub.sam.fastqd .98 100 1 FP 20 1 $WORK_DIR/TempOverlap/$NameStub.1 0 $Threads"
		$OverlapHash "$WORK_DIR/TempOverlap/$NameStub.sam.fastqd" .98 100 1 FP 20 1 "$WORK_DIR/TempOverlap/$NameStub.1" 0 $Threads #> $File.overlap.out
	fi
	
	if [ $( head "$WORK_DIR/TempOverlap/$NameStub.1.fastqd" | wc -l | awk '{print $1}') -eq "0" ]; then
	        echo "ERROR Assembly produce output for $WORK_DIR/TempOverlap/$NameStub.1.fastqd"
	        exit 100
	fi
	
	if [ -s "$WORK_DIR/TempOverlap/$NameStub.2.fastqd" ]
	then 
		echo "skipping second overlap"
	else
		$OverlapHash "$WORK_DIR/TempOverlap/$NameStub.1.fastqd" .98 75 2 FP 20 1 "$WORK_DIR/TempOverlap/$NameStub.2" 1 $Threads #>>  $File.overlap.out
	fi
	
	if [ $( head "$WORK_DIR/TempOverlap/$NameStub.2.fastqd" | wc -l  | awk '{print $1}') -eq "0" ]; then
	        echo "ERROR Assembly produce output for $WORK_DIR/TempOverlap/$NameStub.2.fastqd"
	        exit 100
	fi
	
	if [ -s "$WORK_DIR/TempOverlap/$NameStub.3.fastqd" ]
	then
	        echo "skipping third overlap"
	else
	        $OverlapHash $WORK_DIR/TempOverlap/$NameStub.2.fastqd .98 50 2 $NameStub 20 1 $WORK_DIR/TempOverlap/$NameStub.3 1 $Threads #>>  $File.overlap.out
	fi
	if [ $( head "$WORK_DIR/TempOverlap/$NameStub.3.fastqd" | wc -l | awk '{print $1}') -eq "0" ]; then
	        echo "ERROR Assembly produce output for $WORK_DIR/TempOverlap/$NameStub.3.fastqd"
	        exit 100
	fi
	
	if [ -s "$WORK_DIR/TempOverlap/$NameStub.4.fastqd" ]
	then
	        echo "skipping fourth overlap"
	else
			# This is the original version - every read 
	        $OverlapRebion2 $WORK_DIR/TempOverlap/$NameStub.3.fastqd .98 50 $FinalCoverage  $WORK_DIR/TempOverlap/$NameStub.4 $NameStub 1 $Threads 
	fi
	if [ $( head "$WORK_DIR/TempOverlap/$NameStub.4.fastqd" | wc -l | awk '{print $1}') -eq "0" ]; then 
	        echo "ERROR Assembly produce output for $WORK_DIR/TempOverlap/$NameStub.4.fastqd"
	        exit 100
	fi
 

	if [ -s "$WORK_DIR/$NameStub.overlap.hashcount.fastq" ]
	then 
		echo "skipping final overlap work"
	else
	
		$ReplaceQwithDinFASTQD "$WORK_DIR/TempOverlap/$NameStub.4.fastqd" > "$WORK_DIR/$NameStub.overlap.fastqd"
		$ConvertFASTqD $WORK_DIR/$NameStub.overlap.fastqd > $WORK_DIR/$NameStub.overlap.fastq
	
		echo "$AnnotateOverlap "$HashList" $WORK_DIR/$NameStub.overlap.fastq $WORK_DIR/TempOverlap/$NameStub.overlap.asembly.hash.fastq > $WORK_DIR/$NameStub.overlap.hashcount.fastq"              
		      $AnnotateOverlap "$HashList" $WORK_DIR/$NameStub.overlap.fastq $WORK_DIR/TempOverlap/$NameStub.overlap.asembly.hash.fastq > $WORK_DIR/$NameStub.overlap.hashcount.fastq
	fi
fi

if [ $( head "$WORK_DIR/$NameStub.overlap.hashcount.fastq" | wc -l | awk '{print $1}') -eq "0" ]; then 
        echo "RUFUS could not assemble any contigs from unique reads for the given region. Exiting..."
        exit 0
fi

# Sort fastq file used in subsequence bwa calls for reproducibility
sortedFastq=$WORK_DIR/$NameStub".overlap.hashcount.sorted.fastq" 
if [ -s "$WORK_DIR/$NameStub.overlap.hashcount.fastq" ]
then
    echo "Sorting hashcount fastq file"
    cat "$WORK_DIR/$NameStub.overlap.hashcount.fastq" | paste - - - - | sort -k1 -S 8G | tr "\t" "\n" > $sortedFastq
else
    echo "$NameStub.overlap.hashcount.fastq does not exist, cannot sort"
    echo "Exiting with failure"
    exit 100
fi

if [ -s "$WORK_DIR/$NameStub.overlap.hashcount.fastq.bam" ]
then 
	echo "skipping contig alignment" 
else
    "$bwa" mem -t $Threads -Y  "$humanRefBwa" "$sortedFastq" | samtools view -h - | samtools sort -T $File -O bam - > $WORK_DIR/$NameStub.overlap.hashcount.fastq.bam
	samtools index "$WORK_DIR/$NameStub.overlap.hashcount.fastq.bam"
fi

if [ $( samtools view "$WORK_DIR/$NameStub.overlap.hashcount.fastq.bam" | head | wc -l | awk '{print $1}') -eq "0" ]; then
        echo "ERROR: BWA failed on $WORK_DIR/$NameStub.overlap.hashcount.fastq.bam .  Either the files are exactly the same of something went wrong in previous step" 
        exit 100
fi

echo "string hash lookup"
#############################################################################################################
echo "staring MOB check on sorted fastq"
if [ -s "$WORK_DIR/Intermediates/$NameStub.overlap.hashcount.fastq.MOB.sam" ]
then
	echo "skipping MOB alignemnt check "
else
	echo "$bwa mem -t $Threads -Y -E 0,0 -O 6,6  -d 500 -w 500 -L 0,0 $MOBList "$sortedFastq" | samtools view -h - | samtools sort -T $File -O sam - > "$WORK_DIR/Intermediates/$NameStub.overlap.hashcount.fastq.MOB.sam""
	"$bwa" mem -t $Threads -Y -E 0,0 -O 6,6  -d 500 -w 500 -L 0,0 $MOBList "$sortedFastq" | samtools view -h - | samtools sort -T $File -O sam - > "$WORK_DIR/Intermediates/$NameStub.overlap.hashcount.fastq.MOB.sam"
fi 

if [ -e "$WORK_DIR/Intermediates/$NameStub.overlap.asembly.hash.fastq.ref.fastq" ]
then 
	echo "skipping pull reference sequecnes"
else
	bedtools getfasta -bed <( bedtools bamtobed -i "$WORK_DIR/$NameStub.overlap.hashcount.fastq.bam" |  awk '{s=$2-100; if (s<0) {print $1 "\t" 0  "\t" $3+100} else {print $1 "\t" s  "\t" $3+100}}'  ) -fi "$humanRef" -fo "$WORK_DIR/Intermediates/$NameStub.overlap.asembly.hash.fastq.ref.fastq"

fi 

if [ -e "$WORK_DIR/Intermediates/$NameStub.overlap.hashcount.fastq.Jhash" ]
then 
	echo "skipping var hash generation"
else
	#echo "$JellyFish count -m $HashSize -s 1G -t 20 -o $WORK_DIR/Intermediates/$NameStub.overlap.hashcount.fastq.Jhash $WORK_DIR/$NameStub.overlap.hashcount.fastq"
	$JellyFish count -m $HashSize -s 1G -t 1 -o $WORK_DIR/Intermediates/$NameStub.overlap.hashcount.fastq.Jhash $WORK_DIR/$NameStub.overlap.hashcount.fastq
	#echo "$JellyFish dump  -c $WORK_DIR/Intermediates/$NameStub.overlap.hashcount.fastq.Jhash > $WORK_DIR/Intermediates/$NameStub.overlap.hashcount.fastq.Jhash.tab"
	$JellyFish dump  -c $WORK_DIR/Intermediates/$NameStub.overlap.hashcount.fastq.Jhash > $WORK_DIR/Intermediates/$NameStub.overlap.hashcount.fastq.Jhash.tab
fi 

echo "Creating reference kMer hash..."
if [ -s "$WORK_DIR/Intermediates/$NameStub.overlap.asembly.hash.fastq.ref.fastq.Jhash" ] 
then 
	echo "skipping ref hash generation"
else
	$JellyFish count -m $HashSize -s 1G -t 1 -o $WORK_DIR/Intermediates/$NameStub.overlap.asembly.hash.fastq.ref.fastq.Jhash $WORK_DIR/Intermediates/$NameStub.overlap.asembly.hash.fastq.ref.fastq
	$JellyFish dump -c  $WORK_DIR/Intermediates/$NameStub.overlap.asembly.hash.fastq.ref.fastq.Jhash > $WORK_DIR/Intermediates/$NameStub.overlap.asembly.hash.fastq.ref.fastq.Jhash.tab
fi
 
 echo "Retrieving kMer hashes from sample..." 
if [ -s "$WORK_DIR/Intermediates/$NameStub.overlap.asembly.hash.fastq.sample" ]
then
        echo "skipping  $WORK_DIR/Intermediates/$NameStub.overlap.asembly.hash.fastq.sample file already exitst"
else
	echo "starting hash lookup this one"
        bash $CheckHash $SampleJhash $WORK_DIR/Intermediates/$NameStub.overlap.hashcount.fastq.Jhash.tab 0 $MaxCov> $WORK_DIR/Intermediates/$NameStub.overlap.asembly.hash.fastq.sample & 
		pid=$!
	echo "done with hash lookup"
	wait "$pid" || { echo "ERROR: hash lookup failed"; exit 100; }
fi

echo "Retrieving kMer hashes from control(s)..." 
IFS=' ' read -r -a parents <<< "$ParentsJhash"

# Per-control hash intermediates below concatenate the SUBJECT stub ($NameStub) and the CONTROL name
# ($parent). With long input filenames that basename can exceed the 255-byte NAME_MAX and fail to
# create ("File name too long" -> exit 100 in the overlap stage; hit by ~96-char SMaHT CRAM names).
# RUFUS.interpret derives each control's VCF sample name from the text AFTER the marker
# "overlap.asembly.hash.fastq." (i.e. the $parent portion) and ignores the prefix, so we use a SHORT
# constant prefix here instead of $NameStub. The control name is preserved intact -> output VCF is
# byte-identical; only the on-disk intermediate filename gets shorter. (One subject per WORK_DIR, so a
# constant prefix cannot collide.)
CtrlStub="ctrlhash"
for parent in "${parents[@]}"
        do
			parent_path="$parent"; parent=$(basename "$parent")   # keep full path for the hash query; basename is for on-disk filenames only
            if [ -s "$WORK_DIR/Intermediates/$CtrlStub.overlap.asembly.hash.fastq.$parent" ]
            then
                echo "skiping $WORK_DIR/Intermediates/$CtrlStub.overlap.asembly.hash.fastq.$parent already exists"
            else
				echo "pulling  $WORK_DIR/Intermediates/$CtrlStub".overlap.asembly.hash.fastq."$parent"
                bash $CheckHash "$parent_path" $WORK_DIR/Intermediates/$NameStub.overlap.hashcount.fastq.Jhash.tab 0 $MaxCov> $WORK_DIR/Intermediates/$CtrlStub".overlap.asembly.hash.fastq."$parent &
						pid=$!
				wait "$pid" || { echo "ERROR: hash lookup failed"; exit 100; }
            fi
done

wait

if [ -s "$WORK_DIR/Intermediates/$NameStub.overlap.asembly.hash.fastq.Ref.sample" ]	
then 
	echo "skipping $WORK_DIR/Intermediates/$NameStub.overlap.asembly.hash.fastq.Ref.sample"
else
	bash $CheckHash $SampleJhash $WORK_DIR/Intermediates/$NameStub.overlap.asembly.hash.fastq.ref.fastq.Jhash.tab 0 $MaxCov> $WORK_DIR/Intermediates/$NameStub.overlap.asembly.hash.fastq.Ref.sample&
fi

for parent in "${parents[@]}"
do
	parent_path="$parent"; parent=$(basename "$parent")   # keep full path for the hash query; basename is for on-disk filenames only
    if [ -s "$WORK_DIR/Intermediates/$CtrlStub.overlap.asembly.hash.fastq.Ref.$parent" ]
    then
        echo "skipping $CtrlStub.overlap.asembly.hash.fastq.Ref.$parent already exitst"
    else
	
        #echo "-$parent-"
        #echo "  bash $CheckHash $parent $WORK_DIR/Intermediates/$NameStub.overlap.asembly.hash.fastq.ref.fastq.Jhash.tab 0 $MaxCov> $WORK_DIR/Intermediates/$CtrlStub.overlap.asembly.hash.fastq.Ref.$parent"
        bash $CheckHash "$parent_path" $WORK_DIR/Intermediates/$NameStub.overlap.asembly.hash.fastq.ref.fastq.Jhash.tab 0 $MaxCov> $WORK_DIR/Intermediates/$CtrlStub.overlap.asembly.hash.fastq.Ref.$parent &
    	echo "uncomment this"
    fi
done
wait
parentCRString=""
c="-c"
cr="-cR"
space=" "


######################## BUILDING UP parent c and cR string ##############################
for parent in "${parents[@]}";
do
	parent_path="$parent"; parent=$(basename "$parent")   # keep full path for the hash query; basename is for on-disk filenames only
    parentCRString="$parentCRString -c $WORK_DIR/Intermediates/$CtrlStub.overlap.asembly.hash.fastq.$parent -cR $WORK_DIR/Intermediates/$CtrlStub.overlap.asembly.hash.fastq.Ref.$parent "
done

##########################################################################################
if [ -s "$WORK_DIR/Intermediates/$NameStub.ref.RepRefHash" ]
then
        echo "Exclude already exists"
else
	if [ "$refHash" = "empty" ]
	then 
		echo "refhash not provided, skipping"
		touch  $WORK_DIR/Intermediates/$NameStub.ref.RepRefHash
	else
		
		#echo "this one" 
		#echo "bash $CheckHash $refHash $WORK_DIR/Intermediates/$NameStub.overlap.hashcount.fastq.Jhash.tab 0 $MaxCov > $WORK_DIR/Intermediates/$NameStub.ref.RepRefHash"
		bash $CheckHash $refHash $WORK_DIR/Intermediates/$NameStub.overlap.hashcount.fastq.Jhash.tab 0 $MaxCov> $WORK_DIR/Intermediates/$NameStub.ref.RepRefHash
		#echo "outa this"
	fi
fi
wait

samtools index $WORK_DIR/$NameStub.overlap.hashcount.fastq.bam
echo ""
echo "" 
echo ""
dumbFix=$(awk '{split($1, a, ".V2"); print a[1]}' <<< $NameStub)
#echo "$RUFUSinterpret -mob $WORK_DIR/Intermediates/$NameStub.overlap.hashcount.fastq.MOB.sam -mod $dumbFix.Jhash.histo.7.7.dist -mQ 20 -r $humanRef -hf "$HashList" -o  $WORK_DIR/$NameStub.overlap.hashcount.fastq.bam -m $MaxAlleleSize $(echo $parentCRString) -sR $WORK_DIR/Intermediates/$NameStub.overlap.asembly.hash.fastq.Ref.sample -s $WORK_DIR/Intermediates/$NameStub.overlap.asembly.hash.fastq.sample -e $WORK_DIR/Intermediates/$NameStub.ref.RepRefHash"

samtools view -h $WORK_DIR/$NameStub.overlap.hashcount.fastq.bam | perl $AddSA | grep -v chrUn  | $RUFUSinterpret -mob $WORK_DIR/Intermediates/$NameStub.overlap.hashcount.fastq.MOB.sam -mod $dumbFix.Jhash.histo.7.7.dist -mQ 10 -r "$humanRef" -hf "$HashList" -o $NameStub.overlap.hashcount.fastq.bam -m $MaxAlleleSize $parentCRString -sR $WORK_DIR/Intermediates/$NameStub.overlap.asembly.hash.fastq.Ref.sample -s $WORK_DIR/Intermediates/$NameStub.overlap.asembly.hash.fastq.sample -e $WORK_DIR/Intermediates/$NameStub.ref.RepRefHash -rp "$RUFUS_ROOT" -ip "$invocFilePath" -plct "$ParLowCovThreshold"