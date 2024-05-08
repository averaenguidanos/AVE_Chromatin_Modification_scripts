#$ -S /bin/bash
#$ -N cycle
#$ -o cycle
#$ -cwd
#$ -j y
#$ -q diaria_simple

# Read params

INS_FOLDER=$1 
WORK_DIR=$2 
MAIN_FOLDER=$3
CHIP_GENERAL=$4
INPUT_GENERAL=$5
THR=$6
HISTONE_GENERAL=$7

## Peak Calling
# Use epic2 for broad peak calling (I have to ask to get it installed in CICA)
# (ultraperformant reimplementation of SICER)
# https://github.com/biocore-ntnu/epic2

# Access to the result directory

cd ${WORK_DIR}/${MAIN_FOLDER}/results
echo "Access to results folder successful" >> ../log/blackboard

# Call peaks on either Input or H3, separately

epic2 --treatment ../samples/chip/${CHIP_GENERAL}.bam --control ../samples/input/${INPUT_GENERAL}.bam --chromsizes ../genome/chromsizes > ${CHIP_GENERAL}_over_${INPUT_GENERAL}.peaks.bed

if [ ! -f ${CHIP_GENERAL}_over_${INPUT_GENERAL}.peaks.bed]
then
	echo "$0: File '${CHIP_GENERAL}_over_${INPUT_GENERAL}.peaks.bed' not found." >> ../log/blackboard_3
fi

#epic2 --treatment ../samples/chip/${CHIP_GENERAL}_1.bam ../samples/chip/${CHIP_GENERAL}_2.bam --control ../samples/histone/${HISTONE_GENERAL}_1.bam ../samples/histone/${HISTONE_GENERAL}_2.bam --chromsizes ../genome/chromsizes > ${CHIP_GENERAL}_over_${HISTONE_GENERAL}.peaks.bed

#if [ ! -f ${CHIP_GENERAL}_over_${HISTONE_GENERAL}.peaks.bed]
#then
#    	echo "$0: File '${CHIP_GENERAL}_over_${HISTONE_GENERAL}.peaks.bed' not found." >> ../log/blackboard_3
#fi

# Extract high confidence sets by intersection of
# both peak calling modes (vs Input, vs H3)

../../../../software/bedtools2-2.30.0/bedtools intersect -wa -u -a ${CHIP_GENERAL}_over_${HISTONE_GENERAL}.peaks.bed \
                          -b ${CHIP_GENERAL}_over_${INPUT_GENERAL}.peaks.bed | grep -v "_" > ${CHIP_GENERAL}.peaks.bed

# Collapse contiguous peaks into single regions

../../../../software/bedtools2-2.30.0/bedtools merge -d 5000 -i ${CHIP_GENERAL}.peaks.bed | sort -n -k7 -r - | cut -f1,2,3 > ${CHIP_GENERAL}.regions.bed

# Quick and dirty strategy to create a linkable table of regions
#bedtools slop -b 0.25 -pct -i ${CHIP_GENERAL}.regions.bed -g ${WORK_DIR}/${MAIN_FOLDER}/genome/genome.chrom.sizes | bedtools links -base "https://genome-euro.ucsc.edu" -db "hg38&hgS_doOtherUser=submit&hgS_otherUserName=LucaPandol$


