#!/bin/bash


############################# LOADING PARAMETERS #############################

# Working directory: absolute path to the folder containing all files.
WD=$(grep working_directory: $1 | awk '{ print $2 }')
# Genome folder absolute path.
GENOME=$(grep genome_path: $1 | awk '{ print $2 }')
# Number of processors.
THREADS=$(grep number_processors: $1 | awk '{ print $2 }')



## A. Calculating coverage over the designated bins for later analysis ##

############################# GETTING GENOMIC BINS #############################

cd $GENOME

bedtools makewindows -g hg38.genome -w 50000 > genome_50kb.bed

############################# COVERAGE OVER REGIONS #############################

cd ${WD}/results

multiBamSummary BED-file -p 40 --BED /home/gmzlab/Documents/ricardo/genome/homo_sapiens_unknown/genome_50kb.bed -b ../bam_files/*.bam -o results.npz --outRawCounts output_rawCount.txt

cat output_rawCount.txt | sort -k1,1 -k2,2n > output_rawCount.sort.txt




## B. Preparing and training the HMM model for chromatin states

############################# CONVERTING BAMS TO BEDS #############################

for i in `ls bam_files_merge/h3k*.bam bam_files_merge/input*.bam | cut -d '/' -f 2 | cut -d '.' -f 1` 
do
	echo $i
	bedtools bamtobed -i bam_files_merge/${i}.bam > bed_files/${i}.bed
done

############################# BINARIZING GENOME COVERAGE #############################

cd ${WD}

java -jar /home/gmzlab/opt/ChromHMM/ChromHMM.jar BinarizeBed -b 50000 /home/gmzlab/Documents/ricardo/genome/homo_sapiens_ucsc/chromsizes bed_files metadata/cell_mark_table_bed.tsv results/binarize

############################# CHROMATIN STATES #############################

# in results/chromHMM
java -jar ../../opt/ChromHMM/ChromHMM.jar LearnModel -p 12 -b 50000 binarize_agustin_50kb/ model_50kb 4 hg38 

#############################  #############################



## C. Obtain the correlation between the emission parameters of differente models and a reference ##

############################# CORRELATION BETWEEN MODELS #############################

java -jar ../../opt/ChromHMM/ChromHMM.jar CompareModels

#############################  #############################

