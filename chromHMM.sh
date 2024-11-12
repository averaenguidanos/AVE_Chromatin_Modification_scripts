#!/bin/bash



############################# LOADING PARAMETERS #############################

# Working directory: absolute path to the folder containing all files.
WD=$(grep working_directory: $1 | awk '{ print $2 }')
# Genome folder absolute path.
GENOME=$(grep genome_path: $1 | awk '{ print $2 }')
# Number of processors.
THREADS=$(grep number_processors: $1 | awk '{ print $2 }')

############################# GETTING GENOMIC BINS #############################

cd $GENOME

mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e 'select chrom, size from hg38.chromInfo' > hg38.pre_genome

grep -v 'size' hg38.pre_genome > hg38.genome

rm hg38.pre_genome

bedtools makewindows -g hg38.genome -w 50000 -s 50000 > genome_50kb.bed

############################# COVERAGE OVER REGIONS #############################

cd ${WD}/results

multiBamSummary BED-file -p 40 --BED /home/gmzlab/Documents/ricardo/genome/homo_sapiens_unknown/genome_50kb.bed -b ../bam_files/*.bam -o results.npz --outRawCounts output_rawCount.txt

cat output_rawCount.txt | sort -k1,1 -k2,2n > output_rawCount.sort.txt

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

java -jar /home/gmzlab/opt/ChromHMM/ChromHMM.jar LearnModel -p 40 -b 50000 results/binarize results/model 4 hg38

#############################  #############################















