#$ -S /bin/bash
#$ -cwd

THR=$1
SAMPLE_NAME=$2

## Accession to the samples folders

red='\033[0;31m'
nc='\033[0m'
gr='\033[0;32m'

echo -e "${red}\n"
date
echo "PROCESSING $FILE..."
echo -e "${red}----------------------------------------------${nc}"

#echo -e "${gr}\nDecompressing...${nc}"
#gunzip ${SAMPLE_NAME}.fastq.gz

echo -e "${gr}\nMapping with bwa...${nc}"
#bwa aln -n 3 -k 2 -R 300 -t ${THR} genome/genome.fa ${SAMPLE_NAME}.fastq > ${SAMPLE_NAME}.sai
bwa samse -n 3 genome/genome.fa ${SAMPLE_NAME}.sai ${SAMPLE_NAME}.fastq > ${SAMPLE_NAME}.sam

echo -e "${gr}\nFiltering unique mapping reads...${nc}"
samtools view -F 256 -O BAM -@ $THR ${SAMPLE_NAME}.sam > ${SAMPLE_NAME}.smap.bam

echo -e "${gr}\n Removing PCR duplicates...${nc}"
samtools rmdup -s ${SAMPLE_NAME}.smap.bam ${SAMPLE_NAME}.unsorted.bam

## Deleting fastq because it weights a lot and there's not enough space

echo -e "${gr}\n Sorting .bam file...${nc}"
samtools sort ${SAMPLE_NAME}.unsorted.bam > ${SAMPLE_NAME}.bam

echo -e "${gr}\n Indexing .bam file...${nc}"
samtools index ${SAMPLE_NAME}.bam

echo -e "${gr}\n Cleaning up...${nc}"
#rm ${SAMPLE_NAME}.fastq ${SAMPLE_NAME}.sai ${SAMPLE_NAME}.sam ${SAMPLE_NAME}.smap.bam ${SAMPLE_NAME}.unsorted.bam

