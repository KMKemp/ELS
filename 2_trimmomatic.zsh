#!/bin/zsh

mkdir ./trimmomatic/trimmed

for sample in $(cat samples) 
do
    
trimmomatic PE -threads 4 ./raw_reads/${sample}_R1.fastq.gz ./raw_reads/${sample}_R2.fastq.gz -baseout  ./trimmomatic/trimmed/${sample}.fastq.gz ILLUMINACLIP:NexteraPE-PE.fa:2:30:10 MINLEN:30 \
 >> ./trimmomatic/trimmomatic_trimming_log.txt 2>&1

 
done
