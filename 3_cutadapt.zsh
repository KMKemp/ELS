#!/bin/zsh

mkdir ./cutadapt/trimmed_reads
for sample in $(cat samples)
do
    echo "On sample: $sample"
    
    ~/.local/bin/cutadapt -g GTGYCAGCMGCCGCGGTAA -G GGACTACNVGGGTWTCTAAT \
    -m 220 -M 285 \
    -o ./cutadapt/trimmed_reads/${sample}_R1_trimmed.fq.gz \
    -p ./cutadapt/trimmed_reads/${sample}_R2_trimmed.fq.gz \
    ./trimmomatic/trimmed/${sample}_1P.fastq.gz \
    ./trimmomatic/trimmed/${sample}_2P.fastq.gz \
    >> cutadapt_primer_trimming_stats.txt 2>&1
done