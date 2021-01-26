#!/bin/zsh

cd ./raw_reads

ls *_R1.fastq.gz | cut -f1 -d "_" > ../samples
