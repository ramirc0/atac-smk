#!/bin/bash

mkdir -p test

bowtie2 -X 2000 --mm -k 4 \
    --threads 96 \
    -x results/bowtie2_build/GRCH38/GRCH38.fa \
    --rg-id EA100050-GRCH38 \
    --rg SM:EA100050-GRCH38 \
    -1 results/fastp/EA100001_R1.fastq.gz \
    -2 results/fastp/EA100001_R2.fastq.gz | \
samtools view -@ 96 -1 /dev/stdin | \
samtools sort -@ 96 -T EA100050-GRCH38 -o test/EA100050-GRCH38.bam
samtools index -@ 96 test/EA100050-GRCH38.bam
samtools flagstat -@ 96 --output-fmt json test/EA100050-GRCH38.bam

