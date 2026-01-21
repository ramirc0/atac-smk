#!/bin/bash

sambamba markdup -t 24 results/bowtie2_align/EA100050-GRCH38.bam /dev/stdout | \
samtools view -@ 24 -F 1804 -f 2 -q 30 -h -u /dev/stdin | \
samtools sort -@ 24 -n /dev/stdin | \
samtools fixmate -@ 24 -r /dev/stdin /dev/stdout | \
samtools view -@ 24 -F 1804 -f 2 -h -u /dev/stdin | \
samtools sort -@ 24 /dev/stdin | \
samtools view -@ 24 -h | grep -v "chrM" | samtools view -@ 24 -b > test/EA100050-GRCH38.filt.dedup.bam