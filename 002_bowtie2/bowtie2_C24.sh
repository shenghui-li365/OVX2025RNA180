#!/bin/bash
#DSUB -n bowtie2_C24
#DSUB -o bowtie2C24_%J.out
#DSUB -e bowtie2C24_%J.err
#DSUB -R cpu=4
module load arm/bowtie2/2.5.4
module load arm/samtools/1.21
bowtie2 -p 4 -x /share/home/shli24/db/GRCm39_bt2/index115 -1 /share/home/shli24/OVX2025CUT/data/clean_data/C24_1_clean.fq.gz -2 /share/home/shli24/OVX2025CUT/data/clean_data/C24_2_clean.fq.gz |samtools sort -@ 28 -o /share/home/shli24/2025OVX/data/bam/C24.bam
