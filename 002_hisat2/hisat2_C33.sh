#!/bin/bash
#DSUB -n bowtie2_C33
#DSUB -o bowtie2C33_%J.out
#DSUB -e bowtie2C33_%J.err
#DSUB -R cpu=2
module load arm/bowtie2/2.5.4
module load arm/samtools/1.21
bowtie2 -p 28 -x /share/home/shli24/db/GRCm39_bt2/index115 -1 /share/home/shli24/OVX2025CUT/data/clean_data/C33_1_clean.fq.gz -2 /share/home/shli24/OVX2025CUT/data/clean_data/C33_2_clean.fq.gz |samtools sort -@ 28 -o /share/home/shli24/2025OVX/data/bam/C33.bam
