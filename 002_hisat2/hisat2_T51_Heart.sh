#!/bin/bash
#DSUB -n hisat2_T51_Heart
#DSUB -o hisat2_T51_Heart_%J.out
#DSUB -e hisat2_T51_Heart_%J.err
#DSUB -R cpu=2
module load arm/hisat2/2.1.0
module load arm/samtools/1.21
hisat2 --known-splicesite-infile /share/home/shli24/db/Mus_musculus_GRCm39/Ref_RNASeq/mouse.txt -x /share/home/shli24/db/Mus_musculus_GRCm39/index114/index114 --rna-strandness RF -1 /share/home/shli24/2025OVX/data/cleanData/T51_Heart_1_clean.fq.gz -2 /share/home/shli24/2025OVX/data/cleanData/T51_Heart_2_clean.fq.gz 2>T51_Heart.log  | samtools view - -Sb | samtools sort - -T T51_Heart -o /share/home/shli24/2025OVX/data/bam/T51_Heart.bam
