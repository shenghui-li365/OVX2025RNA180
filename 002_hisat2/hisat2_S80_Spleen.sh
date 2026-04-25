#!/bin/bash
#DSUB -n hisat2_S80_Spleen
#DSUB -o hisat2_S80_Spleen_%J.out
#DSUB -e hisat2_S80_Spleen_%J.err
#DSUB -R cpu=2
module load arm/hisat2/2.1.0
module load arm/samtools/1.21
hisat2 --known-splicesite-infile /share/home/shli24/db/Mus_musculus_GRCm39/Ref_RNASeq/mouse.txt -x /share/home/shli24/db/Mus_musculus_GRCm39/index114/index114 --rna-strandness RF -1 /share/home/shli24/2025OVX/data/cleanData/S80_Spleen_1_clean.fq.gz -2 /share/home/shli24/2025OVX/data/cleanData/S80_Spleen_2_clean.fq.gz 2>S80_Spleen.log  | samtools view - -Sb | samtools sort - -T S80_Spleen -o /share/home/shli24/2025OVX/data/bam/S80_Spleen.bam
