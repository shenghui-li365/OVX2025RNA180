#!/bin/bash
#DSUB -n hisat2_Y6_Liver
#DSUB -o hisat2_Y6_Liver_%J.out
#DSUB -e hisat2_Y6_Liver_%J.err
#DSUB -R cpu=2
module load arm/hisat2/2.1.0
module load arm/samtools/1.21
hisat2 --known-splicesite-infile /share/home/shli24/db/Mus_musculus_GRCm39/Ref_RNASeq/mouse.txt -x /share/home/shli24/db/Mus_musculus_GRCm39/index114/index114 --rna-strandness RF -1 /share/home/shli24/2025OVX/data/cleanData/Y6_Liver_1_clean.fq.gz -2 /share/home/shli24/2025OVX/data/cleanData/Y6_Liver_2_clean.fq.gz 2>Y6_Liver.log  | samtools view - -Sb | samtools sort - -T Y6_Liver -o /share/home/shli24/2025OVX/data/bam/Y6_Liver.bam
