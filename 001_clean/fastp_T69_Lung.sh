#!/bin/bash
#DSUB -n fastp_T69_Lung
#DSUB -o fastp_T69_Lung_%J.out
#DSUB -e fastp_T69_Lung_%J.err
#DSUB -R cpu=2
module load arm/fastp/0.23.4
fastp -i /share/home/shli24/2025OVX/data/org_data/T69_Lung_1.fq.gz -I /share/home/shli24/2025OVX/data/org_data/T69_Lung_2.fq.gz -o /share/home/shli24/Tissue_com/data/E2_RNA/cleanData/T69_Lung_1_clean.fq.gz -O /share/home/shli24/Tissue_com/data/E2_RNA/cleanData/T69_Lung_2_clean.fq.gz
