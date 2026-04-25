#!/bin/bash
#DSUB -n fastp_S87_Spleen
#DSUB -o fastp_S87_Spleen_%J.out
#DSUB -e fastp_S87_Spleen_%J.err
#DSUB -R cpu=2
module load arm/fastp/0.23.4
fastp -i /share/home/shli24/2025OVX/data/org_data/S87_Spleen_1.fq.gz -I /share/home/shli24/2025OVX/data/org_data/S87_Spleen_2.fq.gz -o /share/home/shli24/Tissue_com/data/E2_RNA/cleanData/S87_Spleen_1_clean.fq.gz -O /share/home/shli24/Tissue_com/data/E2_RNA/cleanData/S87_Spleen_2_clean.fq.gz
