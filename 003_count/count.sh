#!/bin/bash
#DSUB -n count
#DSUB -o count_%J.out
#DSUB -e count_%J.err
#DSUB -R cpu=12
#DSUB --label x86
/share/home/shli24/subread-2.1.1-Linux-x86_64/bin/featureCounts -p -t exon -g gene_id -a /share/home/shli24/db/Mus_musculus_GRCm39/Mus_musculus.GRCm39.114.gtf -o /share/home/shli24/2025OVX/data/cleanData/count.txt /share/home/shli24/2025OVX/data/bam/*.bam
