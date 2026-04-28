#!/bin/bash -l
#$ -l h_rt=72:00:00
#$ -N deseq_k562_90
#$ -m e
#$ -j y
#$ -P montilab-p
#$ -l buyin=TRUE
#$ -pe omp 16
#$ -l mem_per_core=16G
#$ -t 1-10

cd /rprojectnb2/montilab-p/projects/brcameta/projects/sig_recon/scripts/02_PerturbSeq/stack
module load R/4.4.3
Rscript --verbose 02_deseq.R k562 90th $SGE_TASK_ID
