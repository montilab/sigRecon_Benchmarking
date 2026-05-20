#!/bin/bash -l
#$ -l h_rt=48:00:00
#$ -N deseq_rpe1_10
#$ -m e
#$ -j y
#$ -P montilab-p
#$ -pe omp 28
#$ -t 10

cd /rprojectnb2/montilab-p/projects/brcameta/projects/sig_recon/scripts/02_PerturbSeq/stack
module load R/4.4.3
Rscript --verbose 02_deseq.R rpe1 10th $SGE_TASK_ID
