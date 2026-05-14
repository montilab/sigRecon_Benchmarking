#!/bin/bash -l
#$ -l h_rt=48:00:00
#$ -N deseq_k562_90
#$ -m e
#$ -j y
#$ -P montilab-p
#$ -l buyin=TRUE
#$ -pe omp 16
#$ -l mem_per_core=16G
#$ -t 1-10

# SELECTED_SPLITS=(2 5 8)
# # Use the task ID (starting at 1) as an index to pick a specific number
# SPLIT=${SELECTED_SPLITS[$((SGE_TASK_ID-1))]}

cd /rprojectnb2/montilab-p/projects/brcameta/projects/sig_recon/scripts/02_PerturbSeq/stack
module load R/4.4.3
Rscript --verbose 02_deseq.R k562 90th $SGE_TASK_ID
