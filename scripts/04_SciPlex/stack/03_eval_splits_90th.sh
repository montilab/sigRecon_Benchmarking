#!/bin/bash -l
#$ -l h_rt=24:00:00
#$ -N scip_stack_eval_90th
#$ -m e
#$ -j y
#$ -P el-studies
#$ -pe omp 16

cd /rprojectnb2/montilab-p/projects/brcameta/projects/sig_recon/scripts/04_SciPlex/stack
module load R/4.4.3
Rscript --verbose 03_eval_splits.R \
  --experiment "90th"
