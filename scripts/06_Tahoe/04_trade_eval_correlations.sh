#!/bin/bash -l
#$ -l h_rt=12:00:00
#$ -N tahoe_trade_eval_cor
#$ -m e
#$ -j y
#$ -P el-studies
#$ -pe omp 16

cd /rprojectnb2/montilab-p/projects/brcameta/projects/sig_recon/scripts/06_Tahoe
module load R/4.4.3

Rscript --verbose 04_trade_eval_correlations.R \
  --source-cell-line NCI-H23 \
  --eval-pattern "*_projectcor_*.rds"
