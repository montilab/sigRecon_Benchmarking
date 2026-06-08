#!/bin/bash -l
#$ -l h_rt=12:00:00
#$ -N perturb_seq_trade_eval_cor
#$ -m e
#$ -j y
#$ -P findthecause
#$ -pe omp 16

cd /rprojectnb2/montilab-p/projects/brcameta/projects/sig_recon/scripts/02_PerturbSeq
module load R/4.4.3

Rscript --verbose 02_trade_eval_correlations.R \
  --method projectcor_gsva
