#!/bin/bash -l
#$ -l h_rt=96:00:00
#$ -N deseq_rpe1_ctrl
#$ -m e
#$ -j y
#$ -P montilab-p
#$ -pe omp 28

cd /rprojectnb2/montilab-p/projects/brcameta/projects/sig_recon/scripts/02_PerturbSeq/stack
module load R/4.4.3
Rscript --verbose 02_deseq_ctrl.R rpe1 
