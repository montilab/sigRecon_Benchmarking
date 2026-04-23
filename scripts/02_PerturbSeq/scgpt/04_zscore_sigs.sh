#!/bin/bash -l
#$ -l h_rt=24:00:00
#$ -N zscore_scgpt
#$ -m e
#$ -j y
#$ -P montilab-p
#$ -l buyin=TRUE
#$ -pe omp 16

cd /rprojectnb2/montilab-p/projects/brcameta/projects/sig_recon/scripts/02_PerturbSeq/scgpt
module load R/4.4.3
Rscript --verbose 04_zscore_sigs.R
