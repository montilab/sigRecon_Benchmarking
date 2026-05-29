#!/bin/bash -l
#$ -l h_rt=12:00:00
#$ -N tahoe_baselines
#$ -m e
#$ -j y
#$ -P el-studies
#$ -pe omp 16

cd /rprojectnb2/montilab-p/projects/brcameta/projects/sig_recon/scripts/06_Tahoe/mean
module load R/4.4.3

Rscript --verbose tahoe_baselines.R
