#!/bin/bash -l
#$ -l h_rt=8:00:00
#$ -N scip_proj_cor
#$ -m e
#$ -j y
#$ -P montilab-p
#$ -hard -l buyin=TRUE
#$ -pe omp 16

cd /rprojectnb2/montilab-p/projects/brcameta/projects/sig_recon/scripts/04_SciPlex/projectCor
module load R/4.4.3
Rscript --verbose 01_project_cor.R
