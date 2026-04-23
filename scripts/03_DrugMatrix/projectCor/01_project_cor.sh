#!/bin/bash -l
#$ -l h_rt=12:00:00
#$ -N dm_projectcor
#$ -m e
#$ -j y
#$ -P montilab-p
#$ -hard -l buyin=TRUE
#$ -pe omp 16

cd /rprojectnb2/montilab-p/projects/brcameta/projects/sig_recon/scripts/03_DrugMatrix/projectCor
module load R/4.4.3
Rscript --verbose 01_project_cor.R