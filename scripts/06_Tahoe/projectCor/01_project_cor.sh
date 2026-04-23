#!/bin/bash -l
#$ -l h_rt=12:00:00
#$ -N tahoe_projectcor
#$ -m e
#$ -j y
#$ -P montilab-p
#$ -pe omp 16
#$ -l mem_per_core=16G

cd /rprojectnb2/montilab-p/projects/brcameta/projects/sig_recon/scripts/06_Tahoe/projectCor
module load R/4.4.3
Rscript --verbose 01_project_cor.R