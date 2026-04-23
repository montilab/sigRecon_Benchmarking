#!/bin/bash -l
#$ -l h_rt=12:00:00
#$ -N ps_projectcor
#$ -m e
#$ -j y
#$ -P montilab-p
#$ -pe omp 28

cd /rprojectnb2/montilab-p/projects/brcameta/projects/sig_recon/scripts/02_PerturbSeq/projectCor
module load R/4.4.3
Rscript --verbose 01_project_cor.R