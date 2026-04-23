#!/bin/bash -l
#$ -l h_rt=12:00:00
#$ -N tahoe_sigs
#$ -m e
#$ -j y
#$ -P montilab-p
#$ -pe omp 16
#$ -l mem_per_core=16G

cd /restricted/projectnb/brcameta/projects/sig_recon/scripts/06_Tahoe
module load R/4.2.1
Rscript --verbose 01_generating_sigs.R
