#!/bin/bash -l
#$ -l h_rt=72:00:00
#$ -N dm_wgcna_recon
#$ -m e
#$ -j y
#$ -P montilab-p
#$ -hard -l buyin=TRUE
#$ -pe omp 16

cd /rprojectnb2/montilab-p/projects/brcameta/projects/sig_recon/scripts/03_DrugMatrix/network_prop
module load R/4.4.3
Rscript --verbose 02_wgcna_recon.R
