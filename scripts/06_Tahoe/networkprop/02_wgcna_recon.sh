#!/bin/bash -l
#$ -l h_rt=24:00:00
#$ -N tahoe_wgcna_recon
#$ -m e
#$ -j y
#$ -P lcproject
#$ -pe omp 28

cd /rprojectnb2/montilab-p/projects/brcameta/projects/sig_recon/scripts/06_Tahoe/network_prop
module load R/4.4.3
Rscript --verbose 02_wgcna_recon.R
