#!/bin/bash -l
#$ -l h_rt=12:00:00
#$ -N drugmatrix_net
#$ -m e
#$ -j y
#$ -P montilab-p
#$ -pe omp 16

cd /rprojectnb2/montilab-p/projects/brcameta/projects/sig_recon/scripts/03_DrugMatrix/network_prop
module load R/4.4.3
Rscript --verbose 01_networks.R
