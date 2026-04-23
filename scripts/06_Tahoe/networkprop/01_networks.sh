#!/bin/bash -l
#$ -l h_rt=12:00:00
#$ -N tahoe_net
#$ -m e
#$ -j y
#$ -P brcameta
#$ -pe omp 16

cd /rprojectnb2/montilab-p/projects/brcameta/projects/sig_recon/scripts/06_Tahoe/network_prop
module load R/4.4.3
Rscript --verbose 01_networks.R
