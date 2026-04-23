#!/bin/bash -l
#$ -l h_rt=72:00:00
#$ -N gtex_recon
#$ -m e
#$ -j y
#$ -P brcameta
#$ -pe omp 36

cd /rprojectnb2/montilab-p/projects/brcameta/projects/sig_recon/scripts/07_GTEX/network_prop
module load R/4.2.1
Rscript --verbose 02_gtex_recon.R
