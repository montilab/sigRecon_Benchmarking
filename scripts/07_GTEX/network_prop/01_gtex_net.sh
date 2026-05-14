#!/bin/bash -l
#$ -l h_rt=24:00:00
#$ -N gtex_net
#$ -m e
#$ -j y
#$ -P el-studies
#$ -pe omp 16

cd /rprojectnb2/montilab-p/projects/brcameta/projects/sig_recon/scripts
module load R/4.2.1
Rscript --verbose 07_GTEX/network_prop/01_gtex_net.R
