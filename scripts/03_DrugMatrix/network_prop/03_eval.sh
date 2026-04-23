#!/bin/bash -l
#$ -l h_rt=6:00:00
#$ -N dm_net_eval
#$ -m e
#$ -j y
#$ -P montilab-p
#$ -hard -l buyin=TRUE
#$ -pe omp 16
#$ -l mem_per_core=16G

cd /rprojectnb2/montilab-p/projects/brcameta/projects/sig_recon/scripts/03_DrugMatrix/network_prop
module load R/4.4.3
Rscript --verbose 03_eval.R
