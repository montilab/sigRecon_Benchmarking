#!/bin/bash -l
#$ -l h_rt=12:00:00
#$ -N ps_net_eval
#$ -m e
#$ -j y
#$ -P montilab-p
#$ -hard -l buyin=TRUE
#$ -pe omp 16

cd /rprojectnb2/montilab-p/projects/brcameta/projects/sig_recon/scripts/02_PerturbSeq/network_prop
module load R/4.4.3
Rscript --verbose 03_eval.R
