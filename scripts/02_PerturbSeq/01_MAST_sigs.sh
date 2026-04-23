#!/bin/bash -l
#$ -l h_rt=48:00:00
#$ -N ps_sigs
#$ -m e
#$ -j y
#$ -P montilab-p
#$ -pe omp 28

cd /rprojectnb2/montilab-p/projects/brcameta/projects/sig_recon/scripts/02_PerturbSeq
module load R/4.2.1
Rscript --verbose 01_MAST_sigs.R
