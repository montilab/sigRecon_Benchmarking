#!/bin/bash -l
#$ -l h_rt=12:00:00
#$ -N sciplex_sigs
#$ -o sciplex_sigs.log
#$ -m e
#$ -j y
#$ -P brcameta
#$ -pe omp 16

cd /rprojectnb2/montilab-p/projects/brcameta/sig_recon/scripts/04_SciPlex
module load R/4.1.2
Rscript --verbose 01_defining_signatures.R