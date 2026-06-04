#!/bin/bash -l
#$ -l h_rt=12:00:00
#$ -N ps_orthos_eval
#$ -m e
#$ -j y
#$ -P montilab-p
#$ -pe omp 16

cd /rprojectnb2/montilab-p/projects/brcameta/projects/sig_recon/scripts/02_PerturbSeq/orthos
module load R/4.4.3
Rscript --verbose 02_eval.R
