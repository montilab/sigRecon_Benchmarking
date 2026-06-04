#!/bin/bash -l
#$ -l h_rt=24:00:00
#$ -N ps_orthos
#$ -m e
#$ -j y
#$ -P montilab-p
#$ -hard -l buyin=TRUE
#$ -pe omp 16
#$ -l mem_per_core=16G

cd /rprojectnb2/montilab-p/projects/brcameta/projects/sig_recon/scripts/02_PerturbSeq/orthos
module load R/4.4.3
Rscript --verbose 01_orthos_sigs.R
