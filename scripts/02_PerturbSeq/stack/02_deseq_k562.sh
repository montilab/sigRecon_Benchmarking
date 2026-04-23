#!/bin/bash -l
#$ -l h_rt=72:00:00
#$ -N stack_deseq_k562
#$ -m e
#$ -j y
#$ -P montilab-p
#$ -l buyin=TRUE
#$ -pe omp 16
#$ -l mem_per_core=16G

cd /rprojectnb2/montilab-p/projects/brcameta/projects/sig_recon/scripts/02_PerturbSeq/stack
module load R/4.4.3
Rscript --verbose 02_deseq_k562.R
