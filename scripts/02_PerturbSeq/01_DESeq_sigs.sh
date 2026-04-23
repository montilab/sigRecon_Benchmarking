#!/bin/bash -l
#$ -l h_rt=12:00:00
#$ -N ps_deseq
#$ -m e
#$ -j y
#$ -P brcameta
#$ -pe omp 16
#$ -l mem_per_core=16G

cd /rprojectnb2/montilab-p/projects/brcameta/projects/sig_recon/scripts/02_PerturbSeq
module load R/4.4.3
Rscript --verbose 01_DESeq_sigs.R
