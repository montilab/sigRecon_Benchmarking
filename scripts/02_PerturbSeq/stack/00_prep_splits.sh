#!/bin/bash -l
#$ -l h_rt=12:00:00
#$ -N ps_prep_splits
#$ -m e
#$ -j y
#$ -P montilab-p
#$ -l buyin=TRUE
#$ -pe omp 16
#$ -l mem_per_core=16G

cd /rprojectnb2/montilab-p/projects/brcameta/projects/sig_recon/scripts/02_PerturbSeq/stack
conda activate stack
python 00_prep_splits.py
