#!/bin/bash -l
#$ -l h_rt=12:00:00
#$ -N sci_prep_splits
#$ -m e
#$ -j y
#$ -P montilab-p
#$ -l buyin=TRUE
#$ -pe omp 16
#$ -l mem_per_core=16G

cd /rprojectnb2/montilab-p/projects/brcameta/projects/sig_recon/scripts/04_SciPlex/stack
conda activate stack
python 002_prep_splits.py
