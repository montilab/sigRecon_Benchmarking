#!/bin/bash -l
#$ -l h_rt=24:00:00
#$ -N tahoe_prep
#$ -m e
#$ -j y
#$ -P montilab-p
#$ -pe omp 16
#$ -l mem_per_core=16G

cd /restricted/projectnb/brcameta/projects/sig_recon/scripts/06_Tahoe/stack
mamba activate rapids_singlecell
python 00_prep_splits.py
