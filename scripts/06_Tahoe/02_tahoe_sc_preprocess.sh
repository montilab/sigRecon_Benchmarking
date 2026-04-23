#!/bin/bash -l
#$ -l h_rt=12:00:00
#$ -N tahoe
#$ -m e
#$ -j y
#$ -P montilab-p
#$ -pe omp 28

cd /restricted/projectnb/brcameta/projects/sig_recon/scripts/06_Tahoe
conda activate rapids_singlecell
python 02_tahoe_sc_preprocess.py
