#!/bin/bash -l
#$ -l h_rt=12:00:00
#$ -N tahoe_pb
#$ -m e
#$ -j y
#$ -P brcameta
#$ -pe omp 36

cd /restricted/projectnb/brcameta/projects/sig_recon/scripts/06_Tahoe
conda activate rapids_singlecell
python 00_tahoe_pseudobulk.py
