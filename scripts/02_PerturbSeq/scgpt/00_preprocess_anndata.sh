#!/bin/bash -l
#$ -l h_rt=48:00:00
#$ -N process_replogle
#$ -m e
#$ -j y
#$ -P brcameta
#$ -pe omp 16
#$ -l mem_per_core=16G

cd /rprojectnb2/montilab-p/projects/brcameta/projects/sig_recon/scripts/02_PerturbSeq/scgpt
conda activate scGPT_env
python 00_preprocess_anndata.py
