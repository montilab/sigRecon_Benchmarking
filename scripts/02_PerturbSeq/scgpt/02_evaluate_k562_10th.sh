#!/bin/bash -l
#$ -l h_rt=12:00:00
#$ -N k562_10_eval
#$ -m e
#$ -j y
#$ -l gpus=1
#$ -l gpu_c=8.0
#$ -P el-studies
#$ -pe omp 16
#$ -l mem_per_core=16G

cd /rprojectnb2/montilab-p/projects/brcameta/projects/sig_recon/scripts/02_PerturbSeq/scgpt
conda activate scGPT_env
python 02_evaluate.py --data_name "k562" --split_type "10th" --runs "1,2,3,4,5,6,7,8,9,10"
