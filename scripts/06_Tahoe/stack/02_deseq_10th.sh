#!/bin/bash -l
#$ -l h_rt=48:00:00
#$ -N deseq_tahoe_10
#$ -m e
#$ -j y
#$ -P montilab-p
#$ -pe omp 28

cd /rprojectnb2/montilab-p/projects/brcameta/projects/sig_recon/scripts/06_Tahoe/stack
module load R/4.4.3
Rscript --verbose 02_deseq.R \
  --input_dir "/restricted/projectnb/agedisease/CBMrepositoryData/perturbational_data/tahoe/stack_pred/10th_perturb" \
  --output_dir "/rprojectnb2/montilab-p/projects/brcameta/projects/sig_recon/data/sigs/tahoe/stack/" \
  --experiment "10th" \
  --ctrl_name "DMSO_TF"
