#!/bin/bash -l
#$ -l h_rt=72:00:00
#$ -N deseq_mcf7
#$ -m e
#$ -j y
#$ -P montilab-p
#$ -l buyin=TRUE
#$ -pe omp 16
#$ -l mem_per_core=16G

cd /rprojectnb2/montilab-p/projects/brcameta/projects/sig_recon/scripts/04_SciPlex/stack
module load R/4.4.3
Rscript --verbose 02_deseq.R \
  --input_dir "/restricted/projectnb/agedisease/CBMrepositoryData/perturbational_data/srivatsan_2019/stack_pred_from_mcf7/" \
  --output_dir "/rprojectnb2/montilab-p/projects/brcameta/projects/sig_recon/data/stack/sciplex/" \
  --cell_line "mcf7" \
  --ctrl_name "Vehicle"
