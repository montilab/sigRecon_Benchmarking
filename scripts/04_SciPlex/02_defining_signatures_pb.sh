#!/bin/bash -l
#$ -l h_rt=12:00:00
#$ -N sciplex_deseq
#$ -m e
#$ -j y
#$ -P montilab-p
#$ -hard -l buyin=TRUE
#$ -pe omp 16
#$ -l mem_per_core=16G

cd /rprojectnb2/montilab-p/projects/brcameta/sig_recon/scripts/04_SciPlex
module load R/4.4.3
Rscript --verbose 02_defining_signatures_pb.R