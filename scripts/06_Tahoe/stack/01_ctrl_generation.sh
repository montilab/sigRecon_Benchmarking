#!/bin/bash -l
#$ -l h_rt=24:00:00
#$ -N ctrl_gen
#$ -m e
#$ -j y
#$ -l gpus=1
#$ -l gpu_c=8.0
#$ -P brcameta
#$ -pe omp 16
#$ -l mem_per_core=16G

cd /restricted/projectnb/brcameta/projects/sig_recon/
module load miniconda
conda activate stack

# Directory containing the data files
DATA_DIR="/restricted/projectnb/agedisease/CBMrepositoryData/perturbational_data/tahoe/stack_subsets/ctrl_comparison"

# Create array of unique cell line names
CELL_ARRAY=($(ls ${DATA_DIR}/*_maxconc_dmso.h5ad | xargs -n1 basename | sed 's/_maxconc_dmso.h5ad//' | sort -u))

# Get the cell line name for this task
# CELLNAME=${CELL_ARRAY[$((SGE_TASK_ID-1))]}
CELLNAME="NCI-H23"

pbs=("Idarubicin hydrochloride" "Irinotecan hydrochloride" "Lenalidomide hemihydrate"
 "Lidocaine hydrochloride" "Lumateperone tosylate" "Mitoxantrone dihydrochloride" 
 "Neratinib maleate" "Norepinephrine hydrochloride" "Palmatine chloride"
 "Pasireotide acetate" "Pentamidine isethionate" "Pravastatin sodium"
 "Pyridoxine hydrochloride" "R-Verapamil hydrochloride" "Ropivacaine hydrochloride monohydrate" 
 "S-Adenosyl-L-methionine disulfate tosylate" "Temsirolimus" "Terfenadine" 
 "Thymol" "Tofacitinib" "Tofacitinib citrate"
 "Tolcapone" "Topotecan hydrochloride" "Trametinib" 
 "Triamcinolone" "Triclosan" "Trifluridine" 
 "Trimetrexate" "Tucidinostat" "Vilanterol"
 "Vinblastine sulfate")

echo "Processing ${CELLNAME}"

# Run stack-generation with the current cell line
stack-generation --checkpoint "scripts/stack/notebooks/tutorial-pred-model/bc_large_aligned.ckpt" \
--base-adata "${DATA_DIR}/${CELLNAME}_maxconc_dmso.h5ad" \
--test-adata "${DATA_DIR}/${CELLNAME}_other_dmso.h5ad" \
--genelist "scripts/stack/notebooks/tutorial-pred-model/basecount_1000per_15000max.pkl" \
--split-column drug \
--split-values "${pbs[@]}" \
--batch-size 16 \
--num-steps 5 \
--output-dir "/restricted/projectnb/agedisease/CBMrepositoryData/perturbational_data/tahoe/stack_pred/ctrl_generation/${CELLNAME}" \
--prompt-ratio 0.25 \
--context-ratio 0.4