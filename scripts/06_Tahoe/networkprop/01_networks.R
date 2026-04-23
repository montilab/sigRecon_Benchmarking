library(tidyverse)
library(Seurat)
library(sigrecon)
library(doParallel)
registerDoParallel(cores=10)

PATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/projects/sig_recon")
DATAPATH <- file.path(Sys.getenv("AGED"), "CBMrepositoryData/perturbational_data/tahoe")
do_save <- FALSE

seurat_obj_f <- readRDS(file.path(PATH, "data/tahoe/pseudobulk/merged_pseudobulk_filtered.rds"))

# 1. Subset to pre-selected cell lines and drugs
celllines <- read.csv(file.path(PATH, "data/sigs/tahoe/cell_lines.csv"))$x
drug_splits <- read.csv(file.path(PATH, "data/sigs/tahoe/drug_splits.csv"))
drug_subset <- c(drug_splits$drug, "DMSO_TF")
seurat_obj_f_sub <- subset(seurat_obj_f, 
                           (cell_name %in% celllines) & (drug_name %in% drug_subset))

# 2. Learning Network on Controls
print("Control Networks")
foreach(cell_line = celllines) %dopar% {
  
  # Subsetting to cell line controls
  seurat_cellline <- subset(seurat_obj_f_sub,
                            (cell_name == cell_line))
  drug_samples <- colnames(seurat_cellline)[seurat_cellline$drug_name == "DMSO_TF"]
  seurat_cellline <- seurat_cellline[,drug_samples]
  
  seurat_cellline <- FindVariableFeatures(seurat_cellline, nfeatures = 10000)
  cellline_var_feats <- VariableFeatures(seurat_cellline)
  seurat_cellline <- NormalizeData(seurat_cellline)
  
  if (length(drug_samples) < 10) {
    message(paste0("Not enough samples, skipping ", cell_line))
    return(NULL)
  }
  
  cellline_mat <- t(as.matrix(seurat_cellline[cellline_var_feats,]@assays$RNA$data))
  cellline_ig <- sigrecon::wgcna.adj(cellline_mat, cor.type = "signed hybrid", diag_zero = TRUE, igraph = TRUE)
  saveRDS(cellline_ig, file.path(PATH, paste0("data/wgcna_networks/tahoe/", cell_line, "_control_wgcna.rds")))
  rm(cellline_mat)
  rm(cellline_ig)
  gc()
}
print("Control Networks Done")

print("1/10th networks")
# 3. Learning Network on Control + 1/10th of perturbed
foreach(cell_line = celllines) %dopar% {

  # Subset to cell line
  seurat_cellline <- subset(seurat_obj_f_sub,
                            (cell_name == cell_line))

  # Get control samples
  control_samples <- colnames(seurat_cellline)[seurat_cellline$drug_name == "DMSO_TF"]

  if (length(control_samples) < 10) {
    message(paste0("Not enough control samples, skipping ", cell_line))
    return(NULL)
  }

  # Loop through each split
  for (split_num in 1:10) {
    split_col <- paste0("split_", split_num)

    # Get drugs that are TRUE for this split (1/10th subset)
    drugs_in_split <- drug_splits$drug[drug_splits[[split_col]] == TRUE]

    # Get perturbed samples for drugs in this split
    perturbed_subset <- colnames(seurat_cellline)[seurat_cellline$drug_name %in% drugs_in_split]

    if (length(perturbed_subset) == 0) {
      message(paste0(cell_line, " split ", split_num, ": No perturbed samples, skipping"))
      next
    }

    # Combine control and 1/10th perturbed samples
    combined_samples <- c(control_samples, perturbed_subset)

    # Subset to combined samples BEFORE finding variable features
    seurat_subset <- seurat_cellline[, combined_samples]

    # Find variable features on the subsetted data
    seurat_subset <- FindVariableFeatures(seurat_subset, nfeatures = 10000)
    cellline_var_feats <- VariableFeatures(seurat_subset)
    seurat_subset <- NormalizeData(seurat_subset)
    
    message(paste0(cell_line, " split ", split_num, ": ",
                   length(control_samples), " control + ",
                   length(perturbed_subset), " perturbed (1/10th) = ",
                   length(combined_samples), " total"))

    cellline_mat <- t(as.matrix(seurat_subset[cellline_var_feats, ]@assays$RNA$data))
    cellline_ig <- sigrecon::wgcna.adj(cellline_mat, cor.type = "signed hybrid", diag_zero = TRUE, igraph = TRUE)
    saveRDS(cellline_ig, file.path(PATH, paste0("data/wgcna_networks/tahoe/",
                                                cell_line, "_control_1_10th_split_", split_num, "_wgcna.rds")))

    rm(cellline_mat)
    rm(cellline_ig)
    rm(seurat_subset)
    gc()
  }
}
print("1/10th networks done")

print("9/10th networks")
# 4. Learning Network on Control + 9/10th of perturbed
foreach(cell_line = celllines) %dopar% {

  # Subset to cell line
  seurat_cellline <- subset(seurat_obj_f_sub,
                            (cell_name == cell_line))

  # Get control samples
  control_samples <- colnames(seurat_cellline)[seurat_cellline$drug_name == "DMSO_TF"]

  if (length(control_samples) < 10) {
    message(paste0("Not enough control samples, skipping ", cell_line))
    return(NULL)
  }

  # Loop through each split
  for (split_num in 1:10) {
    split_col <- paste0("split_", split_num)

    # Get drugs that are FALSE for this split (9/10th subset)
    drugs_in_split <- drug_splits$drug[drug_splits[[split_col]] == FALSE]

    # Get perturbed samples for drugs in this split
    perturbed_subset <- colnames(seurat_cellline)[seurat_cellline$drug_name %in% drugs_in_split]

    if (length(perturbed_subset) == 0) {
      message(paste0(cell_line, " split ", split_num, ": No perturbed samples, skipping"))
      next
    }

    # Combine control and 9/10th perturbed samples
    combined_samples <- c(control_samples, perturbed_subset)

    # Subset to combined samples BEFORE finding variable features
    seurat_subset <- seurat_cellline[, combined_samples]

    # Find variable features on the subsetted data
    seurat_subset <- FindVariableFeatures(seurat_subset, nfeatures = 10000)
    cellline_var_feats <- VariableFeatures(seurat_subset)
    seurat_subset <- NormalizeData(seurat_subset)
    
    message(paste0(cell_line, " split ", split_num, ": ",
                   length(control_samples), " control + ",
                   length(perturbed_subset), " perturbed (9/10th) = ",
                   length(combined_samples), " total"))

    cellline_mat <- t(as.matrix(seurat_subset[cellline_var_feats, ]@assays$RNA$data))
    cellline_ig <- sigrecon::wgcna.adj(cellline_mat, cor.type = "signed hybrid", diag_zero = TRUE, igraph = TRUE)
    saveRDS(cellline_ig, file.path(PATH, paste0("data/wgcna_networks/tahoe/",
                                                cell_line, "_control_9_10th_split_", split_num, "_wgcna.rds")))

    rm(cellline_mat)
    rm(cellline_ig)
    rm(seurat_subset)
    gc()
  }
}
print("9/10th networks done")