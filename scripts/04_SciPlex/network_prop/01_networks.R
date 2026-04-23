library(tidyverse)
library(Seurat)
library(sigrecon)
library(doParallel)
registerDoParallel(cores = 10)

PATH     <- file.path(Sys.getenv("MLAB"), "projects/brcameta/projects/sig_recon")
savepath <- file.path(PATH, "data/sigs/sciplex")
do_save  <- TRUE

# ============================================================================
# 0. Load Sci-Plex pseudobulk objects and drug splits
# ============================================================================

a549_pb <- readRDS(file.path(PATH, "data/sci_plex/a549_filtered_pb.rds"))
k562_pb <- readRDS(file.path(PATH, "data/sci_plex/k562_filtered_pb.rds"))
mcf7_pb <- readRDS(file.path(PATH, "data/sci_plex/mcf7_filtered_pb.rds"))

drug_splits <- read.csv(file.path(savepath, "drug_splits.csv"),
                        stringsAsFactors = FALSE)

# Identify split columns
n_splits   <- sum(grepl("^split_", colnames(drug_splits)))
split_cols <- paste0("split_", seq_len(n_splits))

# ============================================================================
# 1. Standardize metadata: mark controls and per-cell-line labels
# ============================================================================

# Adjust these to your metadata column names / control label
DRUG_COL     <- "product_name"    # column in meta.data that contains drug/perturbation name
CONTROL_NAME <- "Vehicle"    # value in DRUG_COL that indicates control

# For each Seurat object, create is_control and cell_name columns
for (obj_name in c("a549_pb", "k562_pb", "mcf7_pb")) {
  obj <- get(obj_name)
  
  if (!DRUG_COL %in% colnames(obj@meta.data)) {
    stop("Column '", DRUG_COL, "' not found in meta.data of ", obj_name)
  }
  
  obj$cell_name  <- gsub("_pb$", "", obj_name)  # e.g. "a549_pb" -> "a549"
  obj$is_control <- ifelse(obj[[DRUG_COL]] == CONTROL_NAME, "Control", "Perturbed")
  
  assign(obj_name, obj)
}

# ============================================================================
# Helper: run WGCNA adjacency on a Seurat subset for a given cell line
# ============================================================================

build_wgcna_network_sciplx <- function(seurat_obj,
                                       out_path,
                                       label_prefix,
                                       min_control_samples = 10,
                                       nfeatures = 10000) {
  # Get control samples
  control_samples <- colnames(seurat_obj)[seurat_obj$is_control == "Control"]
  
  if (length(control_samples) < min_control_samples) {
    message(label_prefix, ": Not enough control samples (", length(control_samples), "), skipping")
    return(NULL)
  }
  
  # ----------------------------
  # 1a. Control-only network
  # ----------------------------
  message(label_prefix, ": Control network")
  seurat_ctrl <- seurat_obj[, control_samples]
  seurat_ctrl <- FindVariableFeatures(seurat_ctrl, nfeatures = nfeatures)
  var_genes_ctrl <- VariableFeatures(seurat_ctrl)
  seurat_ctrl <- NormalizeData(seurat_ctrl)
  
  ctrl_mat <- t(as.matrix(seurat_ctrl[var_genes_ctrl, ]@assays$RNA$data))
  ctrl_ig  <- sigrecon::wgcna.adj(ctrl_mat,
                                  cor.type = "signed hybrid",
                                  diag_zero = TRUE,
                                  beta = 6,
                                  igraph = TRUE)
  
  if (do_save) {
    dir.create(out_path, recursive = TRUE, showWarnings = FALSE)
    saveRDS(ctrl_ig, file.path(out_path, paste0(label_prefix, "_control_wgcna.rds")))
  }
  rm(ctrl_mat, ctrl_ig, seurat_ctrl)
  gc()
  
  # ----------------------------
  # 1b. 1/10th networks
  # ----------------------------
  message(label_prefix, ": 1/10th networks")
  foreach(split_num = seq_len(n_splits)) %dopar% {
    split_col <- paste0("split_", split_num)
    
    # 1/10th: drugs where split_col == TRUE
    drugs_in_split <- drug_splits$drug[drug_splits[[split_col]]]
    
    # perturbed samples whose drug is in this split
    perturbed_samples <- colnames(seurat_obj)[
      seurat_obj$is_control == "Perturbed" &
        seurat_obj$product_name %in% drugs_in_split
    ]
    
    if (length(perturbed_samples) == 0) {
      message(label_prefix, " split ", split_num, ": No perturbed samples, skipping")
      return(NULL)
    }
    
    combined_samples <- c(control_samples, perturbed_samples)
    seurat_subset <- seurat_obj[, combined_samples]
    seurat_subset <- FindVariableFeatures(seurat_subset, nfeatures = nfeatures)
    var_genes <- VariableFeatures(seurat_subset)
    seurat_subset <- NormalizeData(seurat_subset)
    
    message(label_prefix, " split ", split_num, ": ",
            length(control_samples), " control + ",
            length(perturbed_samples), " perturbed (1/10th) = ",
            length(combined_samples), " total")
    
    mat <- t(as.matrix(seurat_subset[var_genes, ]@assays$RNA$data))
    ig  <- sigrecon::wgcna.adj(mat,
                               cor.type = "signed hybrid",
                               diag_zero = TRUE,
                               igraph = TRUE)
    
    if (do_save) {
      saveRDS(ig, file.path(out_path, paste0(label_prefix, "_control_1_10th_split_", split_num, "_wgcna.rds")))
    }
    
    rm(mat, ig, seurat_subset)
    gc()
  }
  
  # ----------------------------
  # 1c. 9/10th networks
  # ----------------------------
  message(label_prefix, ": 9/10th networks")
  foreach(split_num = seq_len(n_splits)) %dopar% {
    split_col <- paste0("split_", split_num)
    
    # 9/10th: drugs where split_col == FALSE
    drugs_in_split <- drug_splits$drug[!drug_splits[[split_col]]]
    
    perturbed_samples <- colnames(seurat_obj)[
      seurat_obj$is_control == "Perturbed" &
        seurat_obj$product_name %in% drugs_in_split
    ]
    
    if (length(perturbed_samples) == 0) {
      message(label_prefix, " split ", split_num, ": No perturbed samples, skipping")
      return(NULL)
    }
    
    combined_samples <- c(control_samples, perturbed_samples)
    seurat_subset <- seurat_obj[, combined_samples]
    seurat_subset <- FindVariableFeatures(seurat_subset, nfeatures = nfeatures)
    var_genes <- VariableFeatures(seurat_subset)
    seurat_subset <- NormalizeData(seurat_subset)
    
    message(label_prefix, " split ", split_num, ": ",
            length(control_samples), " control + ",
            length(perturbed_samples), " perturbed (9/10th) = ",
            length(combined_samples), " total")
    
    mat <- t(as.matrix(seurat_subset[var_genes, ]@assays$RNA$data))
    ig  <- sigrecon::wgcna.adj(mat,
                               cor.type = "signed hybrid",
                               diag_zero = TRUE,
                               igraph = TRUE)
    
    if (do_save) {
      saveRDS(ig, file.path(out_path, paste0(label_prefix, "_control_9_10th_split_", split_num, "_wgcna.rds")))
    }
    
    rm(mat, ig, seurat_subset)
    gc()
    NULL
  }
  
  message(label_prefix, ": networks done")
  invisible(NULL)
}

# ============================================================================
# 2. Run networks for each Sci-Plex cell line
# ============================================================================

OUTPATH <- file.path(PATH, "data/wgcna_networks/sciplex")

# A549
build_wgcna_network_sciplx(
  seurat_obj   = a549_pb,
  out_path     = OUTPATH,
  label_prefix = "a549",
  min_control_samples = 1
)
gc()

# K562
build_wgcna_network_sciplx(
  seurat_obj   = k562_pb,
  out_path     = OUTPATH,
  label_prefix = "k562",
  min_control_samples = 1
)
gc()

# MCF7
build_wgcna_network_sciplx(
  seurat_obj   = mcf7_pb,
  out_path     = OUTPATH,
  label_prefix = "mcf7",
  min_control_samples = 1
)
gc()

message("All Sci-Plex networks done.")