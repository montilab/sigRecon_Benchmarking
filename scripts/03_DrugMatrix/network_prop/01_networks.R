library(tidyverse)
library(Biobase)
library(sigrecon)
library(doParallel)

registerDoParallel(cores = 10)

PATH      <- file.path(Sys.getenv("MLAB"), "projects/brcameta/projects/sig_recon")
DATA_PATH <- file.path(Sys.getenv("AGED"), "CBMrepositoryData/perturbational_data/drugmatrix/")
savepath  <- file.path(PATH, "data/sigs/drugmatrix")
do_save   <- TRUE

# ============================================================================
# 0. Load DrugMatrix ExpressionSet objects and drug splits
# ============================================================================

liver_eset  <- readRDS(file.path(DATA_PATH, "liver.rds"))
kidney_eset <- readRDS(file.path(DATA_PATH, "kidney.rds"))
# Changing eset features to gene symbols to match signatures
liver_gene_symbols <- make.unique(fData(liver_eset)$`Gene Symbol`)
featureNames(liver_eset) <- liver_gene_symbols
kidney_gene_symbols <- make.unique(fData(kidney_eset)$`Gene Symbol`)
featureNames(kidney_eset) <- kidney_gene_symbols

drug_splits <- read.csv(file.path(savepath, "drug_splits.csv"),
                        stringsAsFactors = FALSE)

# Identify split columns
n_splits   <- sum(grepl("^split_", colnames(drug_splits)))
split_cols <- paste0("split_", seq_len(n_splits))

# ============================================================================
# 1. Standardize metadata: mark controls and per-tissue labels
# ============================================================================

# Metadata columns in pData
DRUG_COL     <- "compound:ch1"   # drug / compound name
DOSE_COL     <- "dose:ch1"       # dose column
CONTROL_DOSE <- "0 mg/kg"        # controls are dose == "0 mg/kg"

# Add is_control and tissue_name to each ExpressionSet
for (eset_name in c("liver_eset", "kidney_eset")) {
  eset <- get(eset_name)
  pd   <- pData(eset)
  
  if (!all(c(DRUG_COL, DOSE_COL) %in% colnames(pd))) {
    stop("Required columns not found in pData of ", eset_name,
         ". Needed: ", DRUG_COL, ", ", DOSE_COL)
  }
  
  pd$sample_id  <- rownames(pd)
  pd$tissue     <- gsub("_eset$", "", eset_name)  # "liver_eset" -> "liver"
  pd$is_control <- ifelse(pd[[DOSE_COL]] == CONTROL_DOSE, "Control", "Perturbed")
  
  pData(eset) <- pd
  assign(eset_name, eset)
}

# ============================================================================
# Helper: run WGCNA adjacency on an ExpressionSet (DrugMatrix) subset
# ============================================================================

build_wgcna_network_drugmatrix <- function(eset,
                                           out_path,
                                           label_prefix,
                                           min_control_samples = 10,
                                           nfeatures = 10000) {
  expr_mat <- exprs(eset)   # genes x samples
  meta     <- pData(eset)
  
  # Get control and perturbed sample IDs
  control_samples <- meta$sample_id[meta$is_control == "Control"]
  
  if (length(control_samples) < min_control_samples) {
    message(label_prefix, ": Not enough control samples (", length(control_samples), "), skipping")
    return(NULL)
  }
  
  # ----------------------------
  # 1a. Control-only network
  # ----------------------------
  message(label_prefix, ": Control network")
  
  ctrl_mat <- expr_mat[, control_samples, drop = FALSE]  # genes x control_samples
  
  # Select variable genes across control samples
  if (ncol(ctrl_mat) > 1) {
    gene_var     <- apply(ctrl_mat, 1, var, na.rm = TRUE)
    n_var_genes  <- min(nfeatures, length(gene_var))
    var_gene_ids <- names(sort(gene_var, decreasing = TRUE))[seq_len(n_var_genes)]
  } else {
    # Only one control sample; use all genes
    var_gene_ids <- rownames(ctrl_mat)
  }
  
  ctrl_mat_var <- t(ctrl_mat[var_gene_ids, , drop = FALSE])  # samples x genes
  
  ctrl_ig <- sigrecon::wgcna.adj(ctrl_mat_var,
                                 cor.type  = "signed hybrid",
                                 diag_zero = TRUE,
                                 igraph    = TRUE)
  
  if (do_save) {
    dir.create(out_path, recursive = TRUE, showWarnings = FALSE)
    saveRDS(ctrl_ig, file.path(out_path, paste0(label_prefix, "_control_wgcna.rds")))
  }
  rm(ctrl_mat, ctrl_mat_var, ctrl_ig)
  gc()
  
  # ----------------------------
  # 1b. 1/10th networks
  # ----------------------------
  message(label_prefix, ": 1/10th networks")
  foreach(split_num = seq_len(n_splits)) %dopar% {
    split_col <- paste0("split_", split_num)
    
    # 1/10th: drugs where split_col == TRUE
    drugs_in_split <- drug_splits$drug[drug_splits[[split_col]]]
    
    # perturbed samples whose compound is in this split
    perturbed_samples <- meta$sample_id[
      meta$is_control == "Perturbed" &
        meta[[DRUG_COL]] %in% drugs_in_split
    ]
    
    if (length(perturbed_samples) == 0) {
      message(label_prefix, " split ", split_num, ": No perturbed samples, skipping")
      return(NULL)
    }
    
    combined_samples <- c(control_samples, perturbed_samples)
    subset_mat <- expr_mat[, combined_samples, drop = FALSE]  # genes x combined_samples
    
    # variable genes on this subset
    if (ncol(subset_mat) > 1) {
      gene_var     <- apply(subset_mat, 1, var, na.rm = TRUE)
      n_var_genes  <- min(nfeatures, length(gene_var))
      var_gene_ids <- names(sort(gene_var, decreasing = TRUE))[seq_len(n_var_genes)]
    } else {
      var_gene_ids <- rownames(subset_mat)
    }
    
    mat <- t(subset_mat[var_gene_ids, , drop = FALSE])  # samples x genes
    
    message(label_prefix, " split ", split_num, ": ",
            length(control_samples), " control + ",
            length(perturbed_samples), " perturbed (1/10th) = ",
            length(combined_samples), " total")
    
    ig <- sigrecon::wgcna.adj(mat,
                              cor.type  = "signed hybrid",
                              diag_zero = TRUE,
                              igraph    = TRUE)
    
    if (do_save) {
      saveRDS(
        ig,
        file.path(out_path,
                  paste0(label_prefix, "_control_1_10th_split_", split_num, "_wgcna.rds"))
      )
    }
    
    rm(mat, ig, subset_mat)
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
    
    perturbed_samples <- meta$sample_id[
      meta$is_control == "Perturbed" &
        meta[[DRUG_COL]] %in% drugs_in_split
    ]
    
    if (length(perturbed_samples) == 0) {
      message(label_prefix, " split ", split_num, ": No perturbed samples, skipping")
      return(NULL)
    }
    
    combined_samples <- c(control_samples, perturbed_samples)
    subset_mat <- expr_mat[, combined_samples, drop = FALSE]  # genes x combined_samples
    
    # variable genes on this subset
    if (ncol(subset_mat) > 1) {
      gene_var     <- apply(subset_mat, 1, var, na.rm = TRUE)
      n_var_genes  <- min(nfeatures, length(gene_var))
      var_gene_ids <- names(sort(gene_var, decreasing = TRUE))[seq_len(n_var_genes)]
    } else {
      var_gene_ids <- rownames(subset_mat)
    }
    
    mat <- t(subset_mat[var_gene_ids, , drop = FALSE])  # samples x genes
    
    message(label_prefix, " split ", split_num, ": ",
            length(control_samples), " control + ",
            length(perturbed_samples), " perturbed (9/10th) = ",
            length(combined_samples), " total")
    
    ig <- sigrecon::wgcna.adj(mat,
                              cor.type  = "signed hybrid",
                              diag_zero = TRUE,
                              igraph    = TRUE)
    
    if (do_save) {
      saveRDS(
        ig,
        file.path(out_path,
                  paste0(label_prefix, "_control_9_10th_split_", split_num, "_wgcna.rds"))
      )
    }
    
    rm(mat, ig, subset_mat)
    gc()
    NULL
  }
  
  message(label_prefix, ": networks done")
  invisible(NULL)
}

# ============================================================================
# 2. Run networks for each DrugMatrix tissue
# ============================================================================

OUTPATH <- file.path(PATH, "data/wgcna_networks/drugmatrix")

# Liver
build_wgcna_network_drugmatrix(
  eset                = liver_eset,
  out_path            = OUTPATH,
  label_prefix        = "liver"
)
gc()

# Kidney
build_wgcna_network_drugmatrix(
  eset                = kidney_eset,
  out_path            = OUTPATH,
  label_prefix        = "kidney"
)
gc()

message("All DrugMatrix networks done.")