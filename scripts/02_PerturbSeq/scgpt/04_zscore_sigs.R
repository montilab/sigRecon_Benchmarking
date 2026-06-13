library(Seurat)
library(tidyverse)
library(anndata)
library(reticulate)
library(doParallel)

# Setup Environment
reticulate::use_condaenv("r-sceasy")
options(Seurat.object.assay.version = "v5")
registerDoParallel(15)

PATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/projects/sig_recon")
DATA_PATH <- file.path(PATH, "data/scgpt/prediction_results/")
SAVE_PATH <- file.path(PATH, "data/sigs/perturb-seq/scgpt/")

# Helper function to calculate signatures
calculate_signatures <- function(ad_obj, output_name) {
  
  # 1. Identify Control vs Perturbed
  # Based on your previous concat logic, 'group' contains 'ctrl_...'
  is_ctrl <- grepl("ctrl", ad_obj$obs$group)
  
  # 2. Calculate Control Distribution (Reference)
  # Convert to matrix for speed if it's not already
  ctrl_mtx <- as.matrix(ad_obj$X[is_ctrl, ])
  ctrl_means <- colMeans(ctrl_mtx)
  ctrl_sds <- apply(ctrl_mtx, 2, sd)
  
  # Add epsilon to SD to avoid division by zero for invariant genes
  ctrl_sds[ctrl_sds == 0] <- 1e-6
  
  # 3. Get unique perturbations (excluding controls)
  all_perts <- ad_obj$obs_names[!is_ctrl]
  unique_perts <- unique(all_perts)
  
  gene_names <- rownames(ad_obj$var)
  
  # 4. Parallel loop through perturbations
  sigs <- foreach(pert = unique_perts, .combine=dplyr::bind_rows, .packages = c("dplyr", "stats")) %dopar% {
    
    # Subset cells for this specific perturbation
    idx <- which(ad_obj$obs_names == pert)
    
    # Calculate Mean Profile if duplicates exist, otherwise take the single row
    if(length(idx) > 1) {
      target_profile <- colMeans(ad_obj$X[idx, , drop = FALSE])
    } else {
      target_profile <- ad_obj$X[idx, ]
    }
    
    # Calculate Z-score
    z_scores <- (target_profile - ctrl_means) / ctrl_sds
    
    # Calculate One-tailed P-value based on standard normal distribution since we only want up-regulated genes
    p_vals <- 1 - pnorm(abs(z_scores))
    
    # Create output dataframe
    data.frame(
      pert_id = pert,
      gene_symbol = gene_names,
      z_score = as.numeric(z_scores),
      p_val = as.numeric(p_vals),
      stringsAsFactors = FALSE
    )
  }
  
  return(sigs)
}

# --- Process K562 10th Data ---
# (Contains K562 perturbations and RPE1 controls)
k562_rpe1_adata <- anndata::read_h5ad(file.path(DATA_PATH, "k562_rpe1_with_rpe1_ctrl_10th.h5ad"))
k562_10th_df <- calculate_signatures(k562_rpe1_adata, "k562")
write_csv(k562_10th_df, file.path(SAVE_PATH, "k562_10th_df.csv")) # Optional CSV for easy viewing

# --- Process RPE1 10th Data ---
# (Contains RPE1 perturbations and K562 controls)
rpe1_k562_adata <- anndata::read_h5ad(file.path(DATA_PATH, "rpe1_k562_with_k562_ctrl_10th.h5ad"))
rpe1_10th_df <- calculate_signatures(rpe1_k562_adata, "rpe1")
write_csv(rpe1_10th_df, file.path(SAVE_PATH, "rpe1_10th_df.csv"))

# --- Process K562 90th Data ---
# (Contains K562 perturbations and RPE1 controls)
k562_rpe1_adata <- anndata::read_h5ad(file.path(DATA_PATH, "k562_rpe1_with_rpe1_ctrl_90th.h5ad"))
k562_90th_df <- calculate_signatures(k562_rpe1_adata, "k562")
write_csv(k562_90th_df, file.path(SAVE_PATH, "k562_90th_df.csv")) # Optional CSV for easy viewing

# --- Process RPE1 90th Data ---
# (Contains RPE1 perturbations and K562 controls)
rpe1_k562_adata <- anndata::read_h5ad(file.path(DATA_PATH, "rpe1_k562_with_k562_ctrl_90th.h5ad"))
rpe1_90th_df <- calculate_signatures(rpe1_k562_adata, "rpe1")
write_csv(rpe1_90th_df, file.path(SAVE_PATH, "rpe1_90th_df.csv"))

k562_10th_df <- read_csv(file.path(SAVE_PATH, "k562_10th_df.csv")) # Optional CSV for easy viewing
rpe1_10th_df <- read_csv(file.path(SAVE_PATH, "rpe1_10th_df.csv"))
k562_90th_df <- read_csv(file.path(SAVE_PATH, "k562_90th_df.csv")) 
rpe1_90th_df <- read_csv(file.path(SAVE_PATH, "rpe1_90th_df.csv"))

# ## Saving as named list
k562_10th_df <- k562_10th_df %>%
  dplyr::group_by(pert_id) %>%
  dplyr::mutate(adj_p = p.adjust(p_val, method="fdr")) %>%
  dplyr::ungroup()
rpe1_10th_df <- rpe1_10th_df %>%
  dplyr::group_by(pert_id) %>%
  dplyr::mutate(adj_p = p.adjust(p_val, method="fdr")) %>%
  dplyr::ungroup()
k562_90th_df <- k562_90th_df %>%
  dplyr::group_by(pert_id) %>%
  dplyr::mutate(adj_p = p.adjust(p_val, method="fdr")) %>%
  dplyr::ungroup()
rpe1_90th_df <- rpe1_90th_df %>%
  dplyr::group_by(pert_id) %>%
  dplyr::mutate(adj_p = p.adjust(p_val, method="fdr")) %>%
  dplyr::ungroup()

# scGPT seems to predict the target gene to decrease in expression, but not many genes that increase in expression post-perturbation
k562_10th_df[k562_10th_df$adj_p < 0.1,]

# When using a p-value filter of 0.05, only 9 perturbations show up with 1-2 genes each. Using a less strict filter of 0.1.
k562_10th_sigs <- sig_filter_fn(k562_10th_df, perts = unique(k562_10th_df$pert_id), pert_col = "pert_id", log2fc_col = "z_score", pval_col = "adj_p", alpha = 1, geneid_col = "gene_symbol")
rpe1_10th_sigs <- sig_filter_fn(rpe1_10th_df, perts = unique(rpe1_10th_df$pert_id), pert_col = "pert_id", log2fc_col = "z_score", pval_col = "adj_p", alpha = 1, geneid_col = "gene_symbol")
saveRDS(k562_10th_sigs, file.path(SAVE_PATH, "k562_10th_sigs.rds"))
saveRDS(rpe1_10th_sigs, file.path(SAVE_PATH, "rpe1_10th_sigs.rds"))

k562_90th_sigs <- sig_filter_fn(k562_90th_df, perts = unique(k562_90th_df$pert_id), pert_col = "pert_id", log2fc_col = "z_score", pval_col = "adj_p", alpha = 1, geneid_col = "gene_symbol")
rpe1_90th_sigs <- sig_filter_fn(rpe1_90th_df, perts = unique(rpe1_90th_df$pert_id), pert_col = "pert_id", log2fc_col = "z_score", pval_col = "adj_p", alpha = 1, geneid_col = "gene_symbol")
saveRDS(k562_90th_sigs, file.path(SAVE_PATH, "k562_90th_sigs.rds"))
saveRDS(rpe1_90th_sigs, file.path(SAVE_PATH, "rpe1_90th_sigs.rds"))

