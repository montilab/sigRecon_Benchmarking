library(tidyverse)
library(Seurat)
library(sigrecon)
library(doParallel)

registerDoParallel(cores = 10)

PATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/projects/sig_recon")
DATAPATH <- file.path(PATH, "data/perturb_seq")
do_save <- TRUE

# ============================================================================
# 0. Load pseudobulk counts & metadata
# ============================================================================

# K562
k562_counts <- read.csv(file.path(DATAPATH, "k562_processed_pb.csv"), row.names = 1)
colnames(k562_counts) <- str_replace_all(colnames(k562_counts),
                                         pattern = "\\.",
                                         replacement = "-")
k562_meta <- read.csv(file.path(DATAPATH, "k562_processed_pb_metadata.csv"), row.names = 1)
stopifnot(all.equal(colnames(k562_counts),rownames(k562_meta)))

# RPE1
rpe1_counts <- read.csv(file.path(DATAPATH, "rpe1_processed_pb.csv"), row.names = 1)
colnames(rpe1_counts) <- str_replace_all(colnames(rpe1_counts),
                                         pattern = "\\.",
                                         replacement = "-")
rpe1_meta <- read.csv(file.path(DATAPATH, "rpe1_processed_pb_metadata.csv"), row.names = 1)
rownames(rpe1_meta) <- str_replace_all(rownames(rpe1_meta),
                                       pattern = "AC118549\\.",
                                       replacement = "AC118549-")
stopifnot(all.equal(colnames(rpe1_counts),rownames(rpe1_meta)))

# ============================================================================
# 1. Convert to Seurat objects (one per cell line)
# ============================================================================

# Using raw counts as in Tahoe script
seurat_k562 <- CreateSeuratObject(counts = k562_counts, meta.data = k562_meta)
seurat_k562$cell_name <- "k562"
seurat_k562$is_control <- ifelse(seurat_k562$gene == "non-targeting", "NTC", "Perturbed")

seurat_rpe1 <- CreateSeuratObject(counts = rpe1_counts, meta.data = rpe1_meta)
seurat_rpe1$cell_name <- "rpe1"
seurat_rpe1$is_control <- ifelse(seurat_rpe1$gene == "non-targeting", "NTC", "Perturbed")

# ============================================================================
# 2. Load perturbation splits (pb_splits)
# ============================================================================

pb_splits <- read.csv(file.path(PATH, "data/sigs/perturb-seq/pb_splits.csv"),
                      stringsAsFactors = FALSE)

n_splits <- sum(grepl("^split_", colnames(pb_splits)))
split_cols <- paste0("split_", seq_len(n_splits))

# ============================================================================
# Helper: run WGCNA adjacency on a Seurat subset
# ============================================================================

build_wgcna_network <- function(seurat_obj, out_path, label_prefix,
                                experiments = c("control", "10th", "90th"),
                                min_control_samples = 10,
                                nfeatures = 10000) {
  # Get control samples
  control_samples <- colnames(seurat_obj)[seurat_obj$is_control == "NTC"]
  
  if (length(control_samples) < min_control_samples) {
    message(label_prefix, ": Not enough control samples (", length(control_samples), "), skipping")
    return(NULL)
  }
  
  # ----------------------------
  # 2a. Control-only network
  # ----------------------------
  if("control" %in% experiments) {
    message(label_prefix, ": Control network")
    seurat_ctrl <- seurat_obj[, control_samples]
    seurat_ctrl <- FindVariableFeatures(seurat_ctrl, nfeatures = nfeatures)
    var_genes_ctrl <- VariableFeatures(seurat_ctrl)
    seurat_ctrl <- NormalizeData(seurat_ctrl)
    
    ctrl_mat <- t(as.matrix(seurat_ctrl[var_genes_ctrl, ]@assays$RNA$data))
    ctrl_ig  <- sigrecon::wgcna.adj(ctrl_mat,
                                    cor.type = "signed hybrid",
                                    diag_zero = TRUE,
                                    igraph = TRUE)
    
    if (do_save) {
      dir.create(out_path, recursive = TRUE, showWarnings = FALSE)
      saveRDS(ctrl_ig, file.path(out_path, paste0(label_prefix, "_control_wgcna.rds")))
    }
    rm(ctrl_mat, ctrl_ig, seurat_ctrl)
    gc()
  }
  
  # ----------------------------
  # 2b. 1/10th networks
  # ----------------------------
  if("10th" %in% experiments) {
    message(label_prefix, ": 1/10th networks")
    foreach(split_num = seq_len(n_splits)) %dopar% {
      split_col <- paste0("split_", split_num)
      
      # 1/10th: drugs where split_col == TRUE
      pbs_in_split <- pb_splits$pb[pb_splits[[split_col]]]
      # perturbed = samples whose "gene" is one of pbs_in_split and not control
      perturbed_samples <- colnames(seurat_obj)[
        seurat_obj$is_control == "Perturbed" &
          seurat_obj$gene %in% pbs_in_split
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
  }
  
  
  # ----------------------------
  # 2c. 9/10th networks
  # ----------------------------
  if("90th" %in% experiments) {
    message(label_prefix, ": 9/10th networks")
    foreach(split_num = seq_len(n_splits)) %dopar% {
      split_col <- paste0("split_", split_num)
      
      # 9/10th: drugs where split_col == FALSE
      pbs_in_split <- pb_splits$pb[!pb_splits[[split_col]]]
      perturbed_samples <- colnames(seurat_obj)[
        seurat_obj$is_control == "Perturbed" &
          seurat_obj$gene %in% pbs_in_split
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
  }
  
  message(label_prefix, ": networks done")
  invisible(NULL)
}

# ============================================================================
# 3. Run networks for each cell line (K562 and RPE1)
# ============================================================================

# Output directory
OUTPATH <- file.path(PATH, "data/wgcna_networks/perturb-seq")

# # K562
# build_wgcna_network(
#   seurat_obj = seurat_k562,
#   out_path   = OUTPATH,
#   label_prefix = "k562",
#   experiments = c("90th"),
# )
# rm(seurat_k562_all); gc()

# RPE1
build_wgcna_network(
  seurat_obj = seurat_rpe1,
  out_path   = OUTPATH,
  label_prefix = "rpe1",
  experiments = c("90th"),
)
rm(seurat_rpe1_all); gc()

message("All perturb-seq networks done.")