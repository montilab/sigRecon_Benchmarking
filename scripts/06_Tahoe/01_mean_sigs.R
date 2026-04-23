library(tidyverse)
library(Seurat)
library(anndata)
library(reticulate)
library(DESeq2)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(sigrecon)
library(doParallel)
registerDoParallel(cores=10)

PATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/projects/sig_recon")
savepath <- file.path(PATH, "data/sigs/tahoe")
do_save <- TRUE

# ============================================================================
# Helper Functions
# ============================================================================

# Map Ensembl IDs to gene symbols
map_ids <- function(genes) {
  is_ensembl <- grepl("^ENSG", genes)
  
  if (any(is_ensembl)) {
    ensembl_ids <- genes[is_ensembl]
    tryCatch(
      {
        symbols <- mapIds(org.Hs.eg.db,
                          keys = ensembl_ids,
                          column = "SYMBOL",
                          keytype = "ENSEMBL",
                          multiVals = "first")
        symbols[is.na(symbols)] <- ensembl_ids[is.na(symbols)]
        genes[is_ensembl] <- symbols
        message(paste("Mapped", length(ensembl_ids), "Ensembl IDs to gene symbols."))
      }, error = function(e) {
        warning("An error occurred: ", conditionMessage(e), call. = FALSE)
      }
    )
  }
  return(genes)
}

# Run DESeq2 analysis on Seurat object
run_deseq2_tahoe <- function(seurat_obj, sample_ids_perturbed, sample_ids_control) {
  meta <- seurat_obj@meta.data
  all_sample_ids <- c(sample_ids_perturbed, sample_ids_control)
  sub_count <- seurat_obj[, all_sample_ids]@assays$RNA$counts
  
  sub_meta <- meta[colnames(sub_count), ]
  stopifnot(all.equal(rownames(sub_meta), colnames(sub_count)))
  sub_meta$condition <- ifelse(rownames(sub_meta) %in% sample_ids_perturbed, 
                               "treatment", "control")
  sub_meta$condition <- as.factor(sub_meta$condition)
  sub_meta$plate <- as.factor(sub_meta$plate)
  
  
  dds <- DESeqDataSetFromMatrix(
    countData = sub_count,
    colData = sub_meta,
    design = ~ plate + condition
  )
  
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, contrast = c("condition", "treatment", "control"))
  
  return(as.data.frame(res))
}

# ============================================================================
# 0. Load Data and Create Drug Splits
# ============================================================================

seurat_obj_f <- readRDS(file.path(PATH, "data/tahoe/pseudobulk/merged_pseudobulk_filtered.rds"))

set.seed(123)  # For reproducibility

# Get cell lines and drugs
celllines <- seurat_obj_f$cell_name %>% unique()
# Ensure we only have the target cell lines
celllines <- sample(celllines, size = 10)
write.csv(celllines, file.path(savepath, "cell_lines.csv"), row.names = FALSE)

drugs <- seurat_obj_f$drug_name %>% unique()
drugs <- drugs[drugs != "DMSO_TF"]  # Exclude control

cat("Total cell lines after filtering:", length(celllines), "\n")
cat("Cell lines:", paste(celllines, collapse = ", "), "\n")
cat("Total drugs (excluding DMSO):", length(drugs), "\n")

# Find drugs present in all cell lines with sufficient replicates
drug_counts_by_cellline <- seurat_obj_f@meta.data %>%
  dplyr::filter(cell_name %in% celllines) %>%
  dplyr::filter(drug_name != "DMSO_TF") %>%
  dplyr::group_by(cell_name, drug_name) %>%
  dplyr::summarise(n_samples = n(), .groups = "drop") %>%
  dplyr::filter(n_samples >= 2)  # At least 2 replicates

shared_drugs <- drug_counts_by_cellline %>%
  group_by(drug_name) %>%
  summarise(n_celllines = n_distinct(cell_name), .groups = "drop") %>%
  filter(n_celllines == length(celllines)) %>%
  pull(drug_name)

cat("Shared drugs across all", length(celllines), "cell lines:", length(shared_drugs), "\n")

# Create 10 random splits of drugs
n_splits <- 10
drugs_shuffled <- sample(shared_drugs)
split_size <- ceiling(length(shared_drugs) / n_splits)

# Create splits dataframe
drug_splits <- data.frame(
  drug = shared_drugs,
  stringsAsFactors = FALSE
)

# Initialize split columns
for (i in 1:n_splits) {
  split_col <- paste0("split_", i)
  drug_splits[[split_col]] <- FALSE
}

# Assign drugs to their primary split
for (i in 1:n_splits) {
  start_idx <- (i - 1) * split_size + 1
  end_idx <- min(i * split_size, length(drugs_shuffled))
  drugs_in_split <- drugs_shuffled[start_idx:end_idx]
  
  drug_splits[[paste0("split_", i)]][drug_splits$drug %in% drugs_in_split] <- TRUE
}

# Save drug splits
if (do_save) {
  write.csv(drug_splits, 
            file.path(savepath, "drug_splits.csv"),
            row.names = FALSE)
}

cat("Drug splits created:\n")
print(colSums(drug_splits[, -1]))  # Show number of drugs per split

# # ============================================================================
# # Analysis 1: Control - Each Cell Line vs All Others
# # ============================================================================
# 
# cat("\n=== Analysis 1: Control vs Control (Each Cell Line vs All Others) ===\n")
# 
# Get control samples for each cell line
control_samples_by_cellline <- list()
for (cell_line in celllines) {
  ctrl_samples <- rownames(seurat_obj_f@meta.data)[
    seurat_obj_f@meta.data$cell_name == cell_line &
      seurat_obj_f@meta.data$drug_name == "DMSO_TF"
  ]
  control_samples_by_cellline[[cell_line]] <- ctrl_samples
  cat(cell_line, "controls:", length(ctrl_samples), "\n")
}
# 
# # Find shared genes across all cell lines
# shared_genes <- rownames(seurat_obj_f@assays$RNA$counts)
# cat("Total genes:", length(shared_genes), "\n")
# 
# # For each cell line, compare against all others using foreach
celllines <- as.character(celllines)
# 
# ctrl_results <- foreach(target_cell = celllines, .packages = c("DESeq2", "tidyverse", "org.Hs.eg.db", "AnnotationDbi", "sigrecon")) %dopar% {
#   
#   cat("\nProcessing", target_cell, "vs all other cell lines\n")
#   
#   # Get control samples for target cell line
#   target_ctrl <- control_samples_by_cellline[[target_cell]]
#   
#   # Get control samples for all other cell lines
#   other_celllines <- setdiff(celllines, target_cell)
#   other_ctrl <- unlist(control_samples_by_cellline[other_celllines])
#   
#   if (length(target_ctrl) < 2 || length(other_ctrl) < 2) {
#     cat("  Skipping: insufficient controls\n")
#     return(NULL)
#   }
#   
#   cat("  Target controls:", length(target_ctrl), 
#       "Other controls:", length(other_ctrl), "\n")
#   
#   # Combine counts and metadata
#   all_samples <- c(target_ctrl, other_ctrl)
#   combined_counts <- seurat_obj_f@assays$RNA$counts[, all_samples]
#   combined_meta <- seurat_obj_f@meta.data[all_samples, ]
#   combined_meta$cell_line_group <- ifelse(rownames(combined_meta) %in% target_ctrl, 
#                                           "target", "other")
#   combined_meta$cell_line_group <- as.factor(combined_meta$cell_line_group)
#   combined_meta$plate <- as.factor(combined_meta$plate)
#   
#   # Run DESeq2
#   dds <- DESeqDataSetFromMatrix(
#     countData = combined_counts,
#     colData = combined_meta,
#     design = ~ plate + cell_line_group
#   )
#   
#   dds <- DESeq(dds, quiet = TRUE)
#   
#   # Target vs Others
#   res_target_vs_others <- results(dds, contrast = c("cell_line_group", "target", "other"))
#   
#   # Clean cell line name
#   target_clean <- str_replace_all(target_cell, c("-" = "_", " " = "", "/" = "_"))
#   
#   pb_name_target <- paste0(target_clean, "_ctrl")
#   
#   res_target_vs_others_df <- as.data.frame(res_target_vs_others) %>%
#     dplyr::mutate(pb = pb_name_target) %>%
#     as_tibble(rownames = "gene") %>%
#     mutate(gene = map_ids(gene)) %>%
#     arrange(padj, desc(abs(log2FoldChange)))
#   
#   sig_target_vs_others <- sig_filter_fn(
#     res_target_vs_others_df, 
#     perts = pb_name_target, 
#     pert_col = "pb", 
#     log2fc_col = "log2FoldChange", 
#     pval_col = "padj", 
#     geneid_col = "gene"
#   )
#   
#   # Save
#   if (do_save) {
#     write.csv(res_target_vs_others_df,
#               file.path(savepath, "mean/ctrl", paste0(target_clean, "_vs_others.csv")),
#               row.names = FALSE)
#     
#     saveRDS(sig_target_vs_others, 
#             file.path(savepath, "mean/ctrl", paste0(target_clean, "_ctrl.rds")))
#   }
#   
#   if (!is.null(sig_target_vs_others[[pb_name_target]])) {
#     cat("  ", target_cell, "vs others - Up:", length(sig_target_vs_others[[pb_name_target]]$up), "\n")
#   }
#   
#   return(list(
#     cell_line = target_cell,
#     results_df = res_target_vs_others_df,
#     signatures = sig_target_vs_others
#   ))
# }

# ============================================================================
# Analysis 2: 1/10th Perturbed vs All Controls (10 splits)
# ============================================================================

cat("\n=== Analysis 2: 1/10th Perturbed vs Controls ===\n")

# Use nested foreach for cell lines and splits
results_1_10th <- foreach(cell_line = celllines, .packages = c("DESeq2", "tidyverse", "org.Hs.eg.db", "AnnotationDbi", "sigrecon")) %:%
  foreach(i = 1:n_splits) %do% {
    
    cat("Processing", cell_line, "split", i, "of", n_splits, "\n")
    
    cell_line_clean <- str_replace_all(cell_line, c("-" = "_", " " = "", "/" = "_"))
    ctrl_samples <- control_samples_by_cellline[[cell_line]]
    
    # Get drugs in this split
    drugs_in_split <- drug_splits$drug[drug_splits[[paste0("split_", i)]]]
    
    # Get drug samples for this cell line
    drug_samples <- rownames(seurat_obj_f@meta.data)[
      seurat_obj_f@meta.data$cell_name == cell_line & 
        seurat_obj_f@meta.data$drug_name %in% drugs_in_split
    ]
    
    if (length(drug_samples) < 3) {
      cat("  Skipping: not enough drug samples\n")
      return(NULL)
    }
    
    cat("  Drugs:", length(drugs_in_split), 
        "Drug samples:", length(drug_samples), 
        "Control samples:", length(ctrl_samples), "\n")
    
    # Run DESeq2
    res <- run_deseq2_tahoe(
      seurat_obj = seurat_obj_f,
      sample_ids_perturbed = drug_samples,
      sample_ids_control = ctrl_samples
    )
    
    pb_name <- paste0(cell_line_clean, "_split_", i)
    
    res_df <- res %>%
      dplyr::mutate(pb = pb_name) %>%
      as_tibble(rownames = "gene") %>%
      mutate(gene = map_ids(gene)) %>%
      arrange(padj, desc(abs(log2FoldChange)))
    
    # Extract signatures
    sig <- sig_filter_fn(
      res_df, 
      perts = pb_name, 
      pert_col = "pb", 
      log2fc_col = "log2FoldChange", 
      pval_col = "padj", 
      geneid_col = "gene"
    )
    
    # Save full results
    if (do_save) {
      write.csv(res_df,
                file.path(savepath, "mean/10th_perturb", 
                          paste0(cell_line_clean, "_10per_", i, ".csv")),
                row.names = FALSE)
    }
    
    if (!is.null(sig[[pb_name]])) {
      cat("  Up:", length(sig[[pb_name]]$up), "\n")
    }
    
    return(list(
      cell_line = cell_line,
      split = i,
      pb_name = pb_name,
      results_df = res_df,
      signature = sig[[pb_name]]
    ))
  }

# Reorganize results by cell line and save
for (idx in 1:length(celllines)) {
  cell_line <- celllines[idx]
  cell_line_clean <- str_replace_all(cell_line, c("-" = "_", " " = "", "/" = "_"))
  
  # Extract signatures for this cell line
  sig_list_1_10th <- list()
  for (split_idx in 1:n_splits) {
    result <- results_1_10th[[idx]][[split_idx]]
    if (!is.null(result)) {
      sig_list_1_10th[[paste0("split_", split_idx)]] <- result$signature
    }
  }
  
  # Save signatures for this cell line
  if (do_save && length(sig_list_1_10th) > 0) {
    saveRDS(sig_list_1_10th,
            file.path(savepath, "mean/10th_perturb", paste0(cell_line_clean, "_10_signatures.rds")))
  }
}

# ============================================================================
# Analysis 3: 9/10th Perturbed vs All Controls (10 splits)
# ============================================================================

cat("\n=== Analysis 3: 9/10th Perturbed vs Controls ===\n")

# Use nested foreach for cell lines and splits
results_9_10th <- foreach(cell_line = celllines, .packages = c("DESeq2", "tidyverse", "org.Hs.eg.db", "AnnotationDbi", "sigrecon")) %:%
  foreach(i = 1:n_splits) %dopar% {
    
    cat("Processing", cell_line, "split", i, "of", n_splits, "(excluding split", i, ")\n")
    
    cell_line_clean <- str_replace_all(cell_line, c("-" = "_", " " = "", "/" = "_"))
    ctrl_samples <- control_samples_by_cellline[[cell_line]]
    
    # Get drugs NOT in this split (9/10th)
    drugs_in_9_10th <- drug_splits$drug[!drug_splits[[paste0("split_", i)]]]
    
    # Get drug samples for this cell line
    drug_samples <- rownames(seurat_obj_f@meta.data)[
      seurat_obj_f@meta.data$cell_name == cell_line & 
        seurat_obj_f@meta.data$drug_name %in% drugs_in_9_10th
    ]
    
    if (length(drug_samples) < 3) {
      cat("  Skipping: not enough drug samples\n")
      return(NULL)
    }
    
    cat("  Drugs:", length(drugs_in_9_10th),
        "Drug samples:", length(drug_samples), 
        "Control samples:", length(ctrl_samples), "\n")
    
    # Run DESeq2
    res <- run_deseq2_tahoe(
      seurat_obj = seurat_obj_f,
      sample_ids_perturbed = drug_samples,
      sample_ids_control = ctrl_samples
    )
    
    pb_name <- paste0(cell_line_clean, "_split_", i, "_excluded")
    
    res_df <- res %>%
      dplyr::mutate(pb = pb_name) %>%
      as_tibble(rownames = "gene") %>%
      mutate(gene = map_ids(gene)) %>%
      arrange(padj, desc(abs(log2FoldChange)))
    
    # Extract signatures
    sig <- sig_filter_fn(
      res_df, 
      perts = pb_name, 
      pert_col = "pb", 
      log2fc_col = "log2FoldChange", 
      pval_col = "padj", 
      geneid_col = "gene"
    )
    
    # Save full results
    if (do_save) {
      write.csv(res_df,
                file.path(savepath, "mean/90th_perturb", 
                          paste0(cell_line_clean, "_90_", i, ".csv")),
                row.names = FALSE)
    }
    
    if (!is.null(sig[[pb_name]])) {
      cat("  Up:", length(sig[[pb_name]]$up), "\n")
    }
    
    return(list(
      cell_line = cell_line,
      split = i,
      pb_name = pb_name,
      results_df = res_df,
      signature = sig[[pb_name]]
    ))
  }

# Reorganize results by cell line and save
for (idx in 1:length(celllines)) {
  cell_line <- celllines[idx]
  cell_line_clean <- str_replace_all(cell_line, c("-" = "_", " " = "", "/" = "_"))
  
  # Extract signatures for this cell line
  sig_list_9_10th <- list()
  for (split_idx in 1:n_splits) {
    result <- results_9_10th[[idx]][[split_idx]]
    if (!is.null(result)) {
      sig_list_9_10th[[paste0("split_", split_idx, "_excluded")]] <- result$signature
    }
  }
  
  # Save signatures for this cell line
  if (do_save && length(sig_list_9_10th) > 0) {
    saveRDS(sig_list_9_10th,
            file.path(savepath, "mean/90th_perturb", paste0(cell_line_clean, "_90_signatures.rds")))
  }
}

cat("\n=== Analysis Complete ===\n")
