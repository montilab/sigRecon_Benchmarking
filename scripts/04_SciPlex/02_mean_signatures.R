library(tidyverse)
library(DESeq2)
library(Seurat)
library(sigrecon)
library(doParallel)
registerDoParallel(cores=15)

PATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/projects/sig_recon")
savepath <- file.path(Sys.getenv("MLAB"), "projects/brcameta/projects/sig_recon/data/sigs/sciplex")
do_save <- TRUE

# ============================================================================
# 0. Load Pseudobulked Data and Create Drug Splits
# ============================================================================

a549_pb <- readRDS(file.path(PATH, "data/sci_plex/a549_filtered_pb.rds"))
k562_pb <- readRDS(file.path(PATH, "data/sci_plex/k562_filtered_pb.rds"))
mcf7_pb <- readRDS(file.path(PATH, "data/sci_plex/mcf7_filtered_pb.rds"))

# Get shared drugs across all three cell lines
a549_drugs <- unique(a549_pb$product_name)
a549_drugs <- a549_drugs[a549_drugs != "Vehicle"]

k562_drugs <- unique(k562_pb$product_name)
k562_drugs <- k562_drugs[k562_drugs != "Vehicle"]

mcf7_drugs <- unique(mcf7_pb$product_name)
mcf7_drugs <- mcf7_drugs[mcf7_drugs != "Vehicle"]

shared_drugs <- Reduce(intersect, list(a549_drugs, k562_drugs, mcf7_drugs))
cat("Total shared drugs:", length(shared_drugs), "\n")
cat("A549 unique drugs:", length(a549_drugs), "\n")
cat("K562 unique drugs:", length(k562_drugs), "\n")
cat("MCF7 unique drugs:", length(mcf7_drugs), "\n")

n_splits <- 10
# Save drug splits
if (do_save) {
  # Create 10 random splits of drugs
  set.seed(123)  # For reproducibility
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
  
  write.csv(drug_splits, 
            file.path(savepath, "drug_splits.csv"),
            row.names = FALSE)
} else {
  drug_splits <- read.csv(file.path(savepath, "drug_splits.csv"))
}

cat("Drug splits created:\n")
print(colSums(drug_splits[, -1]))  # Show number of drugs per split

# ============================================================================
# Helper Functions
# ============================================================================

# Run DESeq2 analysis on Seurat pseudobulk object
run_deseq2_seurat <- function(seurat_obj, sample_ids_perturbed, sample_ids_control) {
  # Get counts and metadata
  counts <- seurat_obj@assays$RNA$counts
  metadata <- seurat_obj@meta.data
  
  # Subset to selected samples
  all_sample_ids <- c(sample_ids_perturbed, sample_ids_control)
  subset_counts <- counts[, all_sample_ids]
  subset_meta <- metadata[all_sample_ids, ]
  
  # Create condition column
  subset_meta$condition <- ifelse(rownames(subset_meta) %in% sample_ids_perturbed, 
                                  "Perturbed", "Control")
  subset_meta$condition <- as.factor(subset_meta$condition)
  
  # Check if we have replicate information
  if ("rep" %in% colnames(subset_meta)) {
    subset_meta$replicate <- as.factor(subset_meta$replicate)
    design_formula <- ~ replicate + condition
  } else {
    design_formula <- ~ condition
  }
  
  # Create DESeq2 object
  dds <- DESeqDataSetFromMatrix(
    countData = subset_counts,
    colData = subset_meta,
    design = design_formula
  )
  
  # Run DESeq
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, contrast = c("condition", "Perturbed", "Control"))
  
  return(as.data.frame(res))
}

# Get drug sample IDs from Seurat object
get_drug_samples <- function(seurat_obj, drugs) {
  sample_ids <- rownames(seurat_obj@meta.data)[
    seurat_obj@meta.data$product_name %in% drugs
  ]
  return(sample_ids)
}

# Get control sample IDs from Seurat object
get_control_samples <- function(seurat_obj) {
  sample_ids <- rownames(seurat_obj@meta.data)[
    seurat_obj@meta.data$product_name == "Vehicle"
  ]
  return(sample_ids)
}

# ============================================================================
# Analysis 1: Control vs Control (Between Cell Lines)
# ============================================================================

cat("\n=== Analysis 1: Control vs Control ===\n")

# Get control samples from each cell line
a549_ctrl <- get_control_samples(a549_pb)
k562_ctrl <- get_control_samples(k562_pb)
mcf7_ctrl <- get_control_samples(mcf7_pb)

cat("A549 controls:", length(a549_ctrl), "\n")
cat("K562 controls:", length(k562_ctrl), "\n")
cat("MCF7 controls:", length(mcf7_ctrl), "\n")

# Find shared genes
shared_genes <- Reduce(intersect, list(
  rownames(a549_pb@assays$RNA$counts),
  rownames(k562_pb@assays$RNA$counts),
  rownames(mcf7_pb@assays$RNA$counts)
))
cat("Shared genes:", length(shared_genes), "\n")

combined_counts <- cbind(
  a549_pb@assays$RNA$counts[shared_genes, a549_ctrl],
  k562_pb@assays$RNA$counts[shared_genes, k562_ctrl],
  mcf7_pb@assays$RNA$counts[shared_genes, mcf7_ctrl]
)

combined_meta <- rbind(
  a549_pb@meta.data[a549_ctrl, c("product_name", "replicate")],
  k562_pb@meta.data[k562_ctrl, c("product_name", "replicate")],
  mcf7_pb@meta.data[mcf7_ctrl, c("product_name", "replicate")]
)
colnames(combined_counts) <- rownames(combined_meta)
combined_meta$cell_line <- c(rep("a549", length(a549_ctrl)),
                             rep("k562", length(k562_ctrl)),
                             rep("mcf7", length(mcf7_ctrl)))
combined_meta$replicate <- as.factor(as.character(combined_meta$replicate))
combined_meta$is_a549 <- with(combined_meta, if_else(cell_line == "a549", "a549", "other"))
combined_meta$is_k562 <- with(combined_meta, if_else(cell_line == "k562", "k562", "other"))
combined_meta$is_mcf7 <- with(combined_meta, if_else(cell_line == "mcf7", "mcf7", "other"))

# DESeq2 A549
dds_a549 <- DESeqDataSetFromMatrix(
  countData = combined_counts,
  colData = combined_meta,
  design = ~ replicate + is_a549
)
dds_a549 <- DESeq(dds_a549, quiet = TRUE)
res_a549 <- results(dds_a549, contrast = c("is_a549", "a549", "other"))

# Convert to data frames and extract signatures
res_a549 <- as.data.frame(res_a549) %>%
  dplyr::mutate(pb = "a549") %>%
  as_tibble(rownames = "gene") %>%
  arrange(padj, desc(abs(log2FoldChange)))
sig_a549 <- sig_filter_fn(res_a549, perts = "a549", pert_col = "pb", log2fc_col = "log2FoldChange", pval_col = "padj", geneid_col = "gene")

# Save A549
if (do_save) {
  write.csv(res_a549,
            file.path(savepath, "mean/ctrl/a549.csv"),
            row.names = FALSE)
  saveRDS(sig_a549, file.path(savepath, "mean/ctrl/a549.rds"))
}

# DESeq2 k562
dds_k562 <- DESeqDataSetFromMatrix(
  countData = combined_counts,
  colData = combined_meta,
  design = ~ replicate + is_k562
)
dds_k562 <- DESeq(dds_k562, quiet = TRUE)
res_k562 <- results(dds_k562, contrast = c("is_k562", "k562", "other"))

# Convert to data frames and extract signatures
res_k562 <- as.data.frame(res_k562) %>%
  dplyr::mutate(pb = "k562") %>%
  as_tibble(rownames = "gene") %>%
  arrange(padj, desc(abs(log2FoldChange)))
sig_k562 <- sig_filter_fn(res_k562, perts = "k562", pert_col = "pb", log2fc_col = "log2FoldChange", pval_col = "padj", geneid_col = "gene")

# Save k562
if (do_save) {
  write.csv(res_k562,
            file.path(savepath, "mean/ctrl/k562.csv"),
            row.names = FALSE)
  saveRDS(sig_k562, file.path(savepath, "mean/ctrl/k562.rds"))
}

# DESeq2 mcf7
dds_mcf7 <- DESeqDataSetFromMatrix(
  countData = combined_counts,
  colData = combined_meta,
  design = ~ replicate + is_mcf7
)
dds_mcf7 <- DESeq(dds_mcf7, quiet = TRUE)
res_mcf7 <- results(dds_mcf7, contrast = c("is_mcf7", "mcf7", "other"))

# Convert to data frames and extract signatures
res_mcf7 <- as.data.frame(res_mcf7) %>%
  dplyr::mutate(pb = "mcf7") %>%
  as_tibble(rownames = "gene") %>%
  arrange(padj, desc(abs(log2FoldChange)))
sig_mcf7 <- sig_filter_fn(res_mcf7, perts = "mcf7", pert_col = "pb", log2fc_col = "log2FoldChange", pval_col = "padj", geneid_col = "gene")

# Save mcf7
if (do_save) {
  write.csv(res_mcf7,
            file.path(savepath, "mean/ctrl/mcf7.csv"),
            row.names = FALSE)
  saveRDS(sig_mcf7, file.path(savepath, "mean/ctrl/mcf7.rds"))
}


# ============================================================================
# Analysis 2: 1/10th Perturbed vs All Controls (10 splits)
# ============================================================================

cat("\n=== Analysis 2: 1/10th Perturbed vs Controls ===\n")

# A549 - 1/10th splits
sig_a549_1_10th <- foreach(i = 1:n_splits, .combine = 'c') %dopar% {
  cat("Processing A549 split", i, "of", n_splits, "\n")

  # Get drugs in this split
  drugs_in_split <- drug_splits$drug[drug_splits[[paste0("split_", i)]]]

  # Get sample IDs
  drug_sample_ids <- get_drug_samples(a549_pb, drugs_in_split)
  ctrl_sample_ids <- a549_ctrl

  if (length(drug_sample_ids) < 3) {
    cat("  Skipping: not enough drug samples\n")
    next
  }

  cat("  Drugs:", length(drugs_in_split),
      "Drug samples:", length(drug_sample_ids),
      "Control samples:", length(ctrl_sample_ids), "\n")

  # Run DESeq2
  res <- run_deseq2_seurat(
    seurat_obj = a549_pb,
    sample_ids_perturbed = drug_sample_ids,
    sample_ids_control = ctrl_sample_ids
  )
  pb_name <- paste0("a549", "_split_", i)
  res_df <- res %>%
    dplyr::mutate(pb = pb_name) %>%
    as_tibble(rownames = "gene") %>%
    arrange(padj, desc(abs(log2FoldChange)))

  # Extract signatures
  sig <- sig_filter_fn(res_df, perts = pb_name, pert_col = "pb", log2fc_col = "log2FoldChange", pval_col = "padj", alpha = 1, geneid_col = "gene")

  # Save full results
  if (do_save) {
    write.csv(res_df,
              file.path(savepath, paste0("mean/10th_perturb/a549_10per_", i, ".csv")),
              row.names = FALSE)
  }

  if (!is.null(sig)) {
    cat("  Up:", length(sig$up), "Down:", length(sig$dn), "\n")
  }
  sig
}

# K562 - 1/10th splits
sig_k562_1_10th <- foreach(i = 1:n_splits, .combine = 'c') %dopar% {
  cat("Processing K562 split", i, "of", n_splits, "\n")

  drugs_in_split <- drug_splits$drug[drug_splits[[paste0("split_", i)]]]
  drug_sample_ids <- get_drug_samples(k562_pb, drugs_in_split)
  ctrl_sample_ids <- k562_ctrl

  if (length(drug_sample_ids) < 3) {
    cat("  Skipping: not enough drug samples\n")
    next
  }

  cat("  Drugs:", length(drugs_in_split),
      "Drug samples:", length(drug_sample_ids),
      "Control samples:", length(ctrl_sample_ids), "\n")

  res <- run_deseq2_seurat(
    seurat_obj = k562_pb,
    sample_ids_perturbed = drug_sample_ids,
    sample_ids_control = ctrl_sample_ids
  )

  pb_name <- paste0("k562", "_split_", i)
  res_df <- res %>%
    dplyr::mutate(pb = pb_name) %>%
    as_tibble(rownames = "gene") %>%
    arrange(padj, desc(abs(log2FoldChange)))

  # Extract signatures
  sig <- sig_filter_fn(res_df, perts = pb_name, pert_col = "pb", log2fc_col = "log2FoldChange", pval_col = "padj", alpha = 1, geneid_col = "gene")

  if (do_save) {
    write.csv(res_df,
              file.path(savepath, paste0("mean/10th_perturb/k562_10per_", i, ".csv")),
              row.names = FALSE)
  }

  if (!is.null(sig)) {
    cat("  Up:", length(sig$up), "Down:", length(sig$dn), "\n")
  }
  sig
}

# MCF7 - 1/10th splits
sig_mcf7_1_10th <- foreach(i = 1:n_splits, .combine = 'c') %dopar% {

  cat("Processing MCF7 split", i, "of", n_splits, "\n")

  drugs_in_split <- drug_splits$drug[drug_splits[[paste0("split_", i)]]]
  drug_sample_ids <- get_drug_samples(mcf7_pb, drugs_in_split)
  ctrl_sample_ids <- mcf7_ctrl

  if (length(drug_sample_ids) < 3) {
    cat("  Skipping: not enough drug samples\n")
    next
  }

  cat("  Drugs:", length(drugs_in_split),
      "Drug samples:", length(drug_sample_ids),
      "Control samples:", length(ctrl_sample_ids), "\n")

  res <- run_deseq2_seurat(
    seurat_obj = mcf7_pb,
    sample_ids_perturbed = drug_sample_ids,
    sample_ids_control = ctrl_sample_ids
  )

  pb_name <- paste0("mcf7", "_split_", i)
  res_df <- res %>%
    dplyr::mutate(pb = pb_name) %>%
    as_tibble(rownames = "gene") %>%
    arrange(padj, desc(abs(log2FoldChange)))

  # Extract signatures
  sig <- sig_filter_fn(res_df, perts = pb_name, pert_col = "pb", log2fc_col = "log2FoldChange", pval_col = "padj", alpha = 1, geneid_col = "gene")

  if (do_save) {
    write.csv(res_df,
              file.path(savepath, paste0("mean/10th_perturb/mcf7_10per_", i, ".csv")),
              row.names = FALSE)
  }

  if (!is.null(sig)) {
    cat("  Up:", length(sig$up), "Down:", length(sig$dn), "\n")
  }

  sig
}

# Save signatures
if (do_save) {
  saveRDS(sig_a549_1_10th,
          file.path(savepath, "mean/a549_10_signatures.rds"))
  saveRDS(sig_k562_1_10th,
          file.path(savepath, "mean/k562_10_signatures.rds"))
  saveRDS(sig_mcf7_1_10th,
          file.path(savepath, "mean/mcf7_10_signatures.rds"))
}

# ============================================================================
# Analysis 3: 9/10th Perturbed vs All Controls (10 splits)
# ============================================================================

cat("\n=== Analysis 3: 9/10th Perturbed vs Controls ===\n")

# A549 - 9/10th splits
sig_a549_9_10th <- foreach(i = 1:n_splits, .combine = 'c') %dopar% {
  cat("Processing A549 split", i, "of", n_splits, "(excluding split", i, ")\n")
  
  # Get drugs NOT in this split (9/10th)
  drugs_in_9_10th <- drug_splits$drug[!drug_splits[[paste0("split_", i)]]]
  
  drug_sample_ids <- get_drug_samples(a549_pb, drugs_in_9_10th)
  ctrl_sample_ids <- a549_ctrl
  
  if (length(drug_sample_ids) < 3) {
    cat("  Skipping: not enough drug samples\n")
    next
  }
  
  cat("  Drugs:", length(drugs_in_9_10th),
      "Drug samples:", length(drug_sample_ids), 
      "Control samples:", length(ctrl_sample_ids), "\n")
  
  res <- run_deseq2_seurat(
    seurat_obj = a549_pb,
    sample_ids_perturbed = drug_sample_ids,
    sample_ids_control = ctrl_sample_ids
  )
  
  pb_name <- paste0("a549", "_split_", i)
  res_df <- res %>%
    dplyr::mutate(pb = pb_name) %>% 
    as_tibble(rownames = "gene") %>%
    arrange(padj, desc(abs(log2FoldChange)))
  
  # Extract signatures
  sig <- sig_filter_fn(res_df, perts = pb_name, pert_col = "pb", log2fc_col = "log2FoldChange", pval_col = "padj", alpha = 1, geneid_col = "gene")
  
  if (do_save) {
    write.csv(res_df,
              file.path(savepath, paste0("mean/90th_perturb/a549_90_", i, ".csv")),
              row.names = FALSE)
  }
  
  if (!is.null(sig)) {
    cat("  Up:", length(sig$up), "Down:", length(sig$dn), "\n")
  }
  sig
}

# K562 - 9/10th splits
sig_k562_9_10th <- foreach(i = 1:n_splits, .combine = 'c') %dopar% {
  cat("Processing K562 split", i, "of", n_splits, "(excluding split", i, ")\n")
  
  drugs_in_9_10th <- drug_splits$drug[!drug_splits[[paste0("split_", i)]]]
  drug_sample_ids <- get_drug_samples(k562_pb, drugs_in_9_10th)
  ctrl_sample_ids <- k562_ctrl
  
  if (length(drug_sample_ids) < 3) {
    cat("  Skipping: not enough drug samples\n")
    next
  }
  
  cat("  Drugs:", length(drugs_in_9_10th),
      "Drug samples:", length(drug_sample_ids), 
      "Control samples:", length(ctrl_sample_ids), "\n")
  
  res <- run_deseq2_seurat(
    seurat_obj = k562_pb,
    sample_ids_perturbed = drug_sample_ids,
    sample_ids_control = ctrl_sample_ids
  )
  
  pb_name <- paste0("k562", "_split_", i)
  res_df <- res %>%
    dplyr::mutate(pb = pb_name) %>% 
    as_tibble(rownames = "gene") %>%
    arrange(padj, desc(abs(log2FoldChange)))
  
  # Extract signatures
  sig <- sig_filter_fn(res_df, perts = pb_name, pert_col = "pb", log2fc_col = "log2FoldChange", pval_col = "padj", alpha = 1, geneid_col = "gene")
  
  if (do_save) {
    write.csv(res_df,
              file.path(savepath, paste0("mean/90th_perturb/k562_90_", i, ".csv")),
              row.names = FALSE)
  }
  
  if (!is.null(sig)) {
    cat("  Up:", length(sig$up), "Down:", length(sig$dn), "\n")
  }
  sig
}

# MCF7 - 9/10th splits
sig_mcf7_9_10th <- foreach(i = 1:n_splits, .combine = 'c') %dopar% {
  cat("Processing MCF7 split", i, "of", n_splits, "(exclßding split", i, ")\n")
  
  drugs_in_9_10th <- drug_splits$drug[!drug_splits[[paste0("split_", i)]]]
  drug_sample_ids <- get_drug_samples(mcf7_pb, drugs_in_9_10th)
  ctrl_sample_ids <- mcf7_ctrl
  
  if (length(drug_sample_ids) < 3) {
    cat("  Skipping: not enough drug samples\n")
    next
  }
  
  cat("  Drugs:", length(drugs_in_9_10th),
      "Drug samples:", length(drug_sample_ids), 
      "Control samples:", length(ctrl_sample_ids), "\n")
  
  res <- run_deseq2_seurat(
    seurat_obj = mcf7_pb,
    sample_ids_perturbed = drug_sample_ids,
    sample_ids_control = ctrl_sample_ids
  )
  
  pb_name <- paste0("mcf7", "_split_", i)
  res_df <- res %>%
    dplyr::mutate(pb = pb_name) %>% 
    as_tibble(rownames = "gene") %>%
    arrange(padj, desc(abs(log2FoldChange)))
  
  # Extract signatures
  sig <- sig_filter_fn(res_df, perts = pb_name, pert_col = "pb", log2fc_col = "log2FoldChange", pval_col = "padj", alpha = 1, geneid_col = "gene")
  
  if (do_save) {
    write.csv(res_df,
              file.path(savepath, paste0("mean/90th_perturb/mcf7_90_", i, ".csv")),
              row.names = FALSE)
  }
  
  if (!is.null(sig)) {
    cat("  Up:", length(sig$up), "Down:", length(sig$dn), "\n")
  }
  sig
}

# Save signatures
if (do_save) {
  saveRDS(sig_a549_9_10th,
          file.path(savepath, "mean/a549_90_signatures.rds"))
  saveRDS(sig_k562_9_10th,
          file.path(savepath, "mean/k562_90_signatures.rds"))
  saveRDS(sig_mcf7_9_10th,
          file.path(savepath, "mean/mcf7_90_signatures.rds"))
}