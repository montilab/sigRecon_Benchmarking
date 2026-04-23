library(DESeq2)
library(tidyverse)
library(doParallel)
library(sigrecon)
registerDoParallel(cores=10)

PATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/projects/sig_recon/data/perturb_seq")
savepath <- file.path(Sys.getenv("MLAB"), "projects/brcameta/projects/sig_recon/data/sigs/perturb-seq")
do_save <- TRUE

# ============================================================================
# 0. Load Data and Create PB Splits
# ============================================================================

k562_counts <- read.csv(file.path(PATH, "k562_processed_pb.csv"), row.names = 1)
colnames(k562_counts) <- str_replace_all(colnames(k562_counts),
                                         pattern = "\\.",
                                         replacement = "-")
k562_meta <- read.csv(file.path(PATH, "k562_processed_pb_metadata.csv"), row.names = 1)
stopifnot(all(colnames(k562_counts) %in% rownames(k562_meta)))
k562_meta$gem_group <- as.factor(k562_meta$gem_group)
k562_meta$is_control <- ifelse(k562_meta$gene == "non-targeting", "NTC", "Perturbed")

rpe1_counts <- read.csv(file.path(PATH, "rpe1_processed_pb.csv"), row.names = 1)
colnames(rpe1_counts) <- str_replace_all(colnames(rpe1_counts),
                                         pattern = "\\.",
                                         replacement = "-")
rpe1_meta <- read.csv(file.path(PATH, "rpe1_processed_pb_metadata.csv"), row.names = 1)
rownames(rpe1_meta) <- str_replace_all(rownames(rpe1_meta),
                                       pattern = "AC118549\\.",
                                       replacement = "AC118549-")
stopifnot(all(colnames(rpe1_counts) %in% rownames(rpe1_meta)))
rpe1_meta$gem_group <- as.factor(rpe1_meta$gem_group)
rpe1_meta$is_control <- ifelse(rpe1_meta$gene == "non-targeting", "NTC", "Perturbed")

# Get shared perturbations between k562 and rpe1
k562_pbs <- unique(k562_meta$gene[k562_meta$gene != "non-targeting"])
rpe1_pbs <- unique(rpe1_meta$gene[rpe1_meta$gene != "non-targeting"])
shared_pbs <- intersect(k562_pbs, rpe1_pbs)

cat("Total shared perturbations:", length(shared_pbs), "\n")
cat("K562 unique pbs:", length(k562_pbs), "\n")
cat("RPE1 unique pbs:", length(rpe1_pbs), "\n")

# Create 10 random splits of perturbations
set.seed(123)  # For reproducibility
n_splits <- 10
pbs_shuffled <- sample(shared_pbs)
split_size <- ceiling(length(shared_pbs) / n_splits)

# Create splits dataframe
pb_splits <- data.frame(
  pb = shared_pbs,
  stringsAsFactors = FALSE
)

# Initialize split columns
for (i in 1:n_splits) {
  split_col <- paste0("split_", i)
  pb_splits[[split_col]] <- FALSE
}

# Assign pbs to their primary split
for (i in 1:n_splits) {
  start_idx <- (i - 1) * split_size + 1
  end_idx <- min(i * split_size, length(pbs_shuffled))
  pbs_in_split <- pbs_shuffled[start_idx:end_idx]
  
  pb_splits[[paste0("split_", i)]][pb_splits$pb %in% pbs_in_split] <- TRUE
}

# Save pb splits
if (do_save) {
  write.csv(pb_splits, 
            file.path(savepath, "pb_splits.csv"),
            row.names = FALSE)
}

cat("Perturbation splits created:\n")
print(colSums(pb_splits[, -1]))  # Show number of pbs per split

# ============================================================================
# Helper Functions
# ============================================================================

# Run DESeq2 analysis
run_deseq2 <- function(counts, metadata, sample_ids_perturbed, sample_ids_control) {
  # Subset to selected samples
  all_sample_ids <- c(sample_ids_perturbed, sample_ids_control)
  subset_counts <- counts[, all_sample_ids]
  subset_meta <- metadata[all_sample_ids, ]
  
  # PRE-FILTER: Remove genes with near-zero expression
  keep <- rowSums(subset_counts >= 10) >= 5  # At least 5 samples with 10+ counts
  subset_counts <- subset_counts[keep, ]
  
  cat("  Filtered to", sum(keep), "genes\n")
  
  # Create group column
  subset_meta$condition <- ifelse(rownames(subset_meta) %in% sample_ids_perturbed, 
                                  "Perturbed", "Control")
  subset_meta$condition <- as.factor(subset_meta$condition)
  
  # Create DESeq2 object
  dds <- DESeqDataSetFromMatrix(
    countData = subset_counts,
    colData = subset_meta,
    design = ~ gem_group + condition
  )
  
  # Run DESeq
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, contrast = c("condition", "Perturbed", "Control"))
  
  return(as.data.frame(res))
}


# ============================================================================
# Analysis 1: Control RPE1 vs Control K562
# ============================================================================

# cat("\n=== Analysis 1: Control vs Control ===\n")
# 
# Get control samples
rpe1_ctrl <- rownames(rpe1_meta)[rpe1_meta$is_control == "NTC"]
k562_ctrl <- rownames(k562_meta)[k562_meta$is_control == "NTC"]
# 
# cat("RPE1 controls:", length(rpe1_ctrl), "\n")
# cat("K562 controls:", length(k562_ctrl), "\n")
# 
# # Find shared genes between datasets
# shared_genes <- intersect(rownames(k562_counts), rownames(rpe1_counts))
# cat("Shared genes:", length(shared_genes), "\n")
# 
# # Combine count matrices
# combined_counts <- cbind(
#   k562_counts[shared_genes, k562_ctrl],
#   rpe1_counts[shared_genes, rpe1_ctrl]
# )
# 
# # Create combined metadata
# combined_meta <- rbind(
#   k562_meta[k562_ctrl, c("gem_group", "gene")],
#   rpe1_meta[rpe1_ctrl, c("gem_group", "gene")]
# )
# # rbind changes rownames to make them unique
# colnames(combined_counts) <- rownames(combined_meta)
# combined_meta$cell_line <- c(rep("k562", length(k562_ctrl)),
#                              rep("rpe1", length(rpe1_ctrl)))
# combined_meta$gem_group <- as.factor(as.character(combined_meta$gem_group))
# 
# # Run DESeq2: RPE1 vs K562
# cat("Running DESeq2 for RPE1 vs K562 controls...\n")
# k562_sample_ids <- colnames(combined_counts)[combined_meta$cell_line == "k562"]
# rpe1_sample_ids <- colnames(combined_counts)[combined_meta$cell_line == "rpe1"]
# 
# dds_ctrl <- DESeqDataSetFromMatrix(
#   countData = combined_counts,
#   colData = combined_meta,
#   design = ~ gem_group + cell_line
# )
# 
# dds_ctrl <- DESeq(dds_ctrl, quiet = TRUE)
# res_rpe1_vs_k562 <- results(dds_ctrl, contrast = c("cell_line", "rpe1", "k562"))
# res_k562_vs_rpe1 <- results(dds_ctrl, contrast = c("cell_line", "k562", "rpe1"))
# 
# # Convert to data frames
# res_rpe1_vs_k562_df <- as.data.frame(res_rpe1_vs_k562) %>%
#   dplyr::mutate(pb = "rpe1") %>%
#   as_tibble(rownames = "gene") %>%
#   arrange(padj, desc(abs(log2FoldChange)))
# 
# res_k562_vs_rpe1_df <- as.data.frame(res_k562_vs_rpe1) %>%
#   dplyr::mutate(pb = "k562") %>%
#   as_tibble(rownames = "gene") %>%
#   arrange(padj, desc(abs(log2FoldChange)))
# 
# # Extract signatures
# sig_rpe1_vs_k562 <- sig_filter_fn(res_rpe1_vs_k562_df, perts = "rpe1", pert_col = "pb", log2fc_col = "log2FoldChange", pval_col = "padj", geneid_col = "gene")
# sig_k562_vs_rpe1 <- sig_filter_fn(res_k562_vs_rpe1_df, perts = "k562", pert_col = "pb", log2fc_col = "log2FoldChange", pval_col = "padj", geneid_col = "gene")
# 
# # Save
# if (do_save) {
#   write.csv(res_rpe1_vs_k562_df,
#             file.path(savepath, "mean/ctrl/ctrl_rpe1_vs_k562.csv"),
#             row.names = FALSE)
#   write.csv(res_k562_vs_rpe1_df,
#             file.path(savepath, "mean/ctrl/ctrl_k562_vs_rpe1.csv"),
#             row.names = FALSE)
#   
#   saveRDS(sig_rpe1_vs_k562, file.path(savepath, "mean/ctrl/rpe1_ctrl.rds"))
#   saveRDS(sig_k562_vs_rpe1, file.path(savepath, "mean/ctrl/k562_ctrl.rds"))
# }
# 
# cat("RPE1 vs K562 - Up:", length(sig_rpe1_vs_k562$rpe1$up), 
#     "Down:", length(sig_rpe1_vs_k562$dn), "\n")
# cat("K562 vs RPE1 - Up:", length(sig_k562_vs_rpe1$k562$up), 
#     "Down:", length(sig_k562_vs_rpe1$dn), "\n")

# ============================================================================
# Analysis 2: 1/10th Perturbed vs All Controls
# ============================================================================

cat("\n=== Analysis 2: 1/10th Perturbed vs Controls ===\n")

# RPE1 - 1/10th splits
sig_rpe1_1_10th <- foreach(i = 1:n_splits, .combine = 'c') %dopar% {
  cat("Processing RPE1 split", i, "of", n_splits, "\n")

  # Get pbs in this split
  pbs_in_split <- pb_splits$pb[pb_splits[[paste0("split_", i)]]]

  # Sample 1/10th of replicates per gene
  set.seed(42)  # For reproducibility
  
  perturbed_samples_list <- lapply(pbs_in_split, function(gene) {
    # Get all samples for this gene
    gene_samples <- rownames(rpe1_meta)[
      rpe1_meta$gene == gene & rpe1_meta$is_control == "Perturbed"
    ]
    
    # Sample 1/10th (at least 1 if possible)
    n_sample <- max(1, floor(length(gene_samples) / 10))
    
    if (length(gene_samples) <= n_sample) {
      return(gene_samples)  # Return all if too few
    } else {
      return(sample(gene_samples, n_sample))
    }
  })
  
  # Combine all sampled replicates
  perturbed_samples <- unlist(perturbed_samples_list)
  
  control_samples <- rpe1_ctrl

  if (length(perturbed_samples) < 3) {
    cat("  Skipping: not enough perturbed samples\n")
    next
  }

  cat("  PBs:", length(pbs_in_split),
      "Perturbed samples:", length(perturbed_samples),
      "Control samples:", length(control_samples), "\n")

  # Run DESeq2
  res <- run_deseq2(
    counts = rpe1_counts,
    metadata = rpe1_meta,
    sample_ids_perturbed = perturbed_samples,
    sample_ids_control = control_samples
  )

  pb_name <- paste0("rpe1", "_split_", i)
  res_df <- res %>%
    as.data.frame() %>%
    dplyr::mutate(pb = pb_name) %>%
    as_tibble(rownames = "gene") %>%
    arrange(padj, desc(abs(log2FoldChange)))

  sig <- sig_filter_fn(res_df, perts = pb_name, pert_col = "pb", log2fc_col = "log2FoldChange", pval_col = "padj", geneid_col = "gene")

  # Save full results
  if (do_save) {
    write.csv(res_df,
              file.path(savepath, paste0("mean/10th_perturb/rpe1_10per_", i, ".csv")),
              row.names = FALSE)
  }

  if (!is.null(sig)) {
    cat("  Up:", length(sig$up), "\n")
  }
  sig
}

# K562 - 1/10th splits
sig_k562_1_10th <- foreach(i = 1:n_splits, .combine = 'c') %dopar% {
  cat("Processing K562 split", i, "of", n_splits, "\n")

  pbs_in_split <- pb_splits$pb[pb_splits[[paste0("split_", i)]]]
  
  # Sample 1/10th of replicates per gene
  set.seed(42)  # For reproducibility
  
  perturbed_samples_list <- lapply(pbs_in_split, function(gene) {
    # Get all samples for this gene
    gene_samples <- rownames(k562_meta)[
      k562_meta$gene == gene & k562_meta$is_control == "Perturbed"
    ]
    
    # Sample 1/10th (at least 1 if possible)
    n_sample <- max(1, floor(length(gene_samples) / 10))
    
    if (length(gene_samples) <= n_sample) {
      return(gene_samples)  # Return all if too few
    } else {
      return(sample(gene_samples, n_sample))
    }
  })
  
  # Combine all sampled replicates
  perturbed_samples <- unlist(perturbed_samples_list)
  control_samples <- k562_ctrl

  if (length(perturbed_samples) < 3) {
    cat("  Skipping: not enough perturbed samples\n")
    next
  }

  cat("  PBs:", length(pbs_in_split),
      "Perturbed samples:", length(perturbed_samples),
      "Control samples:", length(control_samples), "\n")

  res <- run_deseq2(
    counts = k562_counts,
    metadata = k562_meta,
    sample_ids_perturbed = perturbed_samples,
    sample_ids_control = control_samples
  )

  pb_name <- paste0("k562", "_split_", i)
  res_df <- res %>%
    as.data.frame() %>%
    dplyr::mutate(pb = pb_name) %>%
    as_tibble(rownames = "gene") %>%
    arrange(padj, desc(abs(log2FoldChange)))

  sig <- sig_filter_fn(res_df, perts = pb_name, pert_col = "pb", log2fc_col = "log2FoldChange", pval_col = "padj", geneid_col = "gene")

  if (do_save) {
    write.csv(res_df,
              file.path(savepath, paste0("mean/10th_perturb/k562_10per_", i, ".csv")),
              row.names = FALSE)
  }

  if (!is.null(sig)) {
    cat("  Up:", length(sig$up), "\n")
  }
  sig
}

# Save signatures
if (do_save) {
  saveRDS(sig_rpe1_1_10th,
          file.path(savepath, "mean/10th_perturb/rpe1_10_signatures.rds"))
  saveRDS(sig_k562_1_10th,
          file.path(savepath, "mean/10th_perturb/k562_10_signatures.rds"))
}

# ============================================================================
# Analysis 3: 9/10th Perturbed vs All Controls
# ============================================================================

cat("\n=== Analysis 3: 9/10th Perturbed vs Controls ===\n")

# RPE1 - 9/10th splits
sig_rpe1_9_10th <- foreach(i = 1:n_splits, .combine = 'c') %dopar% {
  cat("Processing RPE1 split", i, "of", n_splits, "(excluding split", i, ")\n")
  
  # Get pbs NOT in this split (9/10th)
  pbs_in_9_10th <- pb_splits$pb[!pb_splits[[paste0("split_", i)]]]
  
  # Sample 1/10th of replicates per gene
  set.seed(42)  # For reproducibility
  
  perturbed_samples_list <- lapply(pbs_in_9_10th, function(gene) {
    # Get all samples for this gene
    gene_samples <- rownames(rpe1_meta)[
      rpe1_meta$gene == gene & rpe1_meta$is_control == "Perturbed"
    ]
    
    # Sample 1/10th (at least 1 if possible)
    n_sample <- max(1, floor(length(gene_samples) / 10))
    
    if (length(gene_samples) <= n_sample) {
      return(gene_samples)  # Return all if too few
    } else {
      return(sample(gene_samples, n_sample))
    }
  })
  
  # Combine all sampled replicates
  perturbed_samples <- unlist(perturbed_samples_list)
  control_samples <- rpe1_ctrl
  
  if (length(perturbed_samples) < 3) {
    cat("  Skipping: not enough perturbed samples\n")
    next
  }
  
  cat("  PBs:", length(pbs_in_9_10th),
      "Perturbed samples:", length(perturbed_samples), 
      "Control samples:", length(control_samples), "\n")
  
  res <- run_deseq2(
    counts = rpe1_counts,
    metadata = rpe1_meta,
    sample_ids_perturbed = perturbed_samples,
    sample_ids_control = control_samples
  )
  
  pb_name <- paste0("rpe1", "_split_", i)
  res_df <- res %>%
    as.data.frame() %>%
    dplyr::mutate(pb = pb_name) %>% 
    as_tibble(rownames = "gene") %>%
    arrange(padj, desc(abs(log2FoldChange)))
  
  sig <- sig_filter_fn(res_df, perts = pb_name, pert_col = "pb", log2fc_col = "log2FoldChange", pval_col = "padj", geneid_col = "gene")
  
  if (do_save) {
    write.csv(res_df,
              file.path(savepath, paste0("mean/90th_perturb/rpe1_90_", i, ".csv")),
              row.names = FALSE)
  }
  
  if (!is.null(sig)) {
    cat("  Up:", length(sig[[pb_name]][["up"]]), "\n")
  }
  sig
}

saveRDS(sig_rpe1_9_10th,
        file.path(savepath, "mean/90th_perturb/rpe1_90_signatures.rds"))

# K562 - 9/10th splits
sig_k562_9_10th <- foreach(i = 1:n_splits, .combine = 'c') %dopar% {
  cat("Processing K562 split", i, "of", n_splits, "(excluding split", i, ")\n")
  
  pbs_in_9_10th <- pb_splits$pb[!pb_splits[[paste0("split_", i)]]]
  # Sample 1/10th of replicates per gene
  set.seed(42)  # For reproducibility
  
  perturbed_samples_list <- lapply(pbs_in_9_10th, function(gene) {
    # Get all samples for this gene
    gene_samples <- rownames(k562_meta)[
      k562_meta$gene == gene & k562_meta$is_control == "Perturbed"
    ]
    
    # Sample 1/10th (at least 1 if possible)
    n_sample <- max(1, floor(length(gene_samples) / 10))
    
    if (length(gene_samples) <= n_sample) {
      return(gene_samples)  # Return all if too few
    } else {
      return(sample(gene_samples, n_sample))
    }
  })
  
  # Combine all sampled replicates
  perturbed_samples <- unlist(perturbed_samples_list)
  control_samples <- k562_ctrl
  
  if (length(perturbed_samples) < 3) {
    cat("  Skipping: not enough perturbed samples\n")
    next
  }
  
  cat("  PBs:", length(pbs_in_9_10th),
      "Perturbed samples:", length(perturbed_samples), 
      "Control samples:", length(control_samples), "\n")
  
  res <- run_deseq2(
    counts = k562_counts,
    metadata = k562_meta,
    sample_ids_perturbed = perturbed_samples,
    sample_ids_control = control_samples
  )
  
  pb_name <- paste0("k562", "_split_", i)
  res_df <- res %>%
    as.data.frame() %>%
    dplyr::mutate(pb = pb_name) %>% 
    as_tibble(rownames = "gene") %>%
    arrange(padj, desc(abs(log2FoldChange)))
  
  sig <- sig_filter_fn(res_df, perts = pb_name, pert_col = "pb", log2fc_col = "log2FoldChange", pval_col = "padj", geneid_col = "gene")
  
  if (do_save) {
    write.csv(res_df,
              file.path(savepath, paste0("mean/90th_perturb/k562_90_", i, ".csv")),
              row.names = FALSE)
  }
  
  if (!is.null(sig)) {
    cat("  Up:", length(sig[[pb_name]][["up"]]), "\n")
  }
  sig
}

# Save signatures
saveRDS(sig_k562_9_10th,
        file.path(savepath, "mean/90th_perturb/k562_90_signatures.rds"))

cat("\n=== Analysis Complete ===\n")