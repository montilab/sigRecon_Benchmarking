library(tidyverse)
library(limma)
library(Biobase)
library(sigrecon)

PATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/projects/sig_recon")
DATA_PATH <- file.path(Sys.getenv("AGED"), "CBMrepositoryData/perturbational_data/drugmatrix/")
do_save <- TRUE

# ============================================================================
# 0. Load Data and Create Drug Splits
# ============================================================================

# Load signature data to get drug lists
data("drugmatrix.liver")
data("drugmatrix.kidney")

# Get shared drugs between liver and kidney
shared_drugs <- intersect(names(drugmatrix.liver), names(drugmatrix.kidney))
cat("Total shared drugs:", length(shared_drugs), "\n")

# Create 10 random splits of drugs
set.seed(123)  # For reproducibility
n_splits <- 10
drugs_shuffled <- sample(shared_drugs)
split_size <- ceiling(length(shared_drugs) / n_splits)

# Create splits dataframe
drug_splits <- data.frame(
  drug = shared_drugs,
  stringsAsFactors = FALSE
)

# Assign each drug to splits (a drug can appear in multiple splits for the 9/10 analysis)
for (i in 1:n_splits) {
  split_col <- paste0("split_", i)
  drug_splits[[split_col]] <- FALSE
}

# Assign drugs to their primary split (for 1/10 analysis)
for (i in 1:n_splits) {
  start_idx <- (i - 1) * split_size + 1
  end_idx <- min(i * split_size, length(drugs_shuffled))
  drugs_in_split <- drugs_shuffled[start_idx:end_idx]
  
  drug_splits[[paste0("split_", i)]][drug_splits$drug %in% drugs_in_split] <- TRUE
}

# Save drug splits
if (do_save) {
  write.csv(drug_splits, 
            file.path(PATH, "data/sigs/drugmatrix/drug_splits.csv"),
            row.names = FALSE)
}

cat("Drug splits created:\n")
print(colSums(drug_splits[, -1]))  # Show number of drugs per split

# ============================================================================
# Load Expression Sets
# ============================================================================

liver_eset <- readRDS(file.path(DATA_PATH, "liver.rds"))
kidney_eset <- readRDS(file.path(DATA_PATH, "kidney.rds"))

# ============================================================================
# Helper Functions
# ============================================================================

# Run Limma analysis on ExpressionSet
run_limma <- function(eset, group1_ids, group2_ids, group1_name = "group1", group2_name = "group2") {
  # Check for sufficient samples
  if (length(group1_ids) < 2 || length(group2_ids) < 2) {
    warning("Need at least 2 samples per group")
    return(NULL)
  }
  
  # Subset ExpressionSet to selected samples
  all_sample_ids <- c(group1_ids, group2_ids)
  
  # Check that all samples exist in the eset
  available_samples <- intersect(all_sample_ids, sampleNames(eset))
  if (length(available_samples) != length(all_sample_ids)) {
    warning("Some sample IDs not found in ExpressionSet")
    missing <- setdiff(all_sample_ids, available_samples)
    cat("  Missing samples:", length(missing), "\n")
  }
  
  eset_sub <- eset[, available_samples]
  
  # Create group labels
  group_labels <- c(
    rep(group1_name, sum(available_samples %in% group1_ids)),
    rep(group2_name, sum(available_samples %in% group2_ids))
  )
  
  # Ensure sample names are unique and reflect groups
  sampleNames(eset_sub) <- paste0(group_labels, "_", seq_along(group_labels))
  
  # Create design matrix
  group <- as.factor(group_labels)
  design <- model.matrix(~ 0 + group)
  colnames(design) <- levels(group)
  
  # Fit linear model
  fit <- limma::lmFit(eset_sub, design)
  
  # Build contrast matrix (group1 vs group2)
  cont.matrix <- limma::makeContrasts(
    contrasts = paste0(group1_name, " - ", group2_name),
    levels = design
  )
  
  # Apply contrasts and empirical Bayes
  fit2 <- limma::contrasts.fit(fit, cont.matrix)
  fit2 <- limma::eBayes(fit2)
  
  # Extract results
  de_table <- limma::topTable(fit2, coef = 1, adjust = "fdr", number = Inf)
  
  return(de_table)
}

# Extract Signatures from a DE table
extract_signatures <- function(de_table, gene_symbol_col_in_table, fc_thresh = 0, p_val_thresh = 0.05, max_sig = 100) {
  sigs_result <- list(up = character(0), dn = character(0), up_full = character(0), dn_full = character(0))
  
  de_table$logFC_adjpval <- de_table$logFC * (-log10(de_table$adj.P.Val))
  
  if (nrow(de_table) == 0 || !gene_symbol_col_in_table %in% colnames(de_table)) {
    return(sigs_result)
  }
  
  de_table <- de_table %>%
    filter(!is.na(!!sym(gene_symbol_col_in_table)) & !!sym(gene_symbol_col_in_table) != "")
  
  if (nrow(de_table) == 0) {
    return(sigs_result)
  }
  
  de_table$logFC_adjpval <- de_table$logFC * (-log10(de_table$adj.P.Val))
  
  # Up-regulated signatures
  up_sig <- de_table %>%
    dplyr::filter(logFC > fc_thresh, adj.P.Val < p_val_thresh) %>%
    dplyr::arrange(desc(logFC_adjpval)) %>%
    dplyr::pull(!!sym(gene_symbol_col_in_table))
  
  up_full_sig <- de_table %>%
    dplyr::arrange(desc(logFC_adjpval)) %>%
    dplyr::pull(!!sym(gene_symbol_col_in_table))
  
  # Down-regulated signatures
  dn_sig <- de_table %>%
    dplyr::filter(logFC < -fc_thresh, adj.P.Val < p_val_thresh) %>%
    dplyr::arrange(logFC_adjpval) %>%
    dplyr::pull(!!sym(gene_symbol_col_in_table))
  
  dn_full_sig <- de_table %>%
    dplyr::arrange(logFC_adjpval) %>%
    dplyr::pull(!!sym(gene_symbol_col_in_table))
  
  # Removing duplicates
  up_sig <- up_sig[!duplicated(up_sig)]
  up_full_sig <- up_full_sig[!duplicated(up_full_sig)]
  dn_sig <- dn_sig[!duplicated(dn_sig)]
  dn_full_sig <- dn_full_sig[!duplicated(dn_full_sig)]
  
  # Apply max_signature limit
  up_sig_f <- up_sig[1:min(length(up_sig), max_sig)]
  dn_sig_f <- dn_sig[1:min(length(dn_sig), max_sig)]
  
  # Remove any potential empty strings that might have slipped through
  up_sig_f <- up_sig_f[up_sig_f != ""]
  dn_sig_f <- dn_sig_f[dn_sig_f != ""]
  up_full_sig <- up_full_sig[up_full_sig != ""]
  dn_full_sig <- dn_full_sig[dn_full_sig != ""]
  
  # Return results if significant genes are found for shortlisted signatures
  if (length(up_sig_f) >= 5) {
    sigs_result <- list(up = up_sig_f, dn = dn_sig_f,
                        up_full = up_full_sig, dn_full = dn_full_sig)
  } else {
    message(paste("Warning: Not enough significant genes for shortlist for this table (up:", length(up_sig_f)))
    sigs_result <- NULL
  }
  return(sigs_result)
}

# Helper function to get drug sample IDs
get_drug_samples <- function(eset, drugs) {
  sample_meta <- pData(eset)
  sample_ids <- rownames(sample_meta)[
    sample_meta$`compound:ch1` %in% drugs & 
      sample_meta$`dose:ch1` != "0 mg/kg"
  ]
  return(sample_ids)
}

# ============================================================================
# Analysis 1: Control Liver vs Control Kidney
# ============================================================================

cat("\n=== Analysis 1: Control vs Control ===\n")

# Get all control samples
liver_ctrl <- pData(liver_eset) %>%
  filter(`dose:ch1` == "0 mg/kg") %>%
  rownames()

kidney_ctrl <- pData(kidney_eset) %>%
  filter(`dose:ch1` == "0 mg/kg") %>%
  rownames()

cat("Liver controls:", length(liver_ctrl), "\n")
cat("Kidney controls:", length(kidney_ctrl), "\n")

stopifnot(all.equal(featureNames(liver_eset), featureNames(kidney_eset)))

# Combine expression data into a single ExpressionSet
combined_expr <- cbind(exprs(liver_eset)[, liver_ctrl],
                       exprs(kidney_eset)[, kidney_ctrl])

# Create phenotype data
combined_pdata <- rbind(
  pData(liver_eset)[liver_ctrl, ],
  pData(kidney_eset)[kidney_ctrl, ]
)
combined_pdata$tissue_type <- c(rep("liver", length(liver_ctrl)),
                                rep("kidney", length(kidney_ctrl)))

# Create combined ExpressionSet
combined_eset <- ExpressionSet(
  assayData = combined_expr,
  phenoData = AnnotatedDataFrame(combined_pdata),
  featureData = featureData(liver_eset)
)

# Get sample IDs for each tissue
liver_sample_ids <- sampleNames(combined_eset)[combined_eset$tissue_type == "liver"]
kidney_sample_ids <- sampleNames(combined_eset)[combined_eset$tissue_type == "kidney"]

# Run limma: Liver vs Kidney
cat("Running limma for liver vs kidney controls...\n")
res_liver_vs_kidney <- run_limma(
  eset = combined_eset,
  group1_ids = liver_sample_ids,
  group2_ids = kidney_sample_ids,
  group1_name = "liver",
  group2_name = "kidney"
)

# Run limma: Kidney vs Liver (reverse comparison)
res_kidney_vs_liver <- run_limma(
  eset = combined_eset,
  group1_ids = kidney_sample_ids,
  group2_ids = liver_sample_ids,
  group1_name = "kidney",
  group2_name = "liver"
)

# Convert to tibble
res_liver_vs_kidney_df <- res_liver_vs_kidney %>%
  as_tibble(rownames = "gene_id") %>%
  arrange(adj.P.Val, desc(abs(logFC)))

res_kidney_vs_liver_df <- res_kidney_vs_liver %>%
  as_tibble(rownames = "gene_id") %>%
  arrange(adj.P.Val, desc(abs(logFC)))

# Extract signatures
sig_liver_vs_kidney <- extract_signatures(res_liver_vs_kidney, gene_symbol_col_in_table = "Gene.Symbol")
sig_kidney_vs_liver <- extract_signatures(res_kidney_vs_liver, gene_symbol_col_in_table = "Gene.Symbol")

# Save
if (do_save) {
  write.csv(res_liver_vs_kidney_df,
            file.path(PATH, "data/sigs/drugmatrix/mean/ctrl_liver_vs_kidney.csv"),
            row.names = FALSE)
  write.csv(res_kidney_vs_liver_df,
            file.path(PATH, "data/sigs/drugmatrix/mean/ctrl_kidney_vs_liver.csv"),
            row.names = FALSE)
  
  saveRDS(sig_liver_vs_kidney, file.path(PATH, "data/sigs/drugmatrix/mean/ctrl/liver_ctrl.rds"))
  saveRDS(sig_kidney_vs_liver, file.path(PATH, "data/sigs/drugmatrix/mean/ctrl/kidney_ctrl.rds"))
}

cat("Liver vs Kidney - Up:", length(sig_liver_vs_kidney$up), 
    "Down:", length(sig_liver_vs_kidney$dn), "\n")
cat("Kidney vs Liver - Up:", length(sig_kidney_vs_liver$up), 
    "Down:", length(sig_kidney_vs_liver$dn), "\n")

# ============================================================================
# Analysis 2: 1/10th Perturbed vs All Controls (10 splits)
# ============================================================================

cat("\n=== Analysis 2: 1/10th Perturbed vs Controls ===\n")

# Liver - 1/10th splits
sig_liver_1_10th <- list()

for (i in 1:n_splits) {
  cat("Processing liver split", i, "of", n_splits, "\n")
  
  # Get drugs in this split
  drugs_in_split <- drug_splits$drug[drug_splits[[paste0("split_", i)]]]
  
  # Get sample IDs
  drug_sample_ids <- get_drug_samples(liver_eset, drugs_in_split)
  ctrl_sample_ids <- liver_ctrl
  
  if (length(drug_sample_ids) < 3) {
    cat("  Skipping: not enough drug samples\n")
    next
  }
  
  cat("  Drugs:", length(drugs_in_split), 
      "Drug samples:", length(drug_sample_ids), 
      "Control samples:", length(ctrl_sample_ids), "\n")
  
  # Run limma
  res <- run_limma(
    eset = liver_eset,
    group1_ids = drug_sample_ids,
    group2_ids = ctrl_sample_ids,
    group1_name = "drug",
    group2_name = "control"
  )
  
  if (is.null(res)) {
    cat("  Skipping: limma failed\n")
    next
  }
  
  res_df <- res %>%
    as_tibble(rownames = "gene_id") %>%
    arrange(adj.P.Val, desc(abs(logFC)))
  
  # Extract signatures
  sig <- extract_signatures(res, gene_symbol_col_in_table = "Gene.Symbol")
  sig_liver_1_10th[[paste0("split_", i)]] <- sig
  
  # Save full results
  if (do_save) {
    write.csv(res_df,
              file.path(PATH, paste0("data/sigs/drugmatrix/mean/10th_perturb/liver_10per_", i, ".csv")),
              row.names = FALSE)
  }
  
  cat("  Up:", length(sig$up))
}

# Kidney - 1/10th splits
sig_kidney_1_10th <- list()

for (i in 1:n_splits) {
  cat("Processing kidney split", i, "of", n_splits, "\n")
  
  drugs_in_split <- drug_splits$drug[drug_splits[[paste0("split_", i)]]]
  drug_sample_ids <- get_drug_samples(kidney_eset, drugs_in_split)
  ctrl_sample_ids <- kidney_ctrl
  
  if (length(drug_sample_ids) < 3) {
    cat("  Skipping: not enough drug samples\n")
    next
  }
  
  cat("  Drugs:", length(drugs_in_split),
      "Drug samples:", length(drug_sample_ids), 
      "Control samples:", length(ctrl_sample_ids), "\n")
  
  res <- run_limma(
    eset = kidney_eset,
    group1_ids = drug_sample_ids,
    group2_ids = ctrl_sample_ids,
    group1_name = "drug",
    group2_name = "control"
  )
  
  if (is.null(res)) {
    cat("  Skipping: limma failed\n")
    next
  }
  
  res_df <- res %>%
    as_tibble(rownames = "gene_id") %>%
    arrange(adj.P.Val, desc(abs(logFC)))
  
  sig <- extract_signatures(res, gene_symbol_col_in_table = "Gene.Symbol")
  sig_kidney_1_10th[[paste0("split_", i)]] <- sig
  
  if (do_save) {
    write.csv(res_df,
              file.path(PATH, paste0("data/sigs/drugmatrix/mean/10th_perturb/kidney_10per_", i, ".csv")),
              row.names = FALSE)
  }
}

# Save signatures
if (do_save) {
  saveRDS(sig_liver_1_10th,
          file.path(PATH, "data/sigs/drugmatrix/mean/liver_10_signatures.rds"))
  saveRDS(sig_kidney_1_10th,
          file.path(PATH, "data/sigs/drugmatrix/mean/kidney_10_signatures.rds"))
}

# ============================================================================
# Analysis 3: 9/10th Perturbed vs All Controls (10 splits)
# ============================================================================

cat("\n=== Analysis 3: 9/10th Perturbed vs Controls ===\n")

# Liver - 9/10th splits
sig_liver_9_10th <- list()

for (i in 1:n_splits) {
  cat("Processing liver split", i, "of", n_splits, "(excluding split", i, ")\n")
  
  # Get drugs NOT in this split (i.e., 9/10th of drugs)
  drugs_in_9_10th <- drug_splits$drug[!drug_splits[[paste0("split_", i)]]]
  
  drug_sample_ids <- get_drug_samples(liver_eset, drugs_in_9_10th)
  ctrl_sample_ids <- liver_ctrl
  
  if (length(drug_sample_ids) < 3) {
    cat("  Skipping: not enough drug samples\n")
    next
  }
  
  cat("  Drugs:", length(drugs_in_9_10th),
      "Drug samples:", length(drug_sample_ids), 
      "Control samples:", length(ctrl_sample_ids), "\n")
  
  res <- run_limma(
    eset = liver_eset,
    group1_ids = drug_sample_ids,
    group2_ids = ctrl_sample_ids,
    group1_name = "drug",
    group2_name = "control"
  )
  
  if (is.null(res)) {
    cat("  Skipping: limma failed\n")
    next
  }
  
  res_df <- res %>%
    as_tibble(rownames = "gene_id") %>%
    arrange(adj.P.Val, desc(abs(logFC)))
  
  sig <- extract_signatures(res, gene_symbol_col_in_table = "Gene.Symbol")
  sig_liver_9_10th[[paste0("split_", i, "_excluded")]] <- sig
  
  if (do_save) {
    write.csv(res_df,
              file.path(PATH, paste0("data/sigs/drugmatrix/mean/90th_perturb/liver_90_", i, ".csv")),
              row.names = FALSE)
  }
  
  cat("  Up:", length(sig$up), "Down:", length(sig$dn), "\n")
}

# Kidney - 9/10th splits
sig_kidney_9_10th <- list()

for (i in 1:n_splits) {
  cat("Processing kidney split", i, "of", n_splits, "(excluding split", i, ")\n")
  
  drugs_in_9_10th <- drug_splits$drug[!drug_splits[[paste0("split_", i)]]]
  drug_sample_ids <- get_drug_samples(kidney_eset, drugs_in_9_10th)
  ctrl_sample_ids <- kidney_ctrl
  
  if (length(drug_sample_ids) < 3) {
    cat("  Skipping: not enough drug samples\n")
    next
  }
  
  cat("  Drugs:", length(drugs_in_9_10th),
      "Drug samples:", length(drug_sample_ids), 
      "Control samples:", length(ctrl_sample_ids), "\n")
  
  res <- run_limma(
    eset = kidney_eset,
    group1_ids = drug_sample_ids,
    group2_ids = ctrl_sample_ids,
    group1_name = "drug",
    group2_name = "control"
  )
  
  if (is.null(res)) {
    cat("  Skipping: limma failed\n")
    next
  }
  
  res_df <- res %>%
    as_tibble(rownames = "gene_id") %>%
    arrange(adj.P.Val, desc(abs(logFC)))
  
  sig <- extract_signatures(res, gene_symbol_col_in_table = "Gene.Symbol")
  sig_kidney_9_10th[[paste0("split_", i, "_excluded")]] <- sig
  
  if (do_save) {
    write.csv(res_df,
              file.path(PATH, paste0("data/sigs/drugmatrix/mean/90th_perturb/kidney_90_", i, ".csv")),
              row.names = FALSE)
  }
  
  cat("  Up:", length(sig$up), "Down:", length(sig$dn), "\n")
}

# Save signatures
if (do_save) {
  saveRDS(sig_liver_9_10th,
          file.path(PATH, "data/sigs/drugmatrix/mean/liver_90_signatures.rds"))
  saveRDS(sig_kidney_9_10th,
          file.path(PATH, "data/sigs/drugmatrix/mean/kidney_90_signatures.rds"))
}
