library(tidyverse)
library(limma)
library(Biobase)

PATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/projects/sig_recon")
DATA_PATH <- file.path(Sys.getenv("AGED"), "CBMrepositoryData/perturbational_data/drugmatrix/")
do_save <- TRUE

# 0. Load Expression Sets
liver_eset <- readRDS(file=file.path(DATA_PATH, "liver.rds"))
kidney_eset <- readRDS(file=file.path(DATA_PATH, "kidney.rds"))

# Assay data is already log-transformed
liver_eset@assayData$exprs[1:1000,1:1000] %>% hist
liver_eset@assayData$exprs[1:1000,1:1000] %>% fivenum
cat("  Values > 20:", sum(liver_eset@assayData$exprs > 20), "\n")
# Helper Function 1: Filter Sample Metadata for drug-specific analysis (max dose/time per compound)
filter_sample_metadata <- function(eSet) {
  eSet %>%
    pData %>%
    filter(`dose:ch1` != "0 mg/kg") %>%
    select(title, source_name_ch1, `compound:ch1`,`dose:ch1`, `tissue:ch1`, `time:ch1`, `vehicle:ch1`) %>%
    mutate(
      new_dose = as.numeric(str_extract(`dose:ch1`, "\\d+")),
      new_time = as.numeric(str_extract(`time:ch1`, "\\d+"))
    ) %>%
    rownames_to_column("sample_id") %>%
    group_by(`compound:ch1`) %>%
    filter(new_dose == max(new_dose, na.rm = TRUE)) %>% # Filter for max dose
    filter(new_time == max(new_time, na.rm = TRUE)) %>% # Filter for max time
    ungroup() # Ungroup after filtering
}

# Helper Function 2: Perform Limma DE
run_limma_de <- function(eSet, drug_sample_ids, ctrl_sample_ids) {

  # Check for sufficient samples
  if (length(drug_sample_ids) < 2 || length(ctrl_sample_ids) < 2) {
    message("Skipping DE: Not enough control or drug samples (need at least 2 each).")
    return(data.frame()) # Return an empty tibble
  }

  # Subset the ExpressionSet
  eset_sub <- eSet[, c(ctrl_sample_ids, drug_sample_ids)]

  # Ensure sample names are unique and reflect groups (important for design matrix)
  sampleNames(eset_sub) <- c(paste0("Ctrl_", seq_along(ctrl_sample_ids)),
                             paste0("Drug_", seq_along(drug_sample_ids)))

  # Create design data frame
  design_df <- data.frame(sample = sampleNames(eset_sub),
                          group = c(rep("Ctrl", length(ctrl_sample_ids)),
                                    rep("Drug", length(drug_sample_ids))))
  group <- as.factor(design_df$group)

  # Create design matrix for Limma
  design <- model.matrix(~ 0 + group, data = design_df)
  colnames(design) <- levels(group) # Rename columns for clarity in contrasts

  # Fit linear model for each gene
  fit <- limma::lmFit(eset_sub, design)

  # Build comparison contrast
  cont.matrix <- limma::makeContrasts(Drug - Ctrl, levels = design)

  # Apply the contrast matrix and compute statistics
  fit2 <- limma::contrasts.fit(fit, cont.matrix)
  fit2 <- limma::eBayes(fit2, 0.05) # Empirical Bayes smoothing

  # Extract top differentially expressed genes
  # limma::topTable automatically includes gene annotation from fData(eSet) if available
  de_table <- limma::topTable(fit2, coef = 1, adjust = "fdr", number = Inf)

  return(de_table) # Return as a tibble for consistency with tidyverse
}

# Helper Function 3: Extract Signatures from a DE table
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

# 1. Drug-specific
liver_samples <- filter_sample_metadata(liver_eset)
kidney_samples <- filter_sample_metadata(kidney_eset)

# Shared Drugs
# liver_kidney_drugs <- intersect(liver_samples$`compound:ch1`, kidney_samples$`compound:ch1`)
liver_drugs <- unique(liver_samples$`compound:ch1`)

# Limma - Lung
if(do_save) {
  for(drug in liver_drugs) {
    drug_samples <- liver_samples %>% filter(`compound:ch1`== drug)
    drug_vehicle <- drug_samples %>% pull(`vehicle:ch1`)
    drug_time <- drug_samples %>% pull(`time:ch1`)
    stopifnot(length(unique(drug_time)) == 1)
    drug_vehicle <- unique(drug_vehicle)
    drug_time <- unique(drug_time)
    ctrl_samples <- liver_eset %>% pData %>% filter(`dose:ch1` == "0 mg/kg",
                                                    `vehicle:ch1` == drug_vehicle,
                                                    `time:ch1` == drug_time)
    ctrl_ids <- ctrl_samples %>% rownames
    drug_ids <- drug_samples %>% pull(sample_id)
    # At least 2x2
    de_table <- run_limma_de(eSet = liver_eset, drug_sample_ids = drug_ids, ctrl_sample_ids = ctrl_ids)
    if(nrow(de_table) == 0) {
      next
    }
    de_table <- de_table %>% dplyr::select(c("Gene.Symbol", "ENTREZ_GENE_ID", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val"))
    write.csv(de_table, file.path(PATH, paste0("data/sigs/drugmatrix/tables/", drug, "_liver_sigs.csv")), row.names = TRUE)
  }
}

sig_table_paths <- Sys.glob(file.path(PATH, "data/sigs/drugmatrix/tables/*_liver_sigs.csv"))
drugs <- str_extract(sig_table_paths, "(?<=tables/)(.*?)(?=_liver_sigs\\.csv)")
sig_tables <- lapply(sig_table_paths, read.csv)
names(sig_tables) <- drugs

# Extracting Signatures
sigs <- list()
for(drug in drugs) {
  table <- sig_tables[[drug]]
  print(drug)
  sigs[[drug]] <- extract_signatures(table, "Gene.Symbol")
}

saveRDS(sigs, file.path(PATH, "data/sigs/drugmatrix/liver_shortlist_sigs.rds"))

# Limma - Kidney
kidney_drugs <- unique(kidney_samples$`compound:ch1`)
kidney_drugs <- kidney_drugs[!is.na(kidney_drugs)]
if(do_save) {
  for(drug in kidney_drugs) {
    drug_samples <- kidney_samples %>% filter(`compound:ch1`== drug)
    drug_vehicle <- drug_samples %>% pull(`vehicle:ch1`)
    drug_time <- drug_samples %>% pull(`time:ch1`)
    stopifnot(length(unique(drug_time)) == 1)
    drug_vehicle <- unique(drug_vehicle)
    drug_time <- unique(drug_time)
    ctrl_samples <- kidney_eset %>% pData %>% filter(`dose:ch1` == "0 mg/kg",
                                                     `vehicle:ch1` == drug_vehicle,
                                                     `time:ch1` == drug_time)
    ctrl_ids <- ctrl_samples %>% rownames
    drug_ids <- drug_samples %>% pull(sample_id)
    
    de_table <- run_limma_de(eSet = kidney_eset, drug_sample_ids = drug_ids, ctrl_sample_ids = ctrl_ids)
    if(nrow(de_table) == 0) {
      next
    }
    de_table <- de_table %>% dplyr::select(c("Gene.Symbol", "ENTREZ_GENE_ID", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val"))
    write.csv(de_table, file.path(PATH, paste0("data/sigs/drugmatrix/tables/", drug, "_kidney_sigs.csv")), row.names = TRUE)
  }
}

sig_table_paths <- Sys.glob(file.path(PATH, "data/sigs/drugmatrix/tables/*_kidney_sigs.csv"))
drugs <- str_extract(sig_table_paths, "(?<=tables/)(.*?)(?=_kidney_sigs\\.csv)")
sig_tables <- lapply(sig_table_paths, read.csv)
names(sig_tables) <- drugs

# Extracting Signatures
sigs <- list()
for(drug in drugs) {
  table <- sig_tables[[drug]]
  print(drug)
  sigs[[drug]] <- extract_signatures(table, "Gene.Symbol")
}
saveRDS(sigs, file.path(PATH, "data/sigs/drugmatrix/kidney_shortlist_sigs.rds"))

kidney_sigs <- readRDS(file.path(PATH, "data/sigs/drugmatrix/kidney_shortlist_sigs.rds"))
liver_sigs <- readRDS(file.path(PATH, "data/sigs/drugmatrix/liver_shortlist_sigs.rds"))

kidney_sigs <- lapply(kidney_sigs, function(x) x[c("up", "up_full")])
liver_sigs <- lapply(liver_sigs, function(x) x[c("up", "up_full")])

shared_drugs <- intersect(names(kidney_sigs), names(liver_sigs))
kidney_sigs <- kidney_sigs[shared_drugs]
liver_sigs <- liver_sigs[shared_drugs]

drugmatrix.liver <- liver_sigs
drugmatrix.kidney <- kidney_sigs

usethis::use_data(drugmatrix.liver, overwrite= TRUE)
usethis::use_data(drugmatrix.kidney, overwrite = TRUE)