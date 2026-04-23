library(tidyverse)
library(DESeq2)
library(Seurat)
library(doParallel)
library(sigrecon)
registerDoParallel(cores=15)

PATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/projects/sig_recon")
savepath <- file.path(Sys.getenv("MLAB"), "projects/brcameta/projects/sig_recon/data/sigs/sciplex")
do_save <- FALSE

# 0. Pseudobulks
if(do_save) {
  a549 <- readRDS(file.path(PATH, "data/sci_plex/a549_filtered.rds"))
  k562 <- readRDS(file.path(PATH, "data/sci_plex/k562_filtered.rds"))
  mcf7 <- readRDS(file.path(PATH, "data/sci_plex/mcf7_filtered.rds"))
  
  a549 <- CreateSeuratObject(counts = a549@assays$RNA$data,
                             meta.data = a549@meta.data)
  k562 <- CreateSeuratObject(counts = k562@assays$RNA$data,
                             meta.data = k562@meta.data)
  mcf7 <- CreateSeuratObject(counts = mcf7@assays$RNA$data,
                             meta.data = mcf7@meta.data)
  gc()
  
  a549_drug_ids <- a549@meta.data %>% 
    tibble::rownames_to_column(var = "id") %>% 
    dplyr::group_by(product_name, dose) %>% 
    dplyr::mutate(dose = as.numeric(dose)) %>% 
    dplyr::filter(dose == max(dose)) %>% 
    dplyr::pull(id)
  a549_ctrls <- a549@meta.data %>% rownames_to_column(var= "id") %>% dplyr::filter(product_name == "Vehicle") %>% pull(id)
  a549_ids <- c(a549_ctrls, a549_drug_ids)
  
  k562_drug_ids <- k562@meta.data %>% 
    tibble::rownames_to_column(var = "id") %>% 
    dplyr::group_by(product_name, dose) %>% 
    dplyr::mutate(dose = as.numeric(dose)) %>% 
    dplyr::filter(dose == max(dose)) %>% 
    dplyr::pull(id)
  k562_ctrls <- k562@meta.data %>% rownames_to_column(var= "id") %>% dplyr::filter(product_name == "Vehicle") %>% pull(id)
  k562_ids <- c(k562_ctrls, k562_drug_ids)
  
  mcf7_drug_ids <- mcf7@meta.data %>% 
    tibble::rownames_to_column(var = "id") %>% 
    dplyr::group_by(product_name, dose) %>% 
    dplyr::mutate(dose = as.numeric(dose)) %>% 
    dplyr::filter(dose == max(dose)) %>% 
    dplyr::pull(id)
  mcf7_ctrls <- mcf7@meta.data %>% rownames_to_column(var= "id") %>% dplyr::filter(product_name == "Vehicle") %>% pull(id)
  mcf7_ids <- c(mcf7_ctrls, mcf7_drug_ids)
  
  a549_pb <- AggregateExpression(a549[,a549_ids], group.by = c("product_name", "replicate"), return.seurat = TRUE)
  k562_pb <- AggregateExpression(k562[,k562_ids], group.by = c("product_name", "replicate"), return.seurat = TRUE)
  mcf7_pb <- AggregateExpression(mcf7[,mcf7_ids], group.by = c("product_name", "replicate"), return.seurat = TRUE)
  
  saveRDS(a549_pb, file.path(PATH, "data/sci_plex/a549_filtered_pb.rds"))
  saveRDS(k562_pb, file.path(PATH, "data/sci_plex/k562_filtered_pb.rds"))
  saveRDS(mcf7_pb, file.path(PATH, "data/sci_plex/mcf7_filtered_pb.rds"))
  
} else {
  a549_pb <- readRDS(file.path(PATH, "data/sci_plex/a549_filtered_pb.rds"))
  k562_pb <- readRDS(file.path(PATH, "data/sci_plex/k562_filtered_pb.rds"))
  mcf7_pb <- readRDS(file.path(PATH, "data/sci_plex/mcf7_filtered_pb.rds"))
}

# 0. Mean Sigs
metadata <- a549_pb@meta.data
counts <- a549_pb@assays$RNA$counts
metadata$is_control <- ifelse(metadata$product_name == "Vehicle", "Control", "Perturbed")
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = metadata,
                              design = ~ replicate + is_control)

# Run DESeq
dds <- DESeq(dds)
res <- results(dds, contrast=c("is_control", "Perturbed", "Control"))
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
res_df$pb <- "mean"
rownames(res_df) <- NULL
res_df$cell_line <- "a549"
res_df <- res_df %>%
  tidyr::drop_na() %>%
  dplyr::mutate(logFC_pval = log2FoldChange*(-log10(padj))) %>%
  dplyr::arrange(desc(logFC_pval))
saveRDS(res_df, file.path(savepath, "mean/all/a549_mean_pb_deseq_tables.rds"))

metadata <- k562_pb@meta.data
counts <- k562_pb@assays$RNA$counts
metadata$is_control <- ifelse(metadata$product_name == "Vehicle", "Control", "Perturbed")
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = metadata,
                              design = ~ replicate + is_control)

# Run DESeq
dds <- DESeq(dds)
res <- results(dds, contrast=c("is_control", "Perturbed", "Control"))
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
res_df$pb <- "mean"
rownames(res_df) <- NULL
res_df$cell_line <- "k562"
res_df <- res_df %>%
  tidyr::drop_na() %>%
  dplyr::mutate(logFC_pval = log2FoldChange*(-log10(padj))) %>%
  dplyr::arrange(desc(logFC_pval))
saveRDS(res_df, file.path(savepath, "mean/all/k562_mean_pb_deseq_tables.rds"))

metadata <- mcf7_pb@meta.data
counts <- mcf7_pb@assays$RNA$counts
metadata$is_control <- ifelse(metadata$product_name == "Vehicle", "Control", "Perturbed")
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = metadata,
                              design = ~ replicate + is_control)

# Run DESeq
dds <- DESeq(dds)
res <- results(dds, contrast=c("is_control", "Perturbed", "Control"))
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
res_df$pb <- "mean"
rownames(res_df) <- NULL
res_df$cell_line <- "mcf7"
res_df <- res_df %>%
  tidyr::drop_na() %>%
  dplyr::mutate(logFC_pval = log2FoldChange*(-log10(padj))) %>%
  dplyr::arrange(desc(logFC_pval))
saveRDS(res_df, file.path(savepath, "mean/all/mcf7_mean_pb_deseq_tables.rds"))

# 1. A549 Sigs
a549_drugs <- unique(a549_pb$product_name)
a549_drugs <- a549_drugs[a549_drugs != "Vehicle"]
sig_dfs <- foreach(drug = a549_drugs, .combine = dplyr::bind_rows) %dopar% {
  metadata <- a549_pb@meta.data
  counts <- a549_pb@assays$RNA$counts
  subset_meta <- metadata %>% filter(product_name == drug | product_name == "Vehicle")
  subset_meta$is_control <- ifelse(subset_meta$product_name == "Vehicle", "Control", "Perturbed")
  subset_counts <- counts[, rownames(subset_meta)]

  # ✅ Check if we have both drug and vehicle
  unique_products <- unique(subset_meta$product_name)
  if (length(unique_products) < 2) {
    warning("Skipping ", drug_name, " in ", cell_line_name,
            ": no Vehicle controls or only one condition")
    next
  }

  # ✅ Check minimum sample size
  product_counts <- table(subset_meta$product_name)
  if (any(product_counts < 2)) {
    warning("Skipping ", drug_name, " in ", cell_line_name,
            ": insufficient replicates (need ≥2 per condition)")
    next
  }

  dds <- DESeqDataSetFromMatrix(countData = subset_counts,
                                colData = subset_meta,
                                design = ~ replicate + product_name)

  # Run DESeq
  dds <- DESeq(dds)
  res <- results(dds, contrast=c("product_name", drug, "Vehicle"))
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  res_df$pb <- drug
  rownames(res_df) <- NULL
  res_df$cell_line <- "a549"
  res_df <- res_df %>%
    tidyr::drop_na() %>%
    dplyr::mutate(logFC_pval = log2FoldChange*(-log10(padj))) %>%
    dplyr::arrange(desc(logFC_pval))
  res_df
}
saveRDS(sig_dfs, file.path(savepath, "a549_pb_deseq_tables.rds"))
gc()

# 2. K562 Sigs
k562_drugs <- unique(k562_pb$product_name)
k562_drugs <- k562_drugs[k562_drugs != "Vehicle"]
sig_dfs <- foreach(drug = k562_drugs, .combine = dplyr::bind_rows) %dopar% {
  metadata <- k562_pb@meta.data
  counts <- k562_pb@assays$RNA$counts
  subset_meta <- metadata %>% filter(product_name == drug | product_name == "Vehicle")
  subset_meta$is_control <- ifelse(subset_meta$product_name == "Vehicle", "Control", "Perturbed")
  subset_counts <- counts[, rownames(subset_meta)]

  # ✅ Check if we have both drug and vehicle
  unique_products <- unique(subset_meta$product_name)
  if (length(unique_products) < 2) {
    warning("Skipping ", drug_name, " in ", cell_line_name,
            ": no Vehicle controls or only one condition")
    next
  }

  # ✅ Check minimum sample size
  product_counts <- table(subset_meta$product_name)
  if (any(product_counts < 2)) {
    warning("Skipping ", drug_name, " in ", cell_line_name,
            ": insufficient replicates (need ≥2 per condition)")
    next
  }
  dds <- DESeqDataSetFromMatrix(countData = subset_counts,
                                colData = subset_meta,
                                design = ~ replicate + product_name)

  # Run DESeq
  dds <- DESeq(dds)
  res <- results(dds, contrast=c("product_name", drug, "Vehicle"))
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  res_df$pb <- drug
  rownames(res_df) <- NULL
  res_df$cell_line <- "k562"
  res_df <- res_df %>%
    tidyr::drop_na() %>%
    dplyr::mutate(logFC_pval = log2FoldChange*(-log10(padj))) %>%
    dplyr::arrange(desc(logFC_pval))
  res_df
}
saveRDS(sig_dfs, file.path(savepath, "k562_pb_deseq_tables.rds"))
gc()

# 3. mcf7 Sigs
mcf7_drugs <- unique(mcf7_pb$product_name)
mcf7_drugs <- mcf7_drugs[mcf7_drugs != "Vehicle"]
sig_dfs <- foreach(drug = mcf7_drugs, .combine = dplyr::bind_rows) %dopar% {
  metadata <- mcf7_pb@meta.data
  counts <- mcf7_pb@assays$RNA$counts
  subset_meta <- metadata %>% filter(product_name == drug | product_name == "Vehicle")
  subset_meta$is_control <- ifelse(subset_meta$product_name == "Vehicle", "Control", "Perturbed")
  subset_counts <- counts[, rownames(subset_meta)]
  # ✅ Check if we have both drug and vehicle
  unique_products <- unique(subset_meta$product_name)
  if (length(unique_products) < 2) {
    warning("Skipping ", drug_name, " in ", cell_line_name,
            ": no Vehicle controls or only one condition")
    next
  }

  # ✅ Check minimum sample size
  product_counts <- table(subset_meta$product_name)
  if (any(product_counts < 2)) {
    warning("Skipping ", drug_name, " in ", cell_line_name,
            ": insufficient replicates (need ≥2 per condition)")
    next
  }
  dds <- DESeqDataSetFromMatrix(countData = subset_counts,
                                colData = subset_meta,
                                design = ~ replicate + product_name)

  # Run DESeq
  dds <- DESeq(dds)
  res <- results(dds, contrast=c("product_name", drug, "Vehicle"))
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  res_df$pb <- drug
  rownames(res_df) <- NULL
  res_df$cell_line <- "mcf7"
  res_df <- res_df %>%
    tidyr::drop_na() %>%
    dplyr::mutate(logFC_pval = log2FoldChange*(-log10(padj))) %>%
    dplyr::arrange(desc(logFC_pval))
  res_df
}
saveRDS(sig_dfs, file.path(savepath, "mcf7_pb_deseq_tables.rds"))
gc()

# Filtering for significant signatures
a549_tbls <- readRDS(file.path(savepath, "a549_pb_deseq_tables.rds"))
k562_tbls <- readRDS(file.path(savepath, "k562_pb_deseq_tables.rds"))
mcf7_tbls <- readRDS(file.path(savepath, "mcf7_pb_deseq_tables.rds"))

shared_pbs <- intersect(intersect(mcf7_tbls$pb, k562_tbls$pb), a549_tbls$pb)
mcf7_sig_list <- sig_filter_fn(mcf7_tbls, perts=shared_pbs, pert_col = "pb", log2fc_col = "log2FoldChange", pval_col = "padj", geneid_col = "gene")
k562_sig_list <- sig_filter_fn(k562_tbls, perts=shared_pbs, pert_col = "pb", log2fc_col = "log2FoldChange", pval_col = "padj", geneid_col = "gene")
a549_sig_list <- sig_filter_fn(a549_tbls, perts=shared_pbs, pert_col = "pb", log2fc_col = "log2FoldChange", pval_col = "padj", geneid_col = "gene")
# Removing extra string padding
mcf7_sig_list <- lapply(mcf7_sig_list, function(x) lapply(x, function(y) str_remove(string = y, pattern = "\\.\\.\\..+")))
k562_sig_list <- lapply(k562_sig_list, function(x) lapply(x, function(y) str_remove(string = y, pattern = "\\.\\.\\..+")))
a549_sig_list <- lapply(a549_sig_list, function(x) lapply(x, function(y) str_remove(string = y, pattern = "\\.\\.\\..+")))
# Removing NAs
mcf7_sig_list <- lapply(mcf7_sig_list, function(x) lapply(x, function(y) y[!is.na(y)]))
k562_sig_list <- lapply(k562_sig_list, function(x) lapply(x, function(y) y[!is.na(y)]))
a549_sig_list <- lapply(a549_sig_list, function(x) lapply(x, function(y) y[!is.na(y)]))

# Filtering to perturbations with a reasonable number of DEGS (arbitrarily set as 5)
mcf7_sig_list <- mcf7_sig_list[lapply(mcf7_sig_list, function(x) (length(x$up) >= 5)) %>% unlist]
k562_sig_list <- k562_sig_list[lapply(k562_sig_list, function(x) (length(x$up) >= 5)) %>% unlist]
a549_sig_list <- a549_sig_list[lapply(a549_sig_list, function(x) (length(x$up) >= 5)) %>% unlist]
lapply(mcf7_sig_list, function(x) length(x$up)) %>% unlist %>% unname %>% fivenum
lapply(k562_sig_list, function(x) length(x$up)) %>% unlist %>% unname %>% fivenum
lapply(a549_sig_list, function(x) length(x$up)) %>% unlist %>% unname %>% fivenum
lapply(mcf7_sig_list, function(x) length(x$up_full)) %>% unlist %>% unname %>% fivenum
lapply(k562_sig_list, function(x) length(x$up_full)) %>% unlist %>% unname %>% fivenum
lapply(a549_sig_list, function(x) length(x$up_full)) %>% unlist %>% unname %>% fivenum

saveRDS(mcf7_sig_list, file.path(savepath, "mcf7_sigs_filtered.rds"))
saveRDS(k562_sig_list, file.path(savepath, "k562_sigs_filtered.rds"))
saveRDS(a549_sig_list, file.path(savepath, "a549_sigs_filtered.rds"))

mcf7_sig_list <- readRDS(file.path(savepath, "mcf7_sigs_filtered.rds"))
k562_sig_list <- readRDS(file.path(savepath, "k562_sigs_filtered.rds"))
a549_sig_list <- readRDS(file.path(savepath, "a549_sigs_filtered.rds"))
sciplex.mcf7 <- mcf7_sig_list
sciplex.k562 <- k562_sig_list
sciplex.a549 <- a549_sig_list

usethis::use_data(sciplex.mcf7, overwrite = TRUE)
usethis::use_data(sciplex.k562, overwrite = TRUE)
usethis::use_data(sciplex.a549, overwrite = TRUE)

# Mean Signatures
a549_sigs <- readRDS(file.path(savepath, "mean/all/a549_mean_pb_deseq_tables.rds"))
k562_sigs <- readRDS(file.path(savepath, "mean/all/k562_mean_pb_deseq_tables.rds"))
mcf7_sigs <- readRDS(file.path(savepath, "mean/all/mcf7_mean_pb_deseq_tables.rds"))
# for mean sigs, nothing is significant min padj is 0.99
mcf7_sig_list <- sig_filter_fn(mcf7_sigs, "mean", pert_col = "pb", log2fc_col = "log2FoldChange", pval_col = "padj", alpha = 1, geneid_col = "gene")
k562_sig_list <- sig_filter_fn(k562_sigs, "mean", pert_col = "pb", log2fc_col = "log2FoldChange", pval_col = "padj", alpha = 1, geneid_col = "gene")
a549_sig_list <- sig_filter_fn(a549_sigs, "mean", pert_col = "pb", log2fc_col = "log2FoldChange", pval_col = "padj", alpha = 1, geneid_col = "gene")
# Removing extra string padding
mcf7_sig_list <- lapply(mcf7_sig_list, function(x) lapply(x, function(y) str_remove(string = y, pattern = "\\.\\.\\..+")))
k562_sig_list <- lapply(k562_sig_list, function(x) lapply(x, function(y) str_remove(string = y, pattern = "\\.\\.\\..+")))
a549_sig_list <- lapply(a549_sig_list, function(x) lapply(x, function(y) str_remove(string = y, pattern = "\\.\\.\\..+")))
# Removing NAs
mcf7_sig_list <- lapply(mcf7_sig_list, function(x) lapply(x, function(y) y[!is.na(y)]))
k562_sig_list <- lapply(k562_sig_list, function(x) lapply(x, function(y) y[!is.na(y)]))
a549_sig_list <- lapply(a549_sig_list, function(x) lapply(x, function(y) y[!is.na(y)]))

# Filtering to perturbations with a reasonable number of DEGS (arbitrarily set as 5)
mcf7_sig_list <- mcf7_sig_list[lapply(mcf7_sig_list, function(x) (length(x$up) >= 5)) %>% unlist]
k562_sig_list <- k562_sig_list[lapply(k562_sig_list, function(x) (length(x$up) >= 5)) %>% unlist]
a549_sig_list <- a549_sig_list[lapply(a549_sig_list, function(x) (length(x$up) >= 5)) %>% unlist]
saveRDS(mcf7_sig_list, file.path(savepath, "mean/all/mcf7_sigs_filtered.rds"))
saveRDS(k562_sig_list, file.path(savepath, "mean/all/k562_sigs_filtered.rds"))
saveRDS(a549_sig_list, file.path(savepath, "mean/all/a549_sigs_filtered.rds"))