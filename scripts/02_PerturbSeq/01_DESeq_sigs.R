library(DESeq2)
library(tidyverse)
library(doParallel)
library(sigrecon)
registerDoParallel(cores=15)

PATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/projects/sig_recon/data/perturb_seq")
savepath <- file.path(Sys.getenv("MLAB"), "projects/brcameta/projects/sig_recon/data/sigs/perturb-seq")

k562_counts <- read.csv(file.path(PATH, "k562_processed_pb.csv"), row.names = 1)
colnames(k562_counts) <- str_replace_all(colnames(k562_counts),
                                         pattern = "\\.",
                                         replacement = "-")
k562_meta <- read.csv(file.path(PATH, "k562_processed_pb_metadata.csv"), row.names = 1)
stopifnot(all(colnames(k562_counts) %in% rownames(k562_meta)))
k562_meta$gem_group <- as.factor(k562_meta$gem_group)
k562_meta$is_control <- ifelse(k562_meta$gene == "non-targeting", "NTC", "Perturbed")

# 1. K562 PB sigs
unique_targets <- unique(k562_meta$gene[k562_meta$gene != "non-targeting"])

sig_dfs <- foreach(target_gene = unique_targets, .combine = dplyr::bind_rows) %dopar% {
  subset_meta <- k562_meta %>% filter(gene == target_gene | is_control == "NTC")
  subset_counts <- k562_counts[, rownames(subset_meta)]

  dds <- DESeqDataSetFromMatrix(countData = subset_counts,
                                colData = subset_meta,
                                design = ~ gem_group + is_control)

  # Run DESeq
  dds <- DESeq(dds)
  res <- results(dds, contrast=c("is_control", "Perturbed", "NTC"))
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  res_df$pb <- target_gene
  rownames(res_df) <- NULL
  res_df$cell_line <- "k562"
  res_df <- res_df %>%
    tidyr::drop_na() %>%
    dplyr::mutate(logFC_pval = log2FoldChange*(-log10(padj))) %>%
    dplyr::arrange(desc(logFC_pval))
  res_df
}
saveRDS(sig_dfs, file.path(savepath, "k562_pb_deseq_tables.rds"))
rm(k562_counts)
gc()

# 2. RPE1 PB sigs
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

unique_targets <- unique(rpe1_meta$gene[rpe1_meta$gene != "non-targeting"])

sig_dfs <- foreach(target_gene = unique_targets, .combine = dplyr::bind_rows) %dopar% {
  subset_meta <- rpe1_meta %>% filter(gene == target_gene | is_control == "NTC")
  if(sum(subset_meta$gene == target_gene) < 2) {
    print(paste0("Skipping ", target_gene, "not enough samples."))
    NULL
  } else {
    subset_counts <- rpe1_counts[, rownames(subset_meta)]

    dds <- DESeqDataSetFromMatrix(countData = subset_counts,
                                  colData = subset_meta,
                                  design = ~ gem_group + is_control)

    # Run DESeq
    dds <- DESeq(dds)
    res <- results(dds, contrast=c("is_control", "Perturbed", "NTC"))
    res_df <- as.data.frame(res)
    res_df$gene <- rownames(res_df)
    res_df$pb <- target_gene
    rownames(res_df) <- NULL
    res_df$cell_line <- "rpe1"
    res_df <- res_df %>%
      tidyr::drop_na() %>%
      dplyr::mutate(logFC_pval = log2FoldChange*(-log10(padj))) %>%
      dplyr::arrange(desc(logFC_pval))
    res_df
  }
}
saveRDS(sig_dfs, file.path(savepath, "rpe1_pb_deseq_tables.rds"))
rm(rpe1_counts)
gc()

# Filtering sigs
rpe1_tbl <- readRDS(file.path(savepath, "rpe1_pb_deseq_tables.rds"))
k562_tbl <- readRDS(file.path(savepath, "k562_pb_deseq_tables.rds"))

shared_pbs <- intersect(rpe1_tbl$pb, k562_tbl$pb)
rpe1_sigs <- sig_filter_fn(rpe1_tbl, perts = shared_pbs, pert_col = "pb", log2fc_col = "log2FoldChange", pval_col = "padj", geneid_col = "gene")
k562_sigs <- sig_filter_fn(k562_tbl, perts = shared_pbs, pert_col = "pb", log2fc_col = "log2FoldChange", pval_col = "padj", geneid_col = "gene")
saveRDS(rpe1_sigs, file.path(savepath, "rpe1_filtered_sigs.rds"))
saveRDS(k562_sigs, file.path(savepath, "k562_filtered_sigs.rds"))
stopifnot(all.equal(names(rpe1_sigs), names(k562_sigs)))

perturbseq.k562 <- k562_sigs
perturbseq.rpe1 <- rpe1_sigs

usethis::use_data(perturbseq.k562, overwrite = TRUE)
usethis::use_data(perturbseq.rpe1, overwrite = TRUE)
