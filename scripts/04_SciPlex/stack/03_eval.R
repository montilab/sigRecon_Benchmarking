library(Seurat)
library(tidyverse)
library(sigrecon)
library(stringdist)
DATA_PATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/projects/sig_recon/data/sigs/sciplex/stack")
SAVE_PATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/projects/sig_recon/results/eval/sciplex")
SC_PATH <- file.path(Sys.getenv("AGED"), "CBMrepositoryData/perturbational_data/srivatsan_2019")

# 1. Loading data
a549 <- readRDS(file.path(SC_PATH, "a549.rds"))
# meta.features are all equivalent
# k562 <- readRDS(file.path(SC_PATH, "k562.rds"))
# mcf7 <- readRDS(file.path(SC_PATH, "mcf7.rds"))
a549_pred_df <- readRDS(file.path(DATA_PATH, "a549_ctrl_sigs.rds"))
k562_pred_df <- readRDS(file.path(DATA_PATH, "k562_ctrl_sigs.rds"))
mcf7_pred_df <- readRDS(file.path(DATA_PATH, "mcf7_ctrl_sigs.rds"))

a549_true_sigs <- sigrecon::sciplex.a549
k562_true_sigs <- sigrecon::sciplex.k562
mcf7_true_sigs <- sigrecon::sciplex.mcf7
stopifnot(names(a549_true_sigs) == names(k562_true_sigs))
stopifnot(names(mcf7_true_sigs) == names(k562_true_sigs))

## Converting hgnc to ensembl
# There's duplicate genes here
ensembl_hgnc_tbl <- a549@assays$RNA@meta.features %>%
  tibble::rownames_to_column("ensembl_id")
a549_pred_df <- a549_pred_df %>% 
  dplyr::left_join(ensembl_hgnc_tbl %>% select(feature_name, ensembl_id), by = c("gene" = "feature_name")) %>% 
  dplyr::group_by(pvalue, log2FoldChange, product) %>%
  tidyr::fill(ensembl_id, .direction = "downup") %>%
  dplyr::ungroup() %>%
  dplyr::left_join(ensembl_hgnc_tbl %>% select(feature_name, ensembl_id), by = "ensembl_id")
k562_pred_df <- k562_pred_df %>% 
  dplyr::left_join(ensembl_hgnc_tbl %>% select(feature_name, ensembl_id), by = c("gene" = "feature_name")) %>% 
  dplyr::group_by(pvalue, log2FoldChange, product) %>%
  tidyr::fill(ensembl_id, .direction = "downup") %>%
  dplyr::ungroup() %>%
  dplyr::left_join(ensembl_hgnc_tbl %>% select(feature_name, ensembl_id), by = "ensembl_id")
mcf7_pred_df <- mcf7_pred_df %>% 
  dplyr::left_join(ensembl_hgnc_tbl %>% select(feature_name, ensembl_id), by = c("gene" = "feature_name")) %>% 
  dplyr::group_by(pvalue, log2FoldChange, product) %>%
  tidyr::fill(ensembl_id, .direction = "downup") %>%
  dplyr::ungroup() %>%
  dplyr::left_join(ensembl_hgnc_tbl %>% select(feature_name, ensembl_id), by = "ensembl_id")

# 2. Converting to named genelist
a549_mcf7_sigs <- a549_pred_df %>% 
  dplyr::filter(target_cell_line == "mcf7") %>% 
  sig_filter_fn(., perts = unique(a549_pred_df$product), pert_col = "product", log2fc_col = "log2FoldChange", pval_col = "padj", alpha = 1.1, geneid_col = "ensembl_id")
a549_k562_sigs <- a549_pred_df %>% 
  dplyr::filter(target_cell_line == "k562") %>% 
  sig_filter_fn(., perts = unique(a549_pred_df$product), pert_col = "product", log2fc_col = "log2FoldChange", pval_col = "padj", alpha = 1.1, geneid_col = "ensembl_id")
k562_mcf7_sigs <- k562_pred_df %>% 
  dplyr::filter(target_cell_line == "mcf7") %>% 
  sig_filter_fn(., perts = unique(k562_pred_df$product), pert_col = "product", log2fc_col = "log2FoldChange", pval_col = "padj", alpha = 1.1, geneid_col = "ensembl_id")
k562_a549_sigs <- k562_pred_df %>% 
  dplyr::filter(target_cell_line == "a549") %>% 
  sig_filter_fn(., perts = unique(k562_pred_df$product), pert_col = "product", log2fc_col = "log2FoldChange", pval_col = "padj", alpha = 1.1, geneid_col = "ensembl_id")
mcf7_a549_sigs <- mcf7_pred_df %>% 
  dplyr::filter(target_cell_line == "a549") %>% 
  sig_filter_fn(., perts = unique(mcf7_pred_df$product), pert_col = "product", log2fc_col = "log2FoldChange", pval_col = "padj", alpha = 1.1, geneid_col = "ensembl_id")
mcf7_k562_sigs <- mcf7_pred_df %>% 
  dplyr::filter(target_cell_line == "k562") %>% 
  sig_filter_fn(., perts = unique(mcf7_pred_df$product), pert_col = "product", log2fc_col = "log2FoldChange", pval_col = "padj", alpha = 1.1, geneid_col = "ensembl_id")

a549_mcf7_sigs <- lapply(a549_mcf7_sigs, function(x) lapply(x, function(y) y[!is.na(y)]))
a549_k562_sigs <- lapply(a549_k562_sigs, function(x) lapply(x, function(y) y[!is.na(y)]))
k562_mcf7_sigs <- lapply(k562_mcf7_sigs, function(x) lapply(x, function(y) y[!is.na(y)]))
k562_a549_sigs <- lapply(k562_a549_sigs, function(x) lapply(x, function(y) y[!is.na(y)]))
mcf7_a549_sigs <- lapply(mcf7_a549_sigs, function(x) lapply(x, function(y) y[!is.na(y)]))
mcf7_k562_sigs <- lapply(mcf7_k562_sigs, function(x) lapply(x, function(y) y[!is.na(y)]))
stopifnot(names(a549_mcf7_sigs) == names(a549_k562_sigs))
stopifnot(names(k562_mcf7_sigs) == names(a549_k562_sigs))
stopifnot(names(k562_mcf7_sigs) == names(k562_a549_sigs))
stopifnot(names(k562_a549_sigs) == names(mcf7_a549_sigs))
stopifnot(names(mcf7_k562_sigs) == names(mcf7_a549_sigs))

# fixing names
missing_names <- names(a549_mcf7_sigs)[which(!(names(a549_mcf7_sigs) %in% names(a549_true_sigs)))]
# Find potential matches for review
matches_df <- data.frame(
  stack_name = missing_names,
  best_match = NA,
  distance = NA
)

for (i in seq_along(missing_names)) {
  distances <- stringdist(missing_names[i], names(a549_true_sigs), method = "lv")
  best_idx <- which.min(distances)
  matches_df$best_match[i] <- names(a549_true_sigs)[best_idx]
  matches_df$distance[i] <- min(distances)
}

View(matches_df %>% arrange(distance))
matches_df <- matches_df %>% arrange(distance) %>% dplyr::slice(c(1:8,10:13,15,18,20:21))
# Apply after reviewing
idx <- match(names(a549_mcf7_sigs), matches_df$stack_name)
names(a549_mcf7_sigs)[!is.na(idx)] <- matches_df$best_match[na.omit(idx)]
names(a549_k562_sigs)[!is.na(idx)] <- matches_df$best_match[na.omit(idx)]
names(k562_mcf7_sigs)[!is.na(idx)] <- matches_df$best_match[na.omit(idx)]
names(k562_a549_sigs)[!is.na(idx)] <- matches_df$best_match[na.omit(idx)]
names(mcf7_a549_sigs)[!is.na(idx)] <- matches_df$best_match[na.omit(idx)]
names(mcf7_k562_sigs)[!is.na(idx)] <- matches_df$best_match[na.omit(idx)]

a549_mcf7_sigs <- a549_mcf7_sigs[names(a549_true_sigs)]
a549_k562_sigs <- a549_k562_sigs[names(a549_true_sigs)]
k562_mcf7_sigs <- k562_mcf7_sigs[names(a549_true_sigs)]
k562_a549_sigs <- k562_a549_sigs[names(a549_true_sigs)]
mcf7_a549_sigs <- mcf7_a549_sigs[names(a549_true_sigs)]
mcf7_k562_sigs <- mcf7_k562_sigs[names(a549_true_sigs)]

# 3. eval_table
a549_mcf7_eval_table <- sig_eval_table(source_sigs = lapply(a549_true_sigs, function(x) x$up), 
                                       pred_sigs = lapply(a549_mcf7_sigs, function(x) x$up), 
                                       true_sigs = mcf7_true_sigs, 
                                       source = "a549_mcf7")
a549_k562_eval_table <- sig_eval_table(source_sigs = lapply(a549_true_sigs, function(x) x$up), 
                                       pred_sigs = lapply(a549_k562_sigs, function(x) x$up), 
                                       true_sigs = k562_true_sigs,
                                       source = "a549_k562")
k562_mcf7_eval_table <- sig_eval_table(source_sigs = lapply(k562_true_sigs, function(x) x$up), 
                                       pred_sigs = lapply(k562_mcf7_sigs, function(x) x$up), 
                                       true_sigs = mcf7_true_sigs,
                                       source = "k562_mcf7")
k562_a549_eval_table <- sig_eval_table(source_sigs = lapply(k562_true_sigs, function(x) x$up), 
                                       pred_sigs = lapply(k562_a549_sigs, function(x) x$up), 
                                       true_sigs = a549_true_sigs,
                                       source = "k562_a549")
mcf7_a549_eval_table <- sig_eval_table(source_sigs = lapply(mcf7_true_sigs, function(x) x$up), 
                                       pred_sigs = lapply(mcf7_a549_sigs, function(x) x$up), 
                                       true_sigs = a549_true_sigs,
                                       source = "mcf7_a549")
mcf7_k562_eval_table <- sig_eval_table(source_sigs = lapply(mcf7_true_sigs, function(x) x$up), 
                                       pred_sigs = lapply(mcf7_k562_sigs, function(x) x$up), 
                                       true_sigs = k562_true_sigs,
                                       source = "mcf7_k562")

# no change
a549_mcf7_eval_table_n <- sig_eval_table(source_sigs = lapply(a549_true_sigs, function(x) x$up), 
                                       pred_sigs = lapply(a549_true_sigs, function(x) x$up), 
                                       true_sigs = mcf7_true_sigs, 
                                       source = "a549_mcf7")
a549_k562_eval_table_n <- sig_eval_table(source_sigs = lapply(a549_true_sigs, function(x) x$up), 
                                       pred_sigs = lapply(a549_true_sigs, function(x) x$up), 
                                       true_sigs = k562_true_sigs,
                                       source = "a549_k562")
k562_mcf7_eval_table_n <- sig_eval_table(source_sigs = lapply(k562_true_sigs, function(x) x$up), 
                                       pred_sigs = lapply(k562_true_sigs, function(x) x$up), 
                                       true_sigs = mcf7_true_sigs,
                                       source = "k562_mcf7")
k562_a549_eval_table_n <- sig_eval_table(source_sigs = lapply(k562_true_sigs, function(x) x$up), 
                                       pred_sigs = lapply(k562_true_sigs, function(x) x$up), 
                                       true_sigs = a549_true_sigs,
                                       source = "k562_a549")
mcf7_a549_eval_table_n <- sig_eval_table(source_sigs = lapply(mcf7_true_sigs, function(x) x$up), 
                                       pred_sigs = lapply(mcf7_true_sigs, function(x) x$up), 
                                       true_sigs = a549_true_sigs,
                                       source = "mcf7_a549")
mcf7_k562_eval_table_n <- sig_eval_table(source_sigs = lapply(mcf7_true_sigs, function(x) x$up), 
                                       pred_sigs = lapply(mcf7_true_sigs, function(x) x$up), 
                                       true_sigs = k562_true_sigs,
                                       source = "mcf7_k562")
all_eval_table <- rbind(a549_mcf7_eval_table,
                        a549_k562_eval_table,
                        k562_mcf7_eval_table,
                        k562_a549_eval_table,
                        mcf7_a549_eval_table,
                        mcf7_k562_eval_table)
all_eval_no_change_table <- rbind(a549_mcf7_eval_table_n,
                                  a549_k562_eval_table_n,
                                  k562_mcf7_eval_table_n,
                                  k562_a549_eval_table_n,
                                  mcf7_a549_eval_table_n,
                                  mcf7_k562_eval_table_n)
combined_aggregated_df <- merge(all_eval_no_change_table %>% dplyr::select("source", "gene", "NES", "jacc"), 
                                all_eval_table %>% mutate(kept_alpha = kept/(kept+displaced)), 
                                by = c("source", "gene"), 
                                suffixes = c("_FALSE", "_TRUE")) %>% tibble
saveRDS(combined_aggregated_df, file.path(SAVE_PATH, "stack_ctrl_eval_table.rds"))

# Unpaired 
combined_aggregated_df %>% 
  dplyr::group_by(source) %>%
  dplyr::summarize(NES_mean = mean(NES_TRUE, na.rm = TRUE),
                   jacc_mean = mean(jacc_TRUE)) %>%
  dplyr::summarize(NES_mean = mean(NES_mean),
                   jacc_mean = mean(jacc_mean))

# Plots
combined_aggregated_df %>% 
  filter(source=="k562_mcf7") %>%
  ggplot(aes(x=NES_FALSE, y=NES_TRUE)) +
  geom_point(aes(alpha=kept_alpha)) +
  geom_abline(color="black") +
  scale_alpha_continuous("Proportion of Source Kept") +
  labs(title = paste("KS Distances in", "K562 Up"), 
       x = "KS Distance (Before)",
       y = "KS Distance (After)")

combined_aggregated_df %>% 
  filter(source=="k562_mcf7") %>%
  ggplot(aes(x=jacc_FALSE, y=jacc_TRUE)) +
  geom_point(aes(alpha=kept_alpha)) +
  geom_abline(color="black") +
  scale_alpha_continuous("Proportion of Source Kept") +
  labs(title = paste("Jaccard Similarity in", "K562 Up"), 
       x = "Jaccard (Before)",
       y = "Jaccard (After)")


# Paired Tests
combined_aggregated_df %>% 
  group_by(source) %>% 
  summarize(Jacc_Wilcox_pval = wilcox.test(jacc_TRUE, jacc_FALSE, paired = TRUE, alternative="greater")$p.value,
            KS_Wilcox_pval = wilcox.test(NES_TRUE, NES_FALSE, paired = TRUE, alternative="greater")$p.value) %>%
  summarize(KS_meta_p = fishers_meta_p(KS_Wilcox_pval),
            Jacc_meta_p = fishers_meta_p(Jacc_Wilcox_pval))
