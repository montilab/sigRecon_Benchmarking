library(Seurat)
library(tidyverse)
library(sigrecon)
# library(stringdist)
DATA_PATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/projects/sig_recon/data/stack_splits/perturb_seq")

# 1. Loading data
rpe1_pred_df <- readRDS(file.path(DATA_PATH, "rpe1_sigs.rds"))
k562_pred_df <- readRDS(file.path(DATA_PATH, "k562_sigs.rds"))

rpe1_true_sigs <- sigrecon::perturbseq.rpe1
k562_true_sigs <- sigrecon::perturbseq.k562
stopifnot(names(rpe1_true_sigs) == names(k562_true_sigs))

rpe1_pred_df <- rpe1_pred_df %>% 
  tibble::rownames_to_column(var = "id") %>% 
  dplyr::mutate(hgnc = sub("\\.\\.\\..*$", "", id)) 
k562_pred_df <- k562_pred_df %>% 
  tibble::rownames_to_column(var = "id") %>% 
  dplyr::mutate(hgnc = sub("\\.\\.\\..*$", "", id)) 

# 2. Converting to named genelist
k562_rpe1_sigs <- sig_filter_fn(diff_table = rpe1_pred_df, perts = unique(rpe1_pred_df$product), pert_col = "product", log2fc_col = "avg_log2FC", pval_col = "p_val_adj", geneid_col = "hgnc")
saveRDS(k562_rpe1_sigs, file.path(DATA_PATH, "k562_rpe1_sigs.rds"))
rpe1_k562_sigs <- sig_filter_fn(diff_table = k562_pred_df, perts = unique(k562_pred_df$product), pert_col = "product", log2fc_col = "avg_log2FC", pval_col = "p_val_adj", geneid_col = "hgnc")
saveRDS(rpe1_k562_sigs, file.path(DATA_PATH, "rpe1_k562_sigs.rds"))

# Continue from here
k562_rpe1_sigs <- lapply(k562_rpe1_sigs, function(x) lapply(x, function(y) y[!is.na(y)]))
rpe1_k562_sigs <- lapply(rpe1_k562_sigs, function(x) lapply(x, function(y) y[!is.na(y)]))

k562_rpe1_sigs <- k562_rpe1_sigs[names(rpe1_true_sigs)]
rpe1_k562_sigs <- rpe1_k562_sigs[names(k562_true_sigs)]

# 3. eval_table
k562_rpe1_eval_table <- sig_eval_table(source_sigs = lapply(k562_true_sigs, function(x) x$up), 
                                       pred_sigs = lapply(k562_rpe1_sigs, function(x) x$up), 
                                       true_sigs = rpe1_true_sigs,
                                       source = "k562")
rpe1_k562_eval_table <- sig_eval_table(source_sigs = lapply(rpe1_true_sigs, function(x) x$up), 
                                       pred_sigs = lapply(rpe1_k562_sigs, function(x) x$up), 
                                       true_sigs = k562_true_sigs,
                                       source = "rpe1")
all_eval_table <- rbind(k562_rpe1_eval_table,
                        rpe1_k562_eval_table)
saveRDS(all_eval_table, file.path(DATA_PATH, "eval_table.rds"))

# Unpaired comparison
all_eval_table$NES %>% mean(., na.rm = TRUE)
all_eval_table$jacc %>% mean(., na.rm = TRUE)

k562_rpe1_eval_table_n <- sig_eval_table(source_sigs = lapply(k562_true_sigs, function(x) x$up), 
                                         pred_sigs = lapply(k562_true_sigs, function(x) x$up), 
                                         true_sigs = rpe1_true_sigs, source="k562")
rpe1_k562_eval_table_n <- sig_eval_table(source_sigs = lapply(rpe1_true_sigs, function(x) x$up), 
                                         pred_sigs = lapply(rpe1_true_sigs, function(x) x$up), 
                                         true_sigs = k562_true_sigs, source="rpe1")
all_eval_table_n <- rbind(k562_rpe1_eval_table_n,
                          rpe1_k562_eval_table_n)

combined_aggregated_df <- merge(all_eval_table_n %>% dplyr::select("source", "gene", "NES", "jacc"), 
                                all_eval_table %>% mutate(kept_alpha = kept/(kept+displaced)), 
                                by = c("source", "gene"), 
                                suffixes = c("_FALSE", "_TRUE")) %>% tibble

# Plots
combined_aggregated_df %>% 
  filter(source=="k562") %>%
  ggplot(aes(x=NES_FALSE, y=NES_TRUE)) +
  geom_point(aes(alpha=kept_alpha)) +
  geom_abline(color="black") +
  scale_alpha_continuous("Proportion of Source Kept") +
  labs(title = paste("KS Distances in", "K562 Up"), 
       x = "KS Distance (Before)",
       y = "KS Distance (After)")

combined_aggregated_df %>% 
  filter(source=="k562") %>%
  ggplot(aes(x=jacc_FALSE, y=jacc_TRUE)) +
  geom_point(aes(alpha=kept_alpha)) +
  geom_abline(color="black") +
  scale_alpha_continuous("Proportion of Source Kept") +
  labs(title = paste("Jaccard Similarity in", "K562 Up"), 
       x = "Jaccard (Before)",
       y = "Jaccard (After)")

combined_aggregated_df %>% 
  filter(source=="rpe1") %>%
  ggplot(aes(x=NES_FALSE, y=NES_TRUE)) +
  geom_point(aes(alpha=kept_alpha)) +
  geom_abline(color="black") +
  scale_alpha_continuous("Proportion of Source Kept") +
  labs(title = paste("KS Distances in", "RPE1 Up"), 
       x = "KS Distance (Before)",
       y = "KS Distance (After)")

combined_aggregated_df %>% 
  filter(source=="rpe1") %>%
  ggplot(aes(x=jacc_FALSE, y=jacc_TRUE)) +
  geom_point(aes(alpha=kept_alpha)) +
  geom_abline(color="black") +
  scale_alpha_continuous("Proportion of Source Kept") +
  labs(title = paste("Jaccard Similarity in", "RPE1 Up"), 
       x = "Jaccard (Before)",
       y = "Jaccard (After)")
# Paired Tests
combined_aggregated_df %>% 
  group_by(source) %>% 
  summarize(Jacc_Wilcox_pval = wilcox.test(jacc_TRUE, jacc_FALSE, paired = TRUE, alternative="greater")$p.value,
            KS_Wilcox_pval = wilcox.test(NES_TRUE, NES_FALSE, paired = TRUE, alternative="greater")$p.value) %>%
  summarize(KS_meta_p = fishers_meta_p(KS_Wilcox_pval),
            Jacc_meta_p = fishers_meta_p(Jacc_Wilcox_pval))
