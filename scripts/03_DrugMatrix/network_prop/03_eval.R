library(tidyverse)
library(sigrecon)
library(BiocParallel)
library(SummarizedExperiment)
bp <- make_bpparam(workers = 15, RNGseed = 123, type="multicore")
register(bp)

PATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/projects/sig_recon")
SE_PATH <- file.path(Sys.getenv("AGED"), "CBMrepositoryData/perturbational_data/drugmatrix/")
RECON_PATH <- file.path(PATH, "data/sigs/drugmatrix/networkprop/")
SAVE_PATH <- file.path(PATH, "results/eval/drugmatrix")

# Loading sigs
liver_true_sigs <- drugmatrix.liver
kidney_true_sigs <- drugmatrix.kidney

# # Loading esets
# liver_eset  <- readRDS(file.path(DATA_PATH, "liver.rds"))
# kidney_eset <- readRDS(file.path(DATA_PATH, "kidney.rds"))
# # Changing eset features to gene symbols to match signatures
# liver_gene_symbols <- make.unique(fData(liver_eset)$`Gene Symbol`)
# featureNames(liver_eset) <- liver_gene_symbols
# kidney_gene_symbols <- make.unique(fData(kidney_eset)$`Gene Symbol`)
# featureNames(kidney_eset) <- kidney_gene_symbols
# liver_se <- as(liver_eset, "SummarizedExperiment")
# kidney_se <- as(kidney_eset, "SummarizedExperiment")

# Loading drug splits
drug_splits_path <- file.path(PATH, "data/sigs/drugmatrix/drug_splits.csv")

# # No Change
# liver_kidney_eval_table_n <- sig_eval_table(source_sigs = lapply(liver_true_sigs, function(x) x$up), 
#                                             pred_sigs = lapply(liver_true_sigs, function(x) x$up), 
#                                             true_sigs = kidney_true_sigs, 
#                                             source = "liver_kidney",
#                                             # se = kidney_se,
#                                             # pb_col = "compound.ch1",
#                                             # ridge_benchmark = TRUE,
#                                             BPPARAM = bp)
# kidney_liver_eval_table_n <- sig_eval_table(source_sigs = lapply(kidney_true_sigs, function(x) x$up), 
#                                             pred_sigs = lapply(kidney_true_sigs, function(x) x$up), 
#                                             true_sigs = liver_true_sigs, 
#                                             source = "kidney_liver",
#                                             # se = liver_se,
#                                             # pb_col = "compound.ch1",
#                                             # ridge_benchmark = TRUE,
#                                             BPPARAM = bp)
# all_eval_no_change_table <- rbind(liver_kidney_eval_table_n,
#                                   kidney_liver_eval_table_n)
# saveRDS(all_eval_no_change_table, file.path(SAVE_PATH, "nochange_table.rds"))
all_eval_no_change_table <- readRDS(file.path(SAVE_PATH, "nochange_table.rds"))
# # Ctrl sigs
# liver_pred_sigs <- readRDS(file.path(RECON_PATH, "liver_control.rds"))
# kidney_pred_sigs <- readRDS(file.path(RECON_PATH, "kidney_control.rds"))
# 
# # Eval tables
# liver_kidney_eval_table <- sig_eval_table(source_sigs = lapply(liver_true_sigs, function(x) x$up), 
#                                        pred_sigs = liver_pred_sigs, 
#                                        true_sigs = kidney_true_sigs, 
#                                        source = "liver_kidney",
#                                        # se = kidney_se,
#                                        # pb_col = "compound.ch1",
#                                        # ridge_benchmark = TRUE,
#                                        BPPARAM = bp)
# kidney_liver_eval_table <- sig_eval_table(source_sigs = lapply(kidney_true_sigs, function(x) x$up), 
#                                        pred_sigs = kidney_pred_sigs, 
#                                        true_sigs = liver_true_sigs, 
#                                        source = "kidney_liver",
#                                        # se = liver_se,
#                                        # pb_col = "compound.ch1",
#                                        # ridge_benchmark = TRUE,
#                                        BPPARAM = bp)
# 
# all_eval_table <- rbind(liver_kidney_eval_table,
#                         kidney_liver_eval_table)
# combined_aggregated_df <- merge(all_eval_no_change_table %>% dplyr::select("source", "gene", "NES", "jacc"), 
#                                 all_eval_table %>% mutate(kept_alpha = kept/(kept+displaced)), 
#                                 by = c("source", "gene"), 
#                                 suffixes = c("_FALSE", "_TRUE")) %>% tibble
# saveRDS(combined_aggregated_df, file.path(SAVE_PATH, "ctrl_networkprop_eval_table.rds"))
# 
# combined_aggregated_df <- readRDS(file.path(SAVE_PATH, "ctrl_networkprop_eval_table.rds"))
# 
# 
# # Unpaired 
# combined_aggregated_df %>% 
#   dplyr::group_by(source) %>%
#   dplyr::summarise(NES_mean = mean(NES_TRUE),
#                    jacc_mean = mean(jacc_TRUE)) %>%
#   dplyr::summarize(NES_mean = mean(NES_mean),
#                    jacc_mean = mean(jacc_mean))
# 
# # Plots
# # combined_aggregated_df %>% 
# #   filter(source=="kidney_liver") %>%
# #   ggplot(aes(x=NES_FALSE, y=NES_TRUE)) +
# #   geom_point(aes(alpha=kept_alpha)) +
# #   geom_abline(color="black") +
# #   scale_alpha_continuous("Proportion of Source Kept") +
# #   labs(title = paste("KS Distances in", "kidney Up"), 
# #        x = "KS Distance (Before)",
# #        y = "KS Distance (After)")
# # 
# # combined_aggregated_df %>% 
# #   filter(source=="kidney_liver") %>%
# #   ggplot(aes(x=jacc_FALSE, y=jacc_TRUE)) +
# #   geom_point(aes(alpha=kept_alpha)) +
# #   geom_abline(color="black") +
# #   scale_alpha_continuous("Proportion of Source Kept") +
# #   labs(title = paste("Jaccard Similarity in", "kidney Up"), 
# #        x = "Jaccard (Before)",
# #        y = "Jaccard (After)")
# # 
# combined_aggregated_df %>%
#   filter(source=="liver_kidney") %>%
#   ggplot(aes(x=NES_FALSE, y=NES_TRUE)) +
#   geom_point(aes(alpha=kept_alpha)) +
#   geom_abline(color="black") +
#   scale_alpha_continuous("Proportion of Source Kept") +
#   labs(title = paste("KS Distances in", "liver Up"),
#        x = "KS Distance (Before)",
#        y = "KS Distance (After)")
# 
# combined_aggregated_df %>%
#   filter(source=="liver_kidney") %>%
#   ggplot(aes(x=jacc_FALSE, y=jacc_TRUE)) +
#   geom_point(aes(alpha=kept_alpha)) +
#   geom_abline(color="black") +
#   scale_alpha_continuous("Proportion of Source Kept") +
#   labs(title = paste("Jaccard Similarity in", "liver Up"),
#        x = "Jaccard (Before)",
#        y = "Jaccard (After)")
# 
# 
# # Paired Tests
# combined_aggregated_df %>% 
#   group_by(source) %>% 
#   summarize(Jacc_Wilcox_pval = wilcox.test(jacc_TRUE, jacc_FALSE, paired = TRUE, alternative="greater")$p.value,
#             KS_Wilcox_pval = wilcox.test(NES_TRUE, NES_FALSE, paired = TRUE, alternative="greater")$p.value) %>%
#   summarize(KS_meta_p = fishers_meta_p(KS_Wilcox_pval),
#             Jacc_meta_p = fishers_meta_p(Jacc_Wilcox_pval))

# 10th sigs
kidney_sigs_paths <- Sys.glob(file.path(RECON_PATH, "kidney*1_10th*.rds"))
liver_sigs_paths <- Sys.glob(file.path(RECON_PATH, "liver*1_10th*.rds"))
split_names <- sub(".*(split_[0-9]+).*", "\\1", kidney_sigs_paths)
kidney_sigs <- lapply(kidney_sigs_paths, readRDS)
names(kidney_sigs) <- split_names
split_names <- sub(".*(split_[0-9]+).*", "\\1", liver_sigs_paths)
liver_sigs <- lapply(liver_sigs_paths, readRDS)
names(liver_sigs) <- split_names

# Eval tables
liver_kidney_eval_table <- sig_eval_table(source_sigs = lapply(liver_true_sigs, function(x) x$up), 
                                          pred_sigs = liver_sigs, 
                                          true_sigs = kidney_true_sigs, 
                                          source = "liver_kidney",
                                          splits = TRUE,
                                          split_file = drug_splits_path,
                                          split_pb_col = "drug",
                                          split_type = "10th",
                                          BPPARAM = bp)
kidney_liver_eval_table <- sig_eval_table(source_sigs = lapply(kidney_true_sigs, function(x) x$up), 
                                          pred_sigs = kidney_sigs, 
                                          true_sigs = liver_true_sigs, 
                                          source = "kidney_liver",
                                          splits = TRUE,
                                          split_file = drug_splits_path,
                                          split_pb_col = "drug",
                                          split_type = "10th",
                                          BPPARAM = bp)

all_eval_table <- rbind(liver_kidney_eval_table,
                        kidney_liver_eval_table)
combined_aggregated_df <- merge(all_eval_no_change_table %>% dplyr::select("source", "gene", "NES", "jacc"), 
                                all_eval_table %>% mutate(kept_alpha = kept/(kept+displaced)), 
                                by = c("source", "gene"), 
                                suffixes = c("_FALSE", "_TRUE")) %>% tibble
saveRDS(combined_aggregated_df, file.path(SAVE_PATH, "ctrl_10th_networkprop_eval_table.rds"))

combined_aggregated_df <- readRDS(file.path(SAVE_PATH, "ctrl_10th_networkprop_eval_table.rds"))


# Unpaired 
combined_aggregated_df %>% 
  dplyr::group_by(source) %>%
  dplyr::summarise(NES_mean = mean(NES_TRUE),
                   jacc_mean = mean(jacc_TRUE)) %>%
  dplyr::summarize(NES_mean = mean(NES_mean),
                   jacc_mean = mean(jacc_mean))

# # Plots
combined_aggregated_df %>%
  filter(source=="kidney_liver") %>%
  ggplot(aes(x=NES_FALSE, y=NES_TRUE)) +
  geom_point(aes(alpha=kept_alpha)) +
  geom_abline(color="black") +
  scale_alpha_continuous("Proportion of Source Kept") +
  labs(title = paste("KS Distances in", "kidney Up"),
       x = "KS Distance (Before)",
       y = "KS Distance (After)")

combined_aggregated_df %>%
  filter(source=="kidney_liver") %>%
  ggplot(aes(x=jacc_FALSE, y=jacc_TRUE)) +
  geom_point(aes(alpha=kept_alpha)) +
  geom_abline(color="black") +
  scale_alpha_continuous("Proportion of Source Kept") +
  labs(title = paste("Jaccard Similarity in", "kidney Up"),
       x = "Jaccard (Before)",
       y = "Jaccard (After)")

# combined_aggregated_df %>% 
#   filter(source=="liver_kidney") %>%
#   ggplot(aes(x=NES_FALSE, y=NES_TRUE)) +
#   geom_point(aes(alpha=kept_alpha)) +
#   geom_abline(color="black") +
#   scale_alpha_continuous("Proportion of Source Kept") +
#   labs(title = paste("KS Distances in", "liver Up"), 
#        x = "KS Distance (Before)",
#        y = "KS Distance (After)")
# 
# combined_aggregated_df %>% 
#   filter(source=="liver_kidney") %>%
#   ggplot(aes(x=jacc_FALSE, y=jacc_TRUE)) +
#   geom_point(aes(alpha=kept_alpha)) +
#   geom_abline(color="black") +
#   scale_alpha_continuous("Proportion of Source Kept") +
#   labs(title = paste("Jaccard Similarity in", "liver Up"), 
#        x = "Jaccard (Before)",
#        y = "Jaccard (After)")


# Paired Tests
combined_aggregated_df %>% 
  group_by(source) %>% 
  summarize(Jacc_Wilcox_pval = wilcox.test(jacc_TRUE, jacc_FALSE, paired = TRUE, alternative="greater")$p.value,
            KS_Wilcox_pval = wilcox.test(NES_TRUE, NES_FALSE, paired = TRUE, alternative="greater")$p.value) %>%
  summarize(KS_meta_p = fishers_meta_p(KS_Wilcox_pval),
            Jacc_meta_p = fishers_meta_p(Jacc_Wilcox_pval))


# # 90th sigs
# kidney_sigs_paths <- Sys.glob(file.path(RECON_PATH, "kidney*9_10th*.rds"))
# liver_sigs_paths <- Sys.glob(file.path(RECON_PATH, "liver*9_10th*.rds"))
# split_names <- sub(".*(split_[0-9]+).*", "\\1", kidney_sigs_paths)
# kidney_sigs <- lapply(kidney_sigs_paths, readRDS)
# names(kidney_sigs) <- split_names
# split_names <- sub(".*(split_[0-9]+).*", "\\1", liver_sigs_paths)
# liver_sigs <- lapply(liver_sigs_paths, readRDS)
# names(liver_sigs) <- split_names
# 
# # Eval tables
# liver_kidney_eval_table <- sig_eval_table(source_sigs = lapply(liver_true_sigs, function(x) x$up), 
#                                           pred_sigs = liver_sigs, 
#                                           true_sigs = kidney_true_sigs, 
#                                           source = "liver_kidney",
#                                           splits = TRUE,
#                                           split_file = drug_splits_path,
#                                           split_pb_col = "drug",
#                                           split_type = "90th",
#                                           BPPARAM = bp)
# kidney_liver_eval_table <- sig_eval_table(source_sigs = lapply(kidney_true_sigs, function(x) x$up), 
#                                           pred_sigs = kidney_sigs, 
#                                           true_sigs = liver_true_sigs, 
#                                           source = "kidney_liver",
#                                           splits = TRUE,
#                                           split_file = drug_splits_path,
#                                           split_pb_col = "drug",
#                                           split_type = "90th",
#                                           BPPARAM = bp)
# 
# all_eval_table <- rbind(liver_kidney_eval_table,
#                         kidney_liver_eval_table)
# combined_aggregated_df <- merge(all_eval_no_change_table %>% dplyr::select("source", "gene", "NES", "jacc"), 
#                                 all_eval_table %>% mutate(kept_alpha = kept/(kept+displaced)), 
#                                 by = c("source", "gene"), 
#                                 suffixes = c("_FALSE", "_TRUE")) %>% tibble
# saveRDS(combined_aggregated_df, file.path(SAVE_PATH, "ctrl_90th_networkprop_eval_table.rds"))
# 
# combined_aggregated_df <- readRDS(file.path(SAVE_PATH, "ctrl_90th_networkprop_eval_table.rds"))
# 
# # Unpaired 
# combined_aggregated_df %>% 
#   dplyr::group_by(source) %>%
#   dplyr::summarise(NES_mean = mean(NES_TRUE),
#                    jacc_mean = mean(jacc_TRUE)) %>%
#   dplyr::summarize(NES_mean = mean(NES_mean),
#                    jacc_mean = mean(jacc_mean))
# 
# # # Plots
# combined_aggregated_df %>%
#   filter(source=="kidney_liver") %>%
#   ggplot(aes(x=NES_FALSE, y=NES_TRUE)) +
#   geom_point(aes(alpha=kept_alpha)) +
#   geom_abline(color="black") +
#   scale_alpha_continuous("Proportion of Source Kept") +
#   labs(title = paste("KS Distances in", "kidney Up"),
#        x = "KS Distance (Before)",
#        y = "KS Distance (After)")
# 
# combined_aggregated_df %>%
#   filter(source=="kidney_liver") %>%
#   ggplot(aes(x=jacc_FALSE, y=jacc_TRUE)) +
#   geom_point(aes(alpha=kept_alpha)) +
#   geom_abline(color="black") +
#   scale_alpha_continuous("Proportion of Source Kept") +
#   labs(title = paste("Jaccard Similarity in", "kidney Up"),
#        x = "Jaccard (Before)",
#        y = "Jaccard (After)")
# # 
# # combined_aggregated_df %>% 
# #   filter(source=="liver_kidney") %>%
# #   ggplot(aes(x=NES_FALSE, y=NES_TRUE)) +
# #   geom_point(aes(alpha=kept_alpha)) +
# #   geom_abline(color="black") +
# #   scale_alpha_continuous("Proportion of Source Kept") +
# #   labs(title = paste("KS Distances in", "liver Up"), 
# #        x = "KS Distance (Before)",
# #        y = "KS Distance (After)")
# # 
# # combined_aggregated_df %>% 
# #   filter(source=="liver_kidney") %>%
# #   ggplot(aes(x=jacc_FALSE, y=jacc_TRUE)) +
# #   geom_point(aes(alpha=kept_alpha)) +
# #   geom_abline(color="black") +
# #   scale_alpha_continuous("Proportion of Source Kept") +
# #   labs(title = paste("Jaccard Similarity in", "liver Up"), 
# #        x = "Jaccard (Before)",
# #        y = "Jaccard (After)")
# 
# 
# # Paired Tests
# combined_aggregated_df %>% 
#   group_by(source) %>% 
#   summarize(Jacc_Wilcox_pval = wilcox.test(jacc_TRUE, jacc_FALSE, paired = TRUE, alternative="greater")$p.value,
#             KS_Wilcox_pval = wilcox.test(NES_TRUE, NES_FALSE, paired = TRUE, alternative="greater")$p.value) %>%
#   summarize(KS_meta_p = fishers_meta_p(KS_Wilcox_pval),
#             Jacc_meta_p = fishers_meta_p(Jacc_Wilcox_pval))
# 
# 
