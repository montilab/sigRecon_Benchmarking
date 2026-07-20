library(tidyverse)
library(sigrecon)
library(BiocParallel)
library(SummarizedExperiment)
bp <- make_bpparam(workers = 15, RNGseed = 123, type="multicore")
register(bp)

PATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/projects/sig_recon")
SE_PATH <- file.path(Sys.getenv("AGED"), "CBMrepositoryData/perturbational_data/perturb-seq/")
RECON_PATH <- file.path(PATH, "data/sigs/perturb-seq/networkprop/")
SAVE_PATH <- file.path(PATH, "results/eval/perturb-seq")

# Loading sigs
k562_true_sigs <- perturbseq.k562
rpe1_true_sigs <- perturbseq.rpe1

# Loading drug splits
drug_splits_path <- file.path(PATH, "data/sigs/perturb-seq/pb_splits.csv")

# No Change
all_eval_no_change_table <- readRDS(file.path(SAVE_PATH, "no_change_eval_table.rds"))
all_eval_no_change_table$source <- if_else(all_eval_no_change_table$source == "k562",
                                           "k562_rpe1",
                                           "rpe1_k562")
# paired_eval_table() takes a file path (not a data frame) for the
# no-change baseline, so the relabeled table needs to be written out once
# before it can be used with the merges below.
no_change_path <- tempfile(fileext = ".rds")
saveRDS(all_eval_no_change_table, no_change_path)
# # Ctrl sigs
# k562_pred_sigs <- readRDS(file.path(RECON_PATH, "k562_control.rds"))
# rpe1_pred_sigs <- readRDS(file.path(RECON_PATH, "rpe1_control.rds"))
#  
# # Eval tables
# k562_rpe1_eval_table <- sig_eval_table(source_sigs = lapply(k562_true_sigs, function(x) x$up),
#                                        pred_sigs = k562_pred_sigs,
#                                        true_sigs = rpe1_true_sigs,
#                                        source = "k562_rpe1",
#                                        BPPARAM = bp)
# rpe1_k562_eval_table <- sig_eval_table(source_sigs = lapply(rpe1_true_sigs, function(x) x$up),
#                                        pred_sigs = rpe1_pred_sigs,
#                                        true_sigs = k562_true_sigs,
#                                        source = "rpe1_k562",
#                                        BPPARAM = bp)
# 
# all_eval_table <- rbind(k562_rpe1_eval_table,
#                         rpe1_k562_eval_table)
# combined_aggregated_df <- paired_eval_table(all_eval_table, no_change_path)
# saveRDS(combined_aggregated_df, file.path(SAVE_PATH, "ctrl_networkprop_eval_table.rds"))
#
#
# # 10th sigs
# rpe1_sigs_paths <- Sys.glob(file.path(RECON_PATH, "rpe1*1_10th*.rds"))
# k562_sigs_paths <- Sys.glob(file.path(RECON_PATH, "k562*1_10th*.rds"))
# split_names <- sub(".*(split_[0-9]+).*", "\\1", rpe1_sigs_paths)
# rpe1_sigs <- lapply(rpe1_sigs_paths, readRDS)
# names(rpe1_sigs) <- split_names
# split_names <- sub(".*(split_[0-9]+).*", "\\1", k562_sigs_paths)
# k562_sigs <- lapply(k562_sigs_paths, readRDS)
# names(k562_sigs) <- split_names
# 
# # Eval tables
# k562_rpe1_eval_table <- sig_eval_table(source_sigs = lapply(k562_true_sigs, function(x) x$up), 
#                                           pred_sigs = k562_sigs, 
#                                           true_sigs = rpe1_true_sigs, 
#                                           source = "k562_rpe1",
#                                           splits = TRUE,
#                                           split_file = drug_splits_path,
#                                           split_pb_col = "pb",
#                                           split_type = "10th",
#                                           BPPARAM = bp)
# rpe1_k562_eval_table <- sig_eval_table(source_sigs = lapply(rpe1_true_sigs, function(x) x$up), 
#                                           pred_sigs = rpe1_sigs, 
#                                           true_sigs = k562_true_sigs, 
#                                           source = "rpe1_k562",
#                                           splits = TRUE,
#                                           split_file = drug_splits_path,
#                                           split_pb_col = "pb",
#                                           split_type = "10th",
#                                           BPPARAM = bp)
# 
# all_eval_table <- rbind(k562_rpe1_eval_table,
#                         rpe1_k562_eval_table)
# combined_aggregated_df <- paired_eval_table(all_eval_table, no_change_path)
# saveRDS(combined_aggregated_df, file.path(SAVE_PATH, "ctrl_10th_networkprop_eval_table.rds"))


# 90th sigs
rpe1_sigs_paths <- Sys.glob(file.path(RECON_PATH, "rpe1*9_10th*.rds"))
k562_sigs_paths <- Sys.glob(file.path(RECON_PATH, "k562*9_10th*.rds"))
split_names <- sub(".*(split_[0-9]+).*", "\\1", rpe1_sigs_paths)
rpe1_sigs <- lapply(rpe1_sigs_paths, readRDS)
names(rpe1_sigs) <- split_names
split_names <- sub(".*(split_[0-9]+).*", "\\1", k562_sigs_paths)
k562_sigs <- lapply(k562_sigs_paths, readRDS)
names(k562_sigs) <- split_names

# Eval tables
k562_rpe1_eval_table <- sig_eval_table(source_sigs = lapply(k562_true_sigs, function(x) x$up),
                                          pred_sigs = k562_sigs,
                                          true_sigs = rpe1_true_sigs,
                                          source = "k562_rpe1",
                                          splits = TRUE,
                                          split_file = drug_splits_path,
                                          split_pb_col = "pb",
                                          split_type = "90th",
                                          BPPARAM = bp)
rpe1_k562_eval_table <- sig_eval_table(source_sigs = lapply(rpe1_true_sigs, function(x) x$up),
                                          pred_sigs = rpe1_sigs,
                                          true_sigs = k562_true_sigs,
                                          source = "rpe1_k562",
                                          splits = TRUE,
                                          split_file = drug_splits_path,
                                          split_pb_col = "pb",
                                          split_type = "90th",
                                          BPPARAM = bp)

all_eval_table <- rbind(k562_rpe1_eval_table,
                        rpe1_k562_eval_table)
combined_aggregated_df <- paired_eval_table(all_eval_table, no_change_path)
saveRDS(combined_aggregated_df, file.path(SAVE_PATH, "ctrl_90th_networkprop_eval_table.rds"))

# combined_aggregated_df <- readRDS(file.path(SAVE_PATH, "ctrl_networkprop_eval_table.rds"))
# combined_aggregated_df <- readRDS(file.path(SAVE_PATH, "ctrl_10th_networkprop_eval_table.rds"))
# combined_aggregated_df <- readRDS(file.path(SAVE_PATH, "ctrl_90th_networkprop_eval_table.rds"))
# 
# # Unpaired
# combined_aggregated_df %>%
#   dplyr::group_by(source) %>%
#   dplyr::summarise(NES_mean = mean(NES_TRUE, na.rm = TRUE),
#                    jacc_mean = mean(jacc_TRUE, na.rm = TRUE)) %>%
#   dplyr::summarize(NES_mean = mean(NES_mean),
#                    jacc_mean = mean(jacc_mean))
# combined_aggregated_df %>%
#   group_by(source) %>%
#   summarize(Jacc_Wilcox_pval = wilcox.test(jacc_TRUE, jacc_FALSE, paired = TRUE, alternative="greater")$p.value,
#             KS_Wilcox_pval = wilcox.test(NES_TRUE, NES_FALSE, paired = TRUE, alternative="greater")$p.value) %>%
#   summarize(KS_meta_p = fishers_meta_p(KS_Wilcox_pval),
#             Jacc_meta_p = fishers_meta_p(Jacc_Wilcox_pval))


