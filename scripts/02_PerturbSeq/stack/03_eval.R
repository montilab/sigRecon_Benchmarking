library(tidyverse)
library(sigrecon)
library(BiocParallel)
bp <- make_bpparam(workers = 10, RNGseed = 123, type="multicore")
BiocParallel::register(bp)
PATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/projects/sig_recon")
DATA_PATH <- file.path(Sys.getenv("AGED"), "CBMrepositoryData/perturbational_data/replogle_2022/")
RECON_PATH <- file.path(PATH, "data/sigs/perturb-seq/stack/")
SAVE_PATH <- file.path(PATH, "results/eval/perturb-seq")

# Load signatures
rpe1_sig_up <- lapply(perturbseq.rpe1, function(x) x$up)
k562_sig_up <- lapply(perturbseq.k562, function(x) x$up)

# Loading drug splits
pb_splits_path <- file.path(PATH, "data/sigs/perturb-seq/pb_splits.csv")
pb_splits <- read.csv(pb_splits_path)
all_eval_no_change_table <- readRDS(file.path(SAVE_PATH, "no_change_eval_table.rds"))
all_eval_no_change_table$source <- if_else(all_eval_no_change_table$source == "k562",
                                           "k562_rpe1",
                                           "rpe1_k562")

# # 1. Control experiment
# rpe1_ctrl_sigs_df <- readRDS(file.path(RECON_PATH, "rpe1_ctrl_sigs_df.rds"))
# k562_ctrl_sigs_df <- readRDS(file.path(RECON_PATH, "k562_ctrl_sigs_df.rds"))
# rpe1_ctrl_sigs <- sig_filter_fn(rpe1_ctrl_sigs_df,
#                                 perts = pb_splits$pb,
#                                 pval_col = "pvalue",
#                                 alpha = 1.1,
#                                 log2fc_col = "log2FoldChange",
#                                 geneid_col = "gene")
# k562_ctrl_sigs <- sig_filter_fn(k562_ctrl_sigs_df,
#                                 perts = pb_splits$pb,
#                                 pval_col = "pvalue",
#                                 alpha = 1.1,
#                                 log2fc_col = "log2FoldChange",
#                                 geneid_col = "gene")
# rpe1_ctrl_sigs <- lapply(rpe1_ctrl_sigs, function(x) x$up)
# k562_ctrl_sigs <- lapply(k562_ctrl_sigs, function(x) x$up)
# 
# k562_rpe1_eval_table <- sig_eval_table(source_sigs = k562_sig_up,
#                                           pred_sigs = k562_ctrl_sigs,
#                                           true_sigs = perturbseq.rpe1,
#                                           source = "k562_rpe1",
#                                           BPPARAM = bp)
# rpe1_k562_eval_table <- sig_eval_table(source_sigs = rpe1_sig_up,
#                                           pred_sigs = rpe1_ctrl_sigs,
#                                           true_sigs = perturbseq.k562,
#                                           source = "rpe1_k562",
#                                           BPPARAM = bp)
# all_eval_table <- rbind(k562_rpe1_eval_table,
#                         rpe1_k562_eval_table)
# combined_aggregated_df <- merge(all_eval_no_change_table %>% dplyr::select("source", "gene", "NES", "jacc"),
#                                 all_eval_table %>% mutate(kept_alpha = kept/(kept+displaced)),
#                                 by = c("source", "gene"),
#                                 suffixes = c("_FALSE", "_TRUE")) %>% tibble
# 
# saveRDS(combined_aggregated_df, file.path(SAVE_PATH, "ctrl_stack_eval_table.rds"))

# ------------------------------------------------
# 3. 1/10th experiment
# ------------------------------------------------
do_save <- TRUE
# None of the genes reach significance after stack generation which is why alpha = 1
if(do_save) {
  rpe1_paths <- Sys.glob(file.path(RECON_PATH,
                                   paste0("rpe1_10th_split_*.rds")))
  k562_paths  <- Sys.glob(file.path(RECON_PATH,
                                    paste0("k562_10th_split_*.rds")))

  split_names <- sub(".*(split_[0-9]+).*", "\\1", rpe1_paths)
  rpe1_sigs_df <- lapply(rpe1_paths, readRDS)
  rpe1_sigs <- lapply(rpe1_sigs_df, function(x) sig_filter_fn(x,
                                                              perts = unique(pb_splits$pb),
                                                              log2fc_col = "log2FoldChange",
                                                              pval_col = "pvalue",
                                                              geneid_col = "gene",
                                                              alpha = 1.1))
  rpe1_sigs <- lapply(rpe1_sigs, function(x) lapply(x, function(y) y$up))
  names(rpe1_sigs) <- split_names
  saveRDS(rpe1_sigs, file.path(RECON_PATH, "rpe1_10th_sigs.rds"))

  split_names <- sub(".*(split_[0-9]+).*", "\\1", k562_paths)
  k562_sigs_df <- lapply(k562_paths, readRDS)
  k562_sigs <- lapply(k562_sigs_df, function(x) sig_filter_fn(x,
                                                              perts = unique(pb_splits$pb),
                                                              log2fc_col = "log2FoldChange",
                                                              pval_col = "pvalue",
                                                              geneid_col = "gene",
                                                              alpha = 1.1))
  k562_sigs <- lapply(k562_sigs, function(x) lapply(x, function(y) y$up))
  names(k562_sigs) <- split_names
  saveRDS(k562_sigs, file.path(RECON_PATH, "k562_10th_sigs.rds"))
} else {
  rpe1_sigs <- readRDS(file.path(RECON_PATH, "rpe1_10th_sigs.rds"))
  k562_sigs <- readRDS(file.path(RECON_PATH, "k562_10th_sigs.rds"))
}

k562_rpe1_eval_table <- sig_eval_table(
  source_sigs = k562_sig_up,
  pred_sigs = k562_sigs,
  true_sigs = perturbseq.rpe1,
  source = "k562_rpe1",
  splits = TRUE,
  split_file = pb_splits_path,
  split_pb_col = "pb",
  split_type = "10th",
  BPPARAM = bp
)

rpe1_k562_eval_table <- sig_eval_table(
  source_sigs = rpe1_sig_up,
  pred_sigs = rpe1_sigs,
  true_sigs = perturbseq.k562,
  source = "rpe1_k562",
  splits = TRUE,
  split_file = pb_splits_path,
  split_pb_col = "pb",
  split_type = "10th",
  BPPARAM = bp
)

all_eval_table <- rbind(k562_rpe1_eval_table,
                        rpe1_k562_eval_table)

combined_aggregated_df <- merge(
  all_eval_no_change_table %>% dplyr::select("source","gene","NES","jacc"),
  all_eval_table %>% mutate(kept_alpha = kept/(kept+displaced)),
  by=c("source","gene"),
  suffixes=c("_FALSE","_TRUE")
) %>% tibble()

saveRDS(combined_aggregated_df,
        file.path(SAVE_PATH,
                  paste0("10th_stack_eval_table.rds")))


# # ------------------------------------------------
# # 3. 9/10th experiment
# # ------------------------------------------------
# do_save <- TRUE
# # None of the genes reach significance after stack generation which is why alpha = 1
# if(do_save) {
#   rpe1_paths <- Sys.glob(file.path(RECON_PATH,
#                                    paste0("rpe1_90th_split_*.rds")))
#   k562_paths  <- Sys.glob(file.path(RECON_PATH,
#                                     paste0("k562_90th_split_*.rds")))
#   
#   split_names <- sub(".*(split_[0-9]+).*", "\\1", rpe1_paths)
#   rpe1_sigs_df <- lapply(rpe1_paths, readRDS)
#   rpe1_sigs <- lapply(rpe1_sigs_df, function(x) sig_filter_fn(x, 
#                                                               perts = unique(pb_splits$pb),
#                                                               log2fc_col = "log2FoldChange",
#                                                               pval_col = "pvalue",
#                                                               geneid_col = "gene",
#                                                               alpha = 1.1))
#   rpe1_sigs <- lapply(rpe1_sigs, function(x) lapply(x, function(y) y$up))
#   names(rpe1_sigs) <- split_names
#   saveRDS(rpe1_sigs, file.path(RECON_PATH, "rpe1_90th_sigs.rds"))
#   
#   split_names <- sub(".*(split_[0-9]+).*", "\\1", k562_paths)
#   k562_sigs_df <- lapply(k562_paths, readRDS)
#   k562_sigs <- lapply(k562_sigs_df, function(x) sig_filter_fn(x, 
#                                                               perts = unique(pb_splits$pb),
#                                                               log2fc_col = "log2FoldChange",
#                                                               pval_col = "pvalue",
#                                                               geneid_col = "gene",
#                                                               alpha = 1.1))
#   k562_sigs <- lapply(k562_sigs, function(x) lapply(x, function(y) y$up))
#   names(k562_sigs) <- split_names
#   saveRDS(k562_sigs, file.path(RECON_PATH, "k562_90th_sigs.rds"))
# } else {
#   rpe1_sigs <- readRDS(file.path(RECON_PATH, "rpe1_90th_sigs.rds"))
#   k562_sigs <- readRDS(file.path(RECON_PATH, "k562_90th_sigs.rds"))
# }
# 
# k562_rpe1_eval_table <- sig_eval_table(
#   source_sigs = k562_sig_up,
#   pred_sigs = k562_sigs,
#   true_sigs = perturbseq.rpe1,
#   source = "k562_rpe1",
#   splits = TRUE,
#   split_file = pb_splits_path,
#   split_pb_col = "pb",
#   split_type = "90th",
#   BPPARAM = bp
# )
# 
# rpe1_k562_eval_table <- sig_eval_table(
#   source_sigs = rpe1_sig_up,
#   pred_sigs = rpe1_sigs,
#   true_sigs = perturbseq.k562,
#   source = "rpe1_k562",
#   splits = TRUE,
#   split_file = pb_splits_path,
#   split_pb_col = "pb",
#   split_type = "90th",
#   BPPARAM = bp
# )
# 
# all_eval_table <- rbind(k562_rpe1_eval_table,
#                         rpe1_k562_eval_table)
# 
# combined_aggregated_df <- merge(
#   all_eval_no_change_table %>% dplyr::select("source","gene","NES","jacc"),
#   all_eval_table %>% mutate(kept_alpha = kept/(kept+displaced)),
#   by=c("source","gene"),
#   suffixes=c("_FALSE","_TRUE")
# ) %>% tibble()
# 
# saveRDS(combined_aggregated_df,
#         file.path(SAVE_PATH,
#                   paste0("90th_stack_eval_table.rds")))

