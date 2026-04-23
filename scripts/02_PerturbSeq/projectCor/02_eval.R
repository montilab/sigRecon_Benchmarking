library(tidyverse)
library(sigrecon)
library(BiocParallel)
bp <- make_bpparam(workers = 10, RNGseed = 123, type="multicore")
BiocParallel::register(bp)
PATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/projects/sig_recon")
DATA_PATH <- file.path(Sys.getenv("AGED"), "CBMrepositoryData/perturbational_data/replogle_2022/")
RECON_PATH <- file.path(PATH, "data/sigs/perturb-seq/projectCor/")
SAVE_PATH <- file.path(PATH, "results/eval/perturb-seq")

# Load signatures
rpe1_sig_up <- lapply(perturbseq.rpe1, function(x) x$up)
k562_sig_up <- lapply(perturbseq.k562, function(x) x$up)

# Loading drug splits
drug_splits_path <- file.path(PATH, "data/sigs/perturb-seq/pb_splits.csv")

all_eval_no_change_table <- readRDS(file.path(SAVE_PATH, "no_change_eval_table.rds"))
all_eval_no_change_table$source <- if_else(all_eval_no_change_table$source == "k562",
                                           "k562_rpe1",
                                           "rpe1_k562")
# # 1. Control experiment
# k562_eigen_pred_sigs <- readRDS(file.path(RECON_PATH, "k562_rpe1_ctrl_eigen.rds"))
# k562_gsva_pred_sigs <- readRDS(file.path(RECON_PATH, "k562_rpe1_ctrl_gsva.rds"))
# rpe1_eigen_pred_sigs <- readRDS(file.path(RECON_PATH, "rpe1_k562_ctrl_eigen.rds"))
# rpe1_gsva_pred_sigs <- readRDS(file.path(RECON_PATH, "rpe1_k562_ctrl_gsva.rds"))
# 
# k562_rpe1_eval_table <- sig_eval_table(source_sigs = k562_sig_up,
#                                           pred_sigs = k562_eigen_pred_sigs,
#                                           true_sigs = perturbseq.rpe1,
#                                           source = "k562_rpe1",
#                                           BPPARAM = bp)
# rpe1_k562_eval_table <- sig_eval_table(source_sigs = rpe1_sig_up,
#                                           pred_sigs = rpe1_eigen_pred_sigs,
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
# saveRDS(combined_aggregated_df, file.path(SAVE_PATH, "ctrl_projectcor_eigen_eval_table.rds"))
# 
# k562_rpe1_eval_table <- sig_eval_table(source_sigs = k562_sig_up,
#                                           pred_sigs = k562_gsva_pred_sigs,
#                                           true_sigs = perturbseq.rpe1,
#                                           source = "k562_rpe1",
#                                           BPPARAM = bp)
# rpe1_k562_eval_table <- sig_eval_table(source_sigs = rpe1_sig_up,
#                                           pred_sigs = rpe1_gsva_pred_sigs,
#                                           true_sigs = perturbseq.k562,
#                                           source = "rpe1_k562",
#                                           BPPARAM = bp)
# all_eval_table <- rbind(k562_rpe1_eval_table,
#                         rpe1_k562_eval_table)
# combined_aggregated_df <- merge(all_eval_no_change_table %>% dplyr::select("source", "gene", "NES", "jacc"),
#                                 all_eval_table %>% mutate(kept_alpha = kept/(kept+displaced)),
#                                 by = c("source", "gene"),
#                                 suffixes = c("_FALSE", "_TRUE")) %>% tibble
# saveRDS(combined_aggregated_df, file.path(SAVE_PATH, "ctrl_projectcor_gsva_eval_table.rds"))
# 
# # ------------------------------------------------
# # 2. 1/10th experiment
# # ------------------------------------------------
# 
# for(regime in c("gsva","eigen")){
#   
#   message("Processing 1/10th ", regime)
#   
#   rpe1_paths <- Sys.glob(file.path(RECON_PATH,
#                                      paste0("rpe1_k562_10th_", regime, "_split_*.rds")))
#   k562_paths  <- Sys.glob(file.path(RECON_PATH,
#                                      paste0("k562_rpe1_10th_", regime, "_split_*.rds")))
#   
#   split_names <- sub(".*(split_[0-9]+).*", "\\1", rpe1_paths)
#   rpe1_sigs <- lapply(rpe1_paths, readRDS)
#   names(rpe1_sigs) <- split_names
#   
#   split_names <- sub(".*(split_[0-9]+).*", "\\1", k562_paths)
#   k562_sigs <- lapply(k562_paths, readRDS)
#   names(k562_sigs) <- split_names
#   
#   k562_rpe1_eval_table <- sig_eval_table(
#     source_sigs = k562_sig_up,
#     pred_sigs = k562_sigs,
#     true_sigs = perturbseq.rpe1,
#     source = "k562_rpe1",
#     splits = TRUE,
#     split_file = drug_splits_path,
#     split_pb_col = "pb",
#     split_type = "10th",
#     BPPARAM = bp
#   )
#   
#   rpe1_k562_eval_table <- sig_eval_table(
#     source_sigs = rpe1_sig_up,
#     pred_sigs = rpe1_sigs,
#     true_sigs = perturbseq.k562,
#     source = "rpe1_k562",
#     splits = TRUE,
#     split_file = drug_splits_path,
#     split_pb_col = "pb",
#     split_type = "10th",
#     BPPARAM = bp
#   )
#   
#   all_eval_table <- rbind(k562_rpe1_eval_table,
#                           rpe1_k562_eval_table)
#   
#   combined_aggregated_df <- merge(
#     all_eval_no_change_table %>% dplyr::select("source","gene","NES","jacc"),
#     all_eval_table %>% mutate(kept_alpha = kept/(kept+displaced)),
#     by=c("source","gene"),
#     suffixes=c("_FALSE","_TRUE")
#   ) %>% tibble()
#   
#   saveRDS(combined_aggregated_df,
#           file.path(SAVE_PATH,
#                     paste0("10th_projectcor_",regime,"_eval_table.rds")))
# }
# 

# ------------------------------------------------
# 3. 9/10th experiment
# ------------------------------------------------

for(regime in c("gsva","eigen")){

  message("Processing 9/10th ", regime)

  rpe1_paths <- Sys.glob(file.path(RECON_PATH,
                                     paste0("rpe1_k562_90th_", regime, "_split_*.rds")))
  k562_paths  <- Sys.glob(file.path(RECON_PATH,
                                     paste0("k562_rpe1_90th_", regime, "_split_*.rds")))

  split_names <- sub(".*(split_[0-9]+).*", "\\1", rpe1_paths)
  rpe1_sigs <- lapply(rpe1_paths, readRDS)
  names(rpe1_sigs) <- split_names

  split_names <- sub(".*(split_[0-9]+).*", "\\1", k562_paths)
  k562_sigs <- lapply(k562_paths, readRDS)
  names(k562_sigs) <- split_names

  k562_rpe1_eval_table <- sig_eval_table(
    source_sigs = k562_sig_up,
    pred_sigs = k562_sigs,
    true_sigs = perturbseq.rpe1,
    source = "k562_rpe1",
    splits = TRUE,
    split_file = drug_splits_path,
    split_pb_col = "pb",
    split_type = "90th",
    BPPARAM = bp
  )

  rpe1_k562_eval_table <- sig_eval_table(
    source_sigs = rpe1_sig_up,
    pred_sigs = rpe1_sigs,
    true_sigs = perturbseq.k562,
    source = "rpe1_k562",
    splits = TRUE,
    split_file = drug_splits_path,
    split_pb_col = "pb",
    split_type = "90th",
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
                    paste0("90th_projectcor_",regime,"_eval_table.rds")))
}
# 
combined_aggregated_df <- readRDS(file.path(SAVE_PATH, "ctrl_projectcor_gsva_eval_table.rds"))
combined_aggregated_df <- readRDS(file.path(SAVE_PATH, "ctrl_projectcor_eigen_eval_table.rds"))
combined_aggregated_df <- readRDS(file.path(SAVE_PATH, "10th_projectcor_gsva_eval_table.rds"))
combined_aggregated_df <- readRDS(file.path(SAVE_PATH, "10th_projectcor_eigen_eval_table.rds"))
combined_aggregated_df <- readRDS(file.path(SAVE_PATH, "90th_projectcor_gsva_eval_table.rds"))
combined_aggregated_df <- readRDS(file.path(SAVE_PATH, "90th_projectcor_eigen_eval_table.rds"))
combined_aggregated_df %>%
  dplyr::group_by(source) %>%
  dplyr::summarise(NES_mean = mean(NES_TRUE),
                   jacc_mean = mean(jacc_TRUE)) %>%
  dplyr::summarize(NES_mean = mean(NES_mean),
                   jacc_mean = mean(jacc_mean))
combined_aggregated_df %>%
  group_by(source) %>%
  summarize(Jacc_Wilcox_pval = wilcox.test(jacc_TRUE, jacc_FALSE, paired = TRUE, alternative="greater")$p.value,
            KS_Wilcox_pval = wilcox.test(NES_TRUE, NES_FALSE, paired = TRUE, alternative="greater")$p.value) %>%
  summarize(KS_meta_p = fishers_meta_p(KS_Wilcox_pval),
            Jacc_meta_p = fishers_meta_p(Jacc_Wilcox_pval))
