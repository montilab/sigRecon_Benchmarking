## ----echo=FALSE, message=FALSE----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(tidyverse)
library(sigrecon)
library(SummarizedExperiment)
library(BiocParallel)
PROJECT_PATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/projects/sig_recon/")
DATAPATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/projects/sig_recon/data/sigs/perturb-seq")
do_save <- FALSE
bp <- make_bpparam(workers = 15, RNGseed = 123, progress = TRUE, type = "multicore")
register(bp)


## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#True Sig
rpe1_true_sig <- sigrecon::perturbseq.rpe1
rpe1_true_sig_up <- lapply(rpe1_true_sig, function(x) x$up)
k562_true_sig <- sigrecon::perturbseq.k562
k562_true_sig_up <- lapply(k562_true_sig, function(x) x$up)


## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
if(do_save) {
  k562_recon_eval_table <- sig_eval_table(source_sigs = k562_true_sig_up,
                                         pred_sigs = k562_true_sig_up, 
                                         true_sigs = rpe1_true_sig,
                                         source = "k562")
  
  rpe1_recon_eval_table <- sig_eval_table(source_sigs = rpe1_true_sig_up,
                                         pred_sigs = rpe1_true_sig_up, 
                                         true_sigs = k562_true_sig,
                                         source = "rpe1")
  
  no_change_df <- rbind(k562_recon_eval_table, rpe1_recon_eval_table)

  saveRDS(no_change_df, file.path(PROJECT_PATH, "results/eval/perturb-seq/no_change_eval_table.rds"))
} else {
  no_change_df <- readRDS(file.path(PROJECT_PATH, "results/eval/perturb-seq/no_change_eval_table.rds"))
}

no_change_path <- file.path(PROJECT_PATH, "results/eval/perturb-seq/no_change_eval_table.rds")


## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
no_change_df %>% 
  dplyr::group_by(source) %>%
  summarize(NES_mean = mean(NES, na.rm = TRUE),
            jacc = mean(jacc))


## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
no_change_df %>% 
  dplyr::group_by(source) %>%
  summarize(NES = mean(NES, na.rm = TRUE),
            jacc = mean(jacc)) %>%
  summarize(NES_mean = mean(NES, na.rm = TRUE),
            jacc_mean = mean(jacc))


## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Pred Sig
k562_sig <- readRDS(file.path(DATAPATH, "mean/ctrl/k562_ctrl.rds"))[["k562"]]
rpe1_sig <- readRDS(file.path(DATAPATH, "mean/ctrl/rpe1_ctrl.rds"))[["rpe1"]]
rpe1_pred_sig <- rep(list(rpe1_sig[["up"]]), length(rpe1_true_sig))
names(rpe1_pred_sig) <- names(rpe1_true_sig)
k562_pred_sig <- rep(list(k562_sig[["up"]]), length(k562_true_sig))
names(k562_pred_sig) <- names(k562_true_sig)


## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
if(do_save) {
  k562_recon_eval_table <- sig_eval_table(source_sigs = k562_true_sig_up,
                                           pred_sigs = rpe1_pred_sig, 
                                           true_sigs = rpe1_true_sig,
                                           source = "k562",
                                           target = "rpe1")
  rpe1_recon_eval_table <- sig_eval_table(source_sigs = rpe1_true_sig_up,
                                            pred_sigs = k562_pred_sig, 
                                            true_sigs = k562_true_sig,
                                            source = "rpe1",
                                            target = "k562")

  change_df <- rbind(k562_recon_eval_table,
                     rpe1_recon_eval_table)
  saveRDS(change_df, file.path(PROJECT_PATH, "results/eval/perturb-seq/ctrl_mean_eval_table.rds"))
} else {
  change_df <- readRDS(file.path(PROJECT_PATH, "results/eval/perturb-seq/ctrl_mean_eval_table.rds"))
}



## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
change_df %>% 
  dplyr::group_by(source) %>%
  summarize(NES_mean = mean(NES),
            jacc = mean(jacc))


## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
change_df %>% 
  dplyr::group_by(source) %>%
  summarize(NES = mean(NES),
            jacc = mean(jacc)) %>%
  summarize(NES_mean = mean(NES),
            jacc_mean = mean(jacc))


## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
combined_aggregated_df <- paired_eval_table(change_df, no_change_path)


## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
combined_aggregated_df %>% 
  group_by(source) %>% 
  summarize(Jacc_Wilcox_pval = wilcox.test(jacc_TRUE, jacc_FALSE, paired = TRUE, alternative="greater")$p.value,
            NES_Wilcox_pval = wilcox.test(NES_TRUE, NES_FALSE, paired = TRUE, alternative="greater")$p.value) %>%
  summarize(NES_meta_p = fishers_meta_p(NES_Wilcox_pval),
            Jacc_meta_p = fishers_meta_p(Jacc_Wilcox_pval))


## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Pred Sig
k562_sig <- readRDS(file.path(DATAPATH, "mean/10th_perturb/k562_10_signatures.rds"))
rpe1_sig <- readRDS(file.path(DATAPATH, "mean/10th_perturb/rpe1_10_signatures.rds"))
rpe1_pred_sig <- list()
for(split in 1:length(rpe1_sig)) {
  rpe1_pred_sig[[split]] <- rep(list(rpe1_sig[[split]][["up"]]), length(rpe1_true_sig))
  names(rpe1_pred_sig[[split]]) <- names(rpe1_true_sig)  
}
names(rpe1_pred_sig) <- paste0("split_", 1:10)
k562_pred_sig <- list()
for(split in 1:length(rpe1_sig)) {
  k562_pred_sig[[split]] <- rep(list(k562_sig[[split]][["up"]]), length(k562_true_sig))
  names(k562_pred_sig[[split]]) <- names(k562_true_sig)  
}
names(k562_pred_sig) <- paste0("split_", 1:10)


## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
do_save <- TRUE
if(do_save) {
  k562_recon_eval_table <- sig_eval_table(source_sigs = k562_true_sig_up,
                                          pred_sigs = rpe1_pred_sig, 
                                          true_sigs = rpe1_true_sig,
                                          source = "k562",
                                          target = "rpe1",
                                          splits = TRUE,
                                          split_file = file.path(PROJECT_PATH, "data/sigs/perturb-seq/pb_splits.csv"),
                                          split_pb_col = "pb",
                                          split_type = "10th",
                                          BPPARAM = bp)
  rpe1_recon_eval_table <- sig_eval_table(source_sigs = rpe1_true_sig_up,
                                          pred_sigs = k562_pred_sig, 
                                          true_sigs = k562_true_sig,
                                          source = "rpe1",
                                          target = "k562",
                                          splits = TRUE,
                                          split_file = file.path(PROJECT_PATH, "data/sigs/perturb-seq/pb_splits.csv"),
                                          split_pb_col = "pb",
                                          split_type = "10th",
                                          BPPARAM = bp)

  change_df <- rbind(k562_recon_eval_table,
                     rpe1_recon_eval_table)
  saveRDS(change_df, file.path(PROJECT_PATH, "results/eval/perturb-seq/ctrl_10th_mean_eval_table.rds"))
} else {
  change_df <- readRDS(file.path(PROJECT_PATH, "results/eval/perturb-seq/ctrl_10th_mean_eval_table.rds"))
}



## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
change_df %>% 
  dplyr::group_by(source) %>%
  summarize(NES_mean = mean(NES, na.rm=TRUE),
            jacc = mean(jacc))


## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
change_df %>% 
  dplyr::group_by(source) %>%
  summarize(NES = mean(NES, na.rm=TRUE),
            jacc = mean(jacc)) %>%
  summarize(NES_mean = mean(NES),
            jacc_mean = mean(jacc))


## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
combined_aggregated_df <- paired_eval_table(change_df, no_change_path)

combined_aggregated_df %>% 
  filter(source=="rpe1") %>%
  filter(target=="k562") %>%
  ggplot(aes(x=NES_FALSE, y=NES_TRUE)) +
  geom_point(aes(alpha=kept_alpha)) +
  geom_abline(color="black") +
  scale_alpha_continuous("Proportion of Source Kept") +
  labs(title = paste("NES Scores in", "rpe1 Up"), 
       x = "NES Score (Before)",
       y = "NES Score (After)")

combined_aggregated_df %>% 
  filter(source=="rpe1") %>%
  filter(target=="k562") %>%
  ggplot(aes(x=jacc_FALSE, y=jacc_TRUE)) +
  geom_point(aes(alpha=kept_alpha)) +
  geom_abline(color="black") +
  scale_alpha_continuous("Proportion of Source Kept") +
  labs(title = paste("Jacc. Scores in", "rpe1 Up"), 
       x = "Jac Score (Before)",
       y = "Jac Score (After)")

## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
combined_aggregated_df %>% 
  drop_na() %>%
  group_by(source) %>% 
  summarize(Jacc_Wilcox_pval = wilcox.test(jacc_TRUE, jacc_FALSE, paired = TRUE, alternative="greater")$p.value,
            NES_Wilcox_pval = wilcox.test(NES_TRUE, NES_FALSE, paired = TRUE, alternative="greater")$p.value) %>%
  summarize(NES_meta_p = fishers_meta_p(NES_Wilcox_pval),
            Jacc_meta_p = fishers_meta_p(Jacc_Wilcox_pval))


## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Pred Sig
k562_sig <- readRDS(file.path(DATAPATH, "mean/90th_perturb/k562_90_signatures.rds"))
rpe1_sig <- readRDS(file.path(DATAPATH, "mean/90th_perturb/rpe1_90_signatures.rds"))
rpe1_pred_sig <- list()
for(split in 1:length(rpe1_sig)) {
  rpe1_pred_sig[[split]] <- rep(list(rpe1_sig[[split]][["up"]]), length(rpe1_true_sig))
  names(rpe1_pred_sig[[split]]) <- names(rpe1_true_sig)  
}
names(rpe1_pred_sig) <- paste0("split_", 1:10)
k562_pred_sig <- list()
for(split in 1:length(rpe1_sig)) {
  k562_pred_sig[[split]] <- rep(list(k562_sig[[split]][["up"]]), length(k562_true_sig))
  names(k562_pred_sig[[split]]) <- names(k562_true_sig)  
}
names(k562_pred_sig) <- paste0("split_", 1:10)



## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
do_save <- TRUE
if(do_save) {
  k562_recon_eval_table <- sig_eval_table(source_sigs = k562_true_sig_up,
                                          pred_sigs = rpe1_pred_sig, 
                                          true_sigs = rpe1_true_sig,
                                          source = "k562",
                                          target = "rpe1",
                                          splits = TRUE,
                                          split_file = file.path(PROJECT_PATH, "data/sigs/perturb-seq/pb_splits.csv"),
                                          split_pb_col = "pb",
                                          split_type = "90th",
                                          BPPARAM = bp)
rpe1_recon_eval_table <- sig_eval_table(source_sigs = rpe1_true_sig_up,
                                        pred_sigs = k562_pred_sig, 
                                        true_sigs = k562_true_sig,
                                        source = "rpe1",
                                        target = "k562",
                                        splits = TRUE,
                                        split_file = file.path(PROJECT_PATH, "data/sigs/perturb-seq/pb_splits.csv"),
                                        split_pb_col = "pb",
                                        split_type = "90th",
                                        BPPARAM = bp)

  change_df <- rbind(k562_recon_eval_table,
                     rpe1_recon_eval_table)
  saveRDS(change_df, file.path(PROJECT_PATH, "results/eval/perturb-seq/ctrl_90th_mean_eval_table.rds"))
} else {
  change_df <- readRDS(file.path(PROJECT_PATH, "results/eval/perturb-seq/ctrl_90th_mean_eval_table.rds"))
}


## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
change_df %>% 
  dplyr::group_by(source) %>%
  summarize(NES_mean = mean(NES, na.rm=TRUE),
            jacc = mean(jacc))


## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
change_df %>% 
  dplyr::group_by(source) %>%
  summarize(NES = mean(NES, na.rm=TRUE),
            jacc = mean(jacc)) %>%
  summarize(NES_mean = mean(NES),
            jacc_mean = mean(jacc))


## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
combined_aggregated_df <- paired_eval_table(change_df, no_change_path)


## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
combined_aggregated_df %>% 
  group_by(source) %>% 
  summarize(Jacc_Wilcox_pval = wilcox.test(jacc_TRUE, jacc_FALSE, paired = TRUE, alternative="greater")$p.value,
            NES_Wilcox_pval = wilcox.test(NES_TRUE, NES_FALSE, paired = TRUE, alternative="greater")$p.value) %>%
  summarize(NES_meta_p = fishers_meta_p(NES_Wilcox_pval),
            Jacc_meta_p = fishers_meta_p(Jacc_Wilcox_pval))


## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
combined_aggregated_df %>% 
  filter(source=="rpe1") %>%
  filter(target=="k562") %>%
  ggplot(aes(x=NES_FALSE, y=NES_TRUE)) +
  geom_point(aes(alpha=kept_alpha)) +
  geom_abline(color="black") +
  scale_alpha_continuous("Proportion of Source Kept") +
  labs(title = paste("NES Scores in", "rpe1 Up"), 
       x = "NES Score (Before)",
       y = "NES Score (After)")

