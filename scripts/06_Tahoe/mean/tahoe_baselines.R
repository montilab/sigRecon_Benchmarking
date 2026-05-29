## ----echo=FALSE, message=FALSE------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(tidyverse)
library(sigrecon)
library(foreach)
library(fgsea)
library(doParallel)

registerDoParallel(15)
library(BiocParallel)
do_save <- FALSE
PROJECT_PATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/projects/sig_recon/")


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
source_sig <-lapply(tahoe.nci_h23, function(x) x$up)
tahoe_sigs <- list(a498 = tahoe.a498,
                   hct15 = tahoe.hct15,
                   hec1a = tahoe.hec_1_a,
                   lovo = tahoe.lovo,
                   mia_paca2 = tahoe.miapaca_2,
                   panc_0327 = tahoe.panc03.27,
                   snu1 = tahoe.snu_1,
                   snu423 = tahoe.snu_423,
                   sw48 = tahoe.sw48)


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
if(do_save) {
  combined_no_recon_df <- foreach(sig_name = names(tahoe_sigs), .combine = rbind) %do% {
    print(sig_name)
    sig <- tahoe_sigs[[sig_name]]
    
    # # Need to find shared drugs because source may have more drugs that pass filtering
    # shared_drugs <- intersect(names(source_sig), names(sig))
    # sig <- sig[shared_drugs]
    # nci_h23_sig <- source_sig[shared_drugs]
    recon_eval_table <- sig_eval_table(source_sigs = nci_h23_sig,
                                       pred_sigs = nci_h23_sig,
                                       true_sigs = sig,
                                       source = "nci_h23",
                                       target = sig_name)
    recon_eval_table
  }
  saveRDS(combined_no_recon_df, file.path(PROJECT_PATH, "results/eval/tahoe/no_change_eval_table.rds"))
} else {
  combined_no_recon_df <- readRDS(file.path(PROJECT_PATH, "results/eval/tahoe/no_change_eval_table.rds"))
}



## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
combined_no_recon_df %>% 
  group_by(target) %>%
  summarize(NES = mean(NES),
            jacc = mean(jacc))

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
combined_no_recon_df %>% 
  group_by(target) %>%
  summarize(NES = mean(NES),
            jacc = mean(jacc)) %>%
  summarize(NES = mean(NES),
            jacc = mean(jacc))


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Pred Sig
control_sigs_path <- Sys.glob(file.path(PROJECT_PATH, "data/sigs/tahoe/mean/ctrl/*.rds"))
control_sigs_path <- control_sigs_path[-6] 
control_sigs <- lapply(control_sigs_path, function(x) readRDS(x)[[1]])
names(control_sigs) <- names(tahoe_sigs)

control_pred_sig <- list()
for(cellline in names(control_sigs)) {
  up_sig <- control_sigs[[cellline]][["up"]]
  control_pred_sig[[cellline]] <- rep(list(up_sig), length(source_sig))
  names(control_pred_sig[[cellline]]) <- names(source_sig)  
}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
if(do_save) {
  combined_df <- foreach(cellline = names(control_pred_sig), .combine = rbind) %dopar% {
    print(cellline)
    sig <- control_pred_sig[[cellline]]
    true_sig <- tahoe_sigs[[cellline]]
    nci_h23_sig <- source_sig
    recon_eval_table <- sig_eval_table(source_sigs = nci_h23_sig,
                                       pred_sigs = sig,
                                       true_sigs = true_sig,
                                       source = "nci_h23",
                                       target = cellline)
    recon_eval_table
  }
  saveRDS(combined_df, file.path(PROJECT_PATH, "results/eval/tahoe/ctrl_mean_eval_table.rds"))
} else {
  combined_df <- readRDS(file.path(PROJECT_PATH, "results/eval/tahoe/ctrl_mean_eval_table.rds"))
}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
combined_aggregated_df <- merge(combined_no_recon_df %>% dplyr::select("source", "target", "gene", "NES", "jacc"), 
                                combined_df %>% mutate(kept_alpha = kept/(kept+displaced)), 
                                by = c("source", "target", "gene"), 
                                suffixes = c("_FALSE", "_TRUE")) %>% tibble

saveRDS(combined_aggregated_df, file.path(PROJECT_PATH, "results/eval/tahoe/ctrl_mean_eval_table.rds"))
## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# plotting a negative fgsea run
ref_vec <- tahoe_sigs$a498$`8-Hydroxyquinoline`$up_full
data_vec <- control_pred_sig$a498$`8-Hydroxyquinoline`
n <- length(ref_vec)
stats <- seq(from = 1, to = -1, length.out = n)
names(stats) <- ref_vec
plotEnrichment(data_vec, stats)

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# plotting a negative fgsea run
ref_vec <- tahoe_sigs$a498$Afatinib$up_full
data_vec <- control_pred_sig$a498$Afatinib
n <- length(ref_vec)
stats <- seq(from = 1, to = -1, length.out = n)
names(stats) <- ref_vec
plotEnrichment(data_vec, stats)


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
combined_df %>% 
  group_by(target) %>%
  summarize(NES = mean(NES),
            jacc = mean(jacc))


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
combined_df %>% 
  group_by(target) %>%
  summarize(NES = mean(NES),
            jacc = mean(jacc)) %>%
  summarize(NES = mean(NES),
            jacc = mean(jacc))


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
combined_aggregated_df %>% 
  group_by(source) %>% 
  summarize(Jacc_Wilcox_pval = wilcox.test(jacc_TRUE, jacc_FALSE, paired = TRUE, alternative="greater")$p.value,
            KS_Wilcox_pval = wilcox.test(NES_TRUE, NES_FALSE, paired = TRUE, alternative="greater")$p.value)


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
combined_aggregated_df %>% 
  group_by(source) %>% 
  summarize(Jacc_Wilcox_pval = wilcox.test(jacc_TRUE, jacc_FALSE, paired = TRUE, alternative="greater")$p.value,
            KS_Wilcox_pval = wilcox.test(NES_TRUE, NES_FALSE, paired = TRUE, alternative="greater")$p.value) %>%
  summarize(KS_meta_p = fishers_meta_p(KS_Wilcox_pval),
            Jacc_meta_p = fishers_meta_p(Jacc_Wilcox_pval))


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Pred Sig
control_10th_sigs_path <- Sys.glob(file.path(PROJECT_PATH, "data/sigs/tahoe/mean/10th_perturb/*.rds"))
control_10th_sigs_path <- control_10th_sigs_path[-6] 
control_10th_sigs <- lapply(control_10th_sigs_path, function(x) readRDS(x))
names(control_10th_sigs) <- names(tahoe_sigs)

control_10th_pred_sigs <- list()
for(cellline in names(control_10th_sigs)) {
  for(split in 1:10) {
    split_col <- paste0("split_", split) 
    up_sig <- control_10th_sigs[[cellline]][[split_col]][["up"]]
    control_10th_pred_sigs[[cellline]][[split_col]] <- rep(list(up_sig), length(source_sig))
    names(control_10th_pred_sigs[[cellline]][[split_col]]) <- names(source_sig)  
  }
}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
if(do_save) {
  combined_df <- foreach(cellline = names(control_10th_pred_sigs), .combine = rbind) %dopar% {
    print(cellline)
    split_sig <- control_10th_pred_sigs[[cellline]]
    true_sig <- tahoe_sigs[[cellline]]
    nci_h23_sig <- source_sig
    recon_eval_table <- sig_eval_table(source_sigs = nci_h23_sig,
                                       pred_sigs = split_sig,
                                       true_sigs = true_sig,
                                       source = "nci_h23",
                                       target = cellline,
                                       splits=TRUE,
                                       split_file = file.path(PROJECT_PATH, "data/sigs/tahoe/drug_splits.csv"),
                                       split_pb_col = "drug",
                                       split_type = "10th"
                                       
                                       )
    recon_eval_table
  }
  saveRDS(combined_df, file.path(PROJECT_PATH, "results/eval/tahoe/ctrl_10th_mean_eval_table.rds"))
} else {
  combined_df <- readRDS(file.path(PROJECT_PATH, "results/eval/tahoe/ctrl_10th_mean_eval_table.rds"))
}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
combined_aggregated_df <- merge(combined_no_recon_df %>% dplyr::select("source", "target", "gene", "NES", "jacc"), 
                                combined_df %>% mutate(kept_alpha = kept/(kept+displaced)), 
                                by = c("source", "target", "gene"), 
                                suffixes = c("_FALSE", "_TRUE")) %>% tibble

saveRDS(combined_aggregated_df, file.path(PROJECT_PATH, "results/eval/tahoe/ctrl_10th_mean_eval_table.rds"))

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
combined_df %>% 
  group_by(target) %>%
  summarize(NES = mean(NES),
            jacc = mean(jacc))

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
combined_df %>% 
  group_by(target) %>%
  summarize(NES = mean(NES),
            jacc = mean(jacc)) %>%
  summarize(NES = mean(NES),
            jacc = mean(jacc))


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
combined_aggregated_df %>% 
  group_by(source) %>% 
  summarize(Jacc_Wilcox_pval = wilcox.test(jacc_TRUE, jacc_FALSE, paired = TRUE, alternative="greater")$p.value,
            KS_Wilcox_pval = wilcox.test(NES_TRUE, NES_FALSE, paired = TRUE, alternative="greater")$p.value)


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
combined_aggregated_df %>% 
  group_by(source) %>% 
  summarize(Jacc_Wilcox_pval = wilcox.test(jacc_TRUE, jacc_FALSE, paired = TRUE, alternative="greater")$p.value,
            KS_Wilcox_pval = wilcox.test(NES_TRUE, NES_FALSE, paired = TRUE, alternative="greater")$p.value) %>%
  summarize(KS_meta_p = fishers_meta_p(KS_Wilcox_pval),
            Jacc_meta_p = fishers_meta_p(Jacc_Wilcox_pval))


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Pred Sig
control_90th_sigs_path <- Sys.glob(file.path(PROJECT_PATH, "data/sigs/tahoe/mean/90th_perturb/*.rds"))
control_90th_sigs_path <- control_90th_sigs_path[-6] 
control_90th_sigs <- lapply(control_90th_sigs_path, function(x) readRDS(x))
names(control_90th_sigs) <- names(tahoe_sigs)

control_90th_pred_sigs <- list()
for(cellline in names(control_90th_sigs)) {
  for(split in 1:10) {
    split_col <- paste0("split_", split) 
    split_col_name <- paste0(split_col, "_excluded")
    up_sig <- control_90th_sigs[[cellline]][[split_col_name]][["up"]]
    control_90th_pred_sigs[[cellline]][[split_col]] <- rep(list(up_sig), length(source_sig))
    names(control_90th_pred_sigs[[cellline]][[split_col]]) <- names(source_sig)  
  }
}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
do_save <- TRUE
if(do_save) {
  combined_df <- foreach(cellline = names(control_90th_pred_sigs), .combine = rbind) %dopar% {
    print(cellline)
    split_sig <- control_90th_pred_sigs[[cellline]]
    true_sig <- tahoe_sigs[[cellline]]
    nci_h23_sig <- source_sig
    recon_eval_table <- sig_eval_table(source_sigs = nci_h23_sig,
                                       pred_sigs = split_sig,
                                       true_sigs = true_sig,
                                       source = "nci_h23",
                                       target = cellline,
                                       splits=TRUE,
                                       split_file = file.path(PROJECT_PATH, "data/sigs/tahoe/drug_splits.csv"),
                                       split_pb_col = "drug",
                                       split_type = "90th"
                                       )
    recon_eval_table
  }
  saveRDS(combined_df, file.path(PROJECT_PATH, "results/eval/tahoe/ctrl_90th_mean_eval_table.rds"))
} else {
  combined_df <- readRDS(file.path(PROJECT_PATH, "results/eval/tahoe/ctrl_90th_mean_eval_table.rds"))
}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
combined_aggregated_df <- merge(combined_no_recon_df %>% dplyr::select("source", "target", "gene", "NES", "jacc"), 
                                combined_df %>% mutate(kept_alpha = kept/(kept+displaced)), 
                                by = c("source", "target", "gene"), 
                                suffixes = c("_FALSE", "_TRUE")) %>% tibble

saveRDS(combined_aggregated_df, file.path(PROJECT_PATH, "results/eval/tahoe/ctrl_90th_mean_eval_table.rds"))


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
combined_df %>% 
  group_by(target) %>%
  summarize(NES = mean(NES),
            jacc = mean(jacc))


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
combined_df %>% 
  group_by(target) %>%
  summarize(NES = mean(NES),
            jacc = mean(jacc)) %>%
  summarize(NES = mean(NES),
            jacc = mean(jacc))


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
combined_aggregated_df %>% 
  group_by(source) %>% 
  summarize(Jacc_Wilcox_pval = wilcox.test(jacc_TRUE, jacc_FALSE, paired = TRUE, alternative="greater")$p.value,
            KS_Wilcox_pval = wilcox.test(NES_TRUE, NES_FALSE, paired = TRUE, alternative="greater")$p.value)


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
combined_aggregated_df %>% 
  group_by(source) %>% 
  summarize(Jacc_Wilcox_pval = wilcox.test(jacc_TRUE, jacc_FALSE, paired = TRUE, alternative="greater")$p.value,
            KS_Wilcox_pval = wilcox.test(NES_TRUE, NES_FALSE, paired = TRUE, alternative="greater")$p.value) %>%
  summarize(KS_meta_p = fishers_meta_p(KS_Wilcox_pval),
            Jacc_meta_p = fishers_meta_p(Jacc_Wilcox_pval))


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
combined_aggregated_df %>% 
  filter(target=="a498") %>%
  ggplot(aes(x=NES_FALSE, y=NES_TRUE)) +
  geom_point(aes(alpha=kept_alpha)) +
  geom_abline(color="black") +
  scale_alpha_continuous("Proportion of Source Kept") +
  labs(title = paste("NES Scores in", "a498 Up"), 
       x = "NES Score (Before)",
       y = "NES Score (After)")

