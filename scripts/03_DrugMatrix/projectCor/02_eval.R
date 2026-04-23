library(tidyverse)
library(sigrecon)
library(BiocParallel)
bp <- make_bpparam(workers = 15, RNGseed = 123, type="multicore")
BiocParallel::register(bp)
PATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/projects/sig_recon")
DATA_PATH <- file.path(Sys.getenv("AGED"), "CBMrepositoryData/perturbational_data/drugmatrix/")
RECON_PATH <- file.path(PATH, "data/sigs/drugmatrix/projectCor/")
SAVE_PATH <- file.path(PATH, "results/eval/drugmatrix")

# Load signatures
kidney_sig_up <- lapply(drugmatrix.kidney, function(x) x$up)
liver_sig_up <- lapply(drugmatrix.liver, function(x) x$up)

# Loading drug splits
drug_splits_path <- file.path(PATH, "data/sigs/drugmatrix/drug_splits.csv")

all_eval_no_change_table <- readRDS(file.path(SAVE_PATH, "nochange_table.rds"))

# 1. Control experiment
liver_eigen_pred_sigs <- readRDS(file.path(RECON_PATH, "liver_kidney_ctrl_eigen.rds"))
liver_gsva_pred_sigs <- readRDS(file.path(RECON_PATH, "liver_kidney_ctrl_gsva.rds"))
kidney_eigen_pred_sigs <- readRDS(file.path(RECON_PATH, "kidney_liver_ctrl_eigen.rds"))
kidney_gsva_pred_sigs <- readRDS(file.path(RECON_PATH, "kidney_liver_ctrl_gsva.rds"))

liver_kidney_eval_table <- sig_eval_table(source_sigs = liver_sig_up,
                                          pred_sigs = liver_eigen_pred_sigs,
                                          true_sigs = drugmatrix.kidney,
                                          source = "liver_kidney",
                                          BPPARAM = bp)
kidney_liver_eval_table <- sig_eval_table(source_sigs = kidney_sig_up,
                                          pred_sigs = kidney_eigen_pred_sigs,
                                          true_sigs = drugmatrix.liver,
                                          source = "kidney_liver",
                                          BPPARAM = bp)
all_eval_table <- rbind(liver_kidney_eval_table,
                        kidney_liver_eval_table)
combined_aggregated_df <- merge(all_eval_no_change_table %>% dplyr::select("source", "gene", "NES", "jacc"),
                                all_eval_table %>% mutate(kept_alpha = kept/(kept+displaced)),
                                by = c("source", "gene"),
                                suffixes = c("_FALSE", "_TRUE")) %>% tibble

saveRDS(combined_aggregated_df, file.path(SAVE_PATH, "ctrl_projectcor_eigen_eval_table.rds"))

liver_kidney_eval_table <- sig_eval_table(source_sigs = liver_sig_up,
                                          pred_sigs = liver_gsva_pred_sigs,
                                          true_sigs = drugmatrix.kidney,
                                          source = "liver_kidney",
                                          BPPARAM = bp)
kidney_liver_eval_table <- sig_eval_table(source_sigs = kidney_sig_up,
                                          pred_sigs = kidney_gsva_pred_sigs,
                                          true_sigs = drugmatrix.liver,
                                          source = "kidney_liver",
                                          BPPARAM = bp)
all_eval_table <- rbind(liver_kidney_eval_table,
                        kidney_liver_eval_table)
combined_aggregated_df <- merge(all_eval_no_change_table %>% dplyr::select("source", "gene", "NES", "jacc"),
                                all_eval_table %>% mutate(kept_alpha = kept/(kept+displaced)),
                                by = c("source", "gene"),
                                suffixes = c("_FALSE", "_TRUE")) %>% tibble
saveRDS(combined_aggregated_df, file.path(SAVE_PATH, "ctrl_projectcor_gsva_eval_table.rds"))

# ------------------------------------------------
# 2. 1/10th experiment
# ------------------------------------------------

for(regime in c("gsva","eigen")){
  
  message("Processing 1/10th ", regime)
  
  kidney_paths <- Sys.glob(file.path(RECON_PATH,
                                     paste0("kidney_liver_10th_", regime, "_split_*.rds")))
  liver_paths  <- Sys.glob(file.path(RECON_PATH,
                                     paste0("liver_kidney_10th_", regime, "_split_*.rds")))
  
  split_names <- sub(".*(split_[0-9]+).*", "\\1", kidney_paths)
  kidney_sigs <- lapply(kidney_paths, readRDS)
  names(kidney_sigs) <- split_names
  
  split_names <- sub(".*(split_[0-9]+).*", "\\1", liver_paths)
  liver_sigs <- lapply(liver_paths, readRDS)
  names(liver_sigs) <- split_names
  
  liver_kidney_eval_table <- sig_eval_table(
    source_sigs = liver_sig_up,
    pred_sigs = liver_sigs,
    true_sigs = drugmatrix.kidney,
    source = "liver_kidney",
    splits = TRUE,
    split_file = drug_splits_path,
    split_pb_col = "drug",
    split_type = "10th",
    BPPARAM = bp
  )
  
  kidney_liver_eval_table <- sig_eval_table(
    source_sigs = kidney_sig_up,
    pred_sigs = kidney_sigs,
    true_sigs = drugmatrix.liver,
    source = "kidney_liver",
    splits = TRUE,
    split_file = drug_splits_path,
    split_pb_col = "drug",
    split_type = "10th",
    BPPARAM = bp
  )
  
  all_eval_table <- rbind(liver_kidney_eval_table,
                          kidney_liver_eval_table)
  
  combined_aggregated_df <- merge(
    all_eval_no_change_table %>% dplyr::select("source","gene","NES","jacc"),
    all_eval_table %>% mutate(kept_alpha = kept/(kept+displaced)),
    by=c("source","gene"),
    suffixes=c("_FALSE","_TRUE")
  ) %>% tibble()
  
  saveRDS(combined_aggregated_df,
          file.path(SAVE_PATH,
                    paste0("10th_projectcor_",regime,"_eval_table.rds")))
}


# ------------------------------------------------
# 3. 9/10th experiment
# ------------------------------------------------

for(regime in c("gsva","eigen")){
  
  message("Processing 9/10th ", regime)
  
  kidney_paths <- Sys.glob(file.path(RECON_PATH,
                                     paste0("kidney_liver_90th_", regime, "_split_*.rds")))
  liver_paths  <- Sys.glob(file.path(RECON_PATH,
                                     paste0("liver_kidney_90th_", regime, "_split_*.rds")))
  
  split_names <- sub(".*(split_[0-9]+).*", "\\1", kidney_paths)
  kidney_sigs <- lapply(kidney_paths, readRDS)
  names(kidney_sigs) <- split_names
  
  split_names <- sub(".*(split_[0-9]+).*", "\\1", liver_paths)
  liver_sigs <- lapply(liver_paths, readRDS)
  names(liver_sigs) <- split_names
  
  liver_kidney_eval_table <- sig_eval_table(
    source_sigs = liver_sig_up,
    pred_sigs = liver_sigs,
    true_sigs = drugmatrix.kidney,
    source = "liver_kidney",
    splits = TRUE,
    split_file = drug_splits_path,
    split_pb_col = "drug",
    split_type = "90th",
    BPPARAM = bp
  )
  
  kidney_liver_eval_table <- sig_eval_table(
    source_sigs = kidney_sig_up,
    pred_sigs = kidney_sigs,
    true_sigs = drugmatrix.liver,
    source = "kidney_liver",
    splits = TRUE,
    split_file = drug_splits_path,
    split_pb_col = "drug",
    split_type = "90th",
    BPPARAM = bp
  )
  
  all_eval_table <- rbind(liver_kidney_eval_table,
                          kidney_liver_eval_table)
  
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

combined_aggregated_df <- readRDS(file.path(SAVE_PATH, "ctrl_projectcor_gsva_eval_table.rds"))
combined_aggregated_df <- readRDS(file.path(SAVE_PATH, "ctrl_projectcor_eigen_eval_table.rds"))
combined_aggregated_df <- readRDS(file.path(SAVE_PATH, "10th_projectcor_gsva_eval_table.rds"))
combined_aggregated_df <- readRDS(file.path(SAVE_PATH, "10th_projectcor_eigen_eval_table.rds"))
combined_aggregated_df <- readRDS(file.path(SAVE_PATH, "90th_projectcor_gsva_eval_table.rds"))
combined_aggregated_df <- readRDS(file.path(SAVE_PATH, "90th_projectcor_eigen_eval_table.rds"))
combined_aggregated_df %>%
  dplyr::group_by(source) %>%
  dplyr::summarise(NES_mean = mean(NES_TRUE, na.rm = TRUE),
                   jacc_mean = mean(jacc_TRUE)) %>%
  dplyr::summarize(NES_mean = mean(NES_mean),
                   jacc_mean = mean(jacc_mean))
combined_aggregated_df %>%
  group_by(source) %>%
  summarize(Jacc_Wilcox_pval = wilcox.test(jacc_TRUE, jacc_FALSE, paired = TRUE, alternative="greater")$p.value,
            KS_Wilcox_pval = wilcox.test(NES_TRUE, NES_FALSE, paired = TRUE, alternative="greater")$p.value) %>%
  summarize(KS_meta_p = fishers_meta_p(KS_Wilcox_pval),
            Jacc_meta_p = fishers_meta_p(Jacc_Wilcox_pval))
