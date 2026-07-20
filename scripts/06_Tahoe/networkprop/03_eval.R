library(tidyverse)
library(sigrecon)
library(doParallel)
registerDoParallel(15)

PATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/projects/sig_recon")

RECON_PATH <- file.path(PATH,"data/sigs/tahoe/networkprop")
SAVE_PATH  <- file.path(PATH,"results/eval/tahoe")

drug_splits <- read.csv(file.path(PATH,"data/sigs/tahoe/drug_splits.csv"))

combined_no_recon_df <- readRDS(
  file.path(PATH,"results/eval/tahoe/no_change_eval_table.rds")
)

# ------------------------------------------------
# Signatures
# ------------------------------------------------

nci_h23_sig <- lapply(tahoe.nci_h23,function(x)x$up)

tahoe_sigs <- list(
  "A498" = tahoe.a498,
  "HCT15" = tahoe.hct15,
  "HEC-1-A" = tahoe.hec_1_a,
  "LoVo" = tahoe.lovo,
  "MIA PaCa-2" = tahoe.miapaca_2,
  "Panc 03.27" = tahoe.panc03.27,
  "SNU-1" = tahoe.snu_1,
  "SNU-423" = tahoe.snu_423,
  "SW48" = tahoe.sw48
)

targets <- names(tahoe_sigs)

# # 1. Ctrl experiment
# 
# combined_df <- foreach(
#   cellline = targets,
#   .combine = bind_rows
# ) %dopar% {
#   
#   true_sig <- tahoe_sigs[[cellline]]
#   
#   file <- file.path(
#     RECON_PATH,
#     paste0(cellline,"_control.rds")
#   )
#   
#   sig <- readRDS(file)
# 
#   recon_eval_table <- sig_eval_table(
#     source_sigs = nci_h23_sig,
#     pred_sigs = sig,
#     true_sigs = true_sig,
#     source = "nci_h23",
#     target = cellline
#   )
#  
#   recon_eval_table
# }
# 
# combined_aggregated_df <- merge(
#   combined_no_recon_df %>% dplyr::select(source,gene,NES,jacc),
#   combined_df %>% mutate(kept_alpha = kept/(kept+displaced)),
#   by=c("source","gene"),
#   suffixes=c("_FALSE","_TRUE")
# ) %>% tibble()
# 
# saveRDS(
#   combined_aggregated_df,
#   file.path(
#     SAVE_PATH,
#     "ctrl_networkprop_.rds")
# )

# 1. Split experiments
pred_sigs_10th <- readRDS(file.path(RECON_PATH, "control_1_10th.rds"))
pred_sigs_90th <- readRDS(file.path(RECON_PATH, "control_9_10th.rds"))

fill_missing_splits <- function(pred_sigs, all_splits = paste0("split_", 1:10)) {
  lapply(pred_sigs, function(cellline_sigs) {
    out <- vector("list", length(all_splits))
    names(out) <- all_splits
    
    for (split in all_splits) {
      out[split] <- if (split %in% names(cellline_sigs)) {
        list(cellline_sigs[[split]])
      } else {
        list(NULL)
      }
    }
    
    out
  })
}

pred_sigs_10th <- fill_missing_splits(pred_sigs_10th)
pred_sigs_90th <- fill_missing_splits(pred_sigs_90th)

## 10th
combined_df <- foreach(
  cellline = targets,
  .combine = bind_rows
) %dopar% {
  
  true_sig <- tahoe_sigs[[cellline]]
  
  sig <- pred_sigs_10th[[cellline]]
  
  recon_eval_table <- sig_eval_table(
    source_sigs = nci_h23_sig,
    pred_sigs = sig,
    true_sigs = true_sig,
    source = "nci_h23",
    target = cellline,
    splits = TRUE,
    split_file = file.path(PATH,"data/sigs/tahoe/drug_splits.csv"),
    split_pb_col = "drug",
    split_type = "10th"
  )
  
  recon_eval_table
}

combined_aggregated_df <- merge(
  combined_no_recon_df %>% dplyr::select(source,gene,NES,jacc),
  combined_df %>% mutate(kept_alpha = kept/(kept+displaced)),
  by=c("source","gene"),
  suffixes=c("_FALSE","_TRUE")
) %>% tibble()

saveRDS(
  combined_aggregated_df,
  file.path(
    SAVE_PATH,
    "network_prop_10th.rds")
)

## 90th
combined_df <- foreach(
  cellline = targets,
  .combine = bind_rows
) %dopar% {
  
  true_sig <- tahoe_sigs[[cellline]]
  
  sig <- pred_sigs_90th[[cellline]]
  
  recon_eval_table <- sig_eval_table(
    source_sigs = nci_h23_sig,
    pred_sigs = sig,
    true_sigs = true_sig,
    source = "nci_h23",
    target = cellline,
    splits = TRUE,
    split_file = file.path(PATH,"data/sigs/tahoe/drug_splits.csv"),
    split_pb_col = "drug",
    split_type = "90th"
  )
  
  recon_eval_table
}

combined_aggregated_df <- merge(
  combined_no_recon_df %>% dplyr::select(source,gene,NES,jacc),
  combined_df %>% mutate(kept_alpha = kept/(kept+displaced)),
  by=c("source","gene"),
  suffixes=c("_FALSE","_TRUE")
) %>% tibble()

saveRDS(
  combined_aggregated_df,
  file.path(
    SAVE_PATH,
    "network_prop_90th.rds")
)