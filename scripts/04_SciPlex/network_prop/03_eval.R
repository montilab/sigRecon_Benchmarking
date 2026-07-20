library(tidyverse)
library(sigrecon)
library(BiocParallel)

bp <- make_bpparam(workers = 15, RNGseed = 123, type="multicore")
register(bp)

PATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/projects/sig_recon")
RECON_PATH <- file.path(PATH, "data/sigs/sciplex/networkprop/")
SAVE_PATH <- file.path(PATH, "results/eval/sciplex")

drug_splits_path <- file.path(PATH, "data/sigs/sciplex/drug_splits.csv")

# ------------------------------------------------
# Load true signatures
# ------------------------------------------------

a549_true_sigs <- sciplex.a549
k562_true_sigs <- sciplex.k562
mcf7_true_sigs <- sciplex.mcf7

a549_sig_up <- lapply(a549_true_sigs, function(x) x$up)
k562_sig_up <- lapply(k562_true_sigs, function(x) x$up)
mcf7_sig_up <- lapply(mcf7_true_sigs, function(x) x$up)

true_sigs <- list(
  a549 = a549_true_sigs,
  k562 = k562_true_sigs,
  mcf7 = mcf7_true_sigs
)

source_sigs <- list(
  a549 = a549_sig_up,
  k562 = k562_sig_up,
  mcf7 = mcf7_sig_up
)

cell_lines <- c("a549","k562","mcf7")

pairs <- expand.grid(source = cell_lines,
                     target = cell_lines,
                     stringsAsFactors = FALSE) %>%
  filter(source != target)

# ------------------------------------------------
# No Change baseline
# ------------------------------------------------

no_change_path <- file.path(SAVE_PATH, "no_change_eval_table.rds")

# ------------------------------------------------
# Helper to run evaluation
# ------------------------------------------------

eval_experiment <- function(file_pattern, split_type = NULL){
  
  eval_tables <- list()
  
  for(i in seq_len(nrow(pairs))){
    
    src <- pairs$source[i]
    tgt <- pairs$target[i]
    
    pattern <- paste0(src,"_",tgt,"_",file_pattern)
    sigs_path <- Sys.glob(file.path(RECON_PATH, pattern))
    pred_sigs <- readRDS(sigs_path)
    
    if(length(sigs_path) == 0) next
    
    # control experiment
    if(is.null(split_type)){
      
      eval_df <- sig_eval_table(
        source_sigs = source_sigs[[src]],
        pred_sigs   = pred_sigs,
        true_sigs   = true_sigs[[tgt]],
        source      = paste0(src, "_", tgt),
        BPPARAM     = bp
      )
      
    } else {
      
      eval_df <- sig_eval_table(
        source_sigs = source_sigs[[src]],
        pred_sigs   = pred_sigs,
        true_sigs   = true_sigs[[tgt]],
        source      = paste0(src, "_", tgt),
        splits      = TRUE,
        split_file  = drug_splits_path,
        split_pb_col = "drug",
        split_type  = split_type,
        BPPARAM     = bp
      )
    }
    
    eval_tables[[paste(src,tgt,sep="_")]] <- eval_df
  }
  
  bind_rows(eval_tables)
} 

# # ------------------------------------------------
# # 1. CONTROL experiment
# # ------------------------------------------------
#
# ctrl_eval <- eval_experiment("control.rds")
#
# combined_ctrl <- paired_eval_table(ctrl_eval, no_change_path)
#
# saveRDS(combined_ctrl,
#         file.path(SAVE_PATH,"ctrl_networkprop_eval_table.rds"))

# # ------------------------------------------------
# # 2. 1/10th experiment
# # ------------------------------------------------
#
# ten_eval <- eval_experiment("*1_10th.rds","10th")
#
# combined_10th <- paired_eval_table(ten_eval, no_change_path)
#
# saveRDS(combined_10th,
#         file.path(SAVE_PATH,"10th_networkprop_eval_table.rds"))

# ------------------------------------------------
# 3. 9/10th experiment
# ------------------------------------------------

ninety_eval <- eval_experiment("*9_10th.rds","90th")

combined_90th <- paired_eval_table(ninety_eval, no_change_path)

saveRDS(combined_90th,
        file.path(SAVE_PATH,"90th_networkprop_eval_table.rds"))