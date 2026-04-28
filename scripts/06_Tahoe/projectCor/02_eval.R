library(tidyverse)
library(sigrecon)
library(BiocParallel)
library(doParallel)

cores <- 15
cl <- makeCluster(cores)
registerDoParallel(cl)

PATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/projects/sig_recon")

RECON_PATH <- file.path(PATH,"data/sigs/tahoe/projectCor")
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

# ------------------------------------------------
# Evaluation function
# ------------------------------------------------

run_eval_parallel <- function(regime, score){
  
  message("Running ", regime," ",score)
  
  combined_df <- foreach(
    cellline = targets,
    .combine = bind_rows,
    .packages = c("sigrecon","tidyverse"),
    .export = c("tahoe_sigs","targets","nci_h23_sig","PATH","RECON_PATH")
  ) %dopar% {
    
    true_sig <- tahoe_sigs[[cellline]]
    
    file <- file.path(
      RECON_PATH,
      paste0("nci_h23_",cellline,"_",regime,"_",score,".rds")
    )
    
    if(!file.exists(file)){
      return(NULL)
    }
    
    sig <- readRDS(file)
    if(regime == "ctrl") {
      recon_eval_table <- sig_eval_table(
        source_sigs = nci_h23_sig,
        pred_sigs = sig,
        true_sigs = true_sig,
        source = "nci_h23",
        target = cellline
      )
    } else {
      recon_eval_table <- sig_eval_table(
        source_sigs = nci_h23_sig,
        pred_sigs = sig,
        true_sigs = true_sig,
        source = "nci_h23",
        target = cellline,
        splits = regime %in% c("10th","90th"),
        split_file = file.path(PATH,"data/sigs/tahoe/drug_splits.csv"),
        split_pb_col = "drug",
        split_type = regime
      )
    }
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
      paste0(regime,"_projectcor_",score,".rds")
    )
  )
  
  combined_aggregated_df
}

# ------------------------------------------------
# Run experiments
# ------------------------------------------------

run_eval_parallel("ctrl","gsva")
run_eval_parallel("ctrl","eigen")

run_eval_parallel("10th","gsva")
run_eval_parallel("10th","eigen")

run_eval_parallel("90th","gsva")
run_eval_parallel("90th","eigen")

stopCluster(cl)