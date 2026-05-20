library(argparse)
library(Seurat)
library(tidyverse)
library(sigrecon)
library(stringdist)
library(foreach)
library(doParallel)

registerDoParallel(15)

DATA_PATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/projects/sig_recon/data/sigs/tahoe/stack")
SAVE_PATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/projects/sig_recon/results/eval/tahoe")
SC_PATH <- file.path(Sys.getenv("AGED"), "CBMrepositoryData/perturbational_data/tahoe")

# parser <- ArgumentParser()
# parser$add_argument("--split_type", required=TRUE)
# args <- parser$parse_args()
# split_type <- args$split_type
split_type <- "10th"
nci_sig <- sigrecon::tahoe.nci_h23

target_sig_paths <- Sys.glob(file.path(DATA_PATH, "tahoe*.rds"))

split_path <- file.path(Sys.getenv("MLAB"),
                        "projects/brcameta/projects/sig_recon/data/sigs/tahoe/drug_splits.csv")

all_eval_table <- foreach(
  target_sig_path = target_sig_paths,
  .combine = dplyr::bind_rows,
  .packages = c("tidyverse","sigrecon")
) %dopar% {
  
  pred_sigs <- readRDS(target_sig_path)
  
  cell_line <- str_match(basename(target_sig_path),
                         "^tahoe_(.*?)_")[,2]
  
  message("Processing ", cell_line)
  
  if(split_type == "ctrl"){
    
    eval_tbl <- sig_eval_table(
      source_sigs = lapply(nci_sig, \(x) x$up),
      pred_sigs   = lapply(pred_sigs, \(x) x$up),
      true_sigs   = nci_sig,
      source      = paste0("nci_h23_", cell_line)
    )
    
    eval_tbl
    
  } else {
    
    split_tables <- lapply(names(pred_sigs), function(split_name){
      
      split_pred <- pred_sigs[[split_name]]
      
      sig_eval_table(
        source_sigs = lapply(nci_sig, \(x) x$up),
        pred_sigs   = lapply(split_pred, \(x) x$up),
        true_sigs   = nci_sig,
        source      = paste0("nci_h23_", cell_line),
        splits = TRUE,
        split_file  = split_path,
        split_type  = split_type
      )
      
    })
    
    bind_rows(split_tables)
  }
}

saveRDS(all_eval_table,
        file.path(SAVE_PATH,
                  paste0("stack_",split_type,"_eval_table.rds")))

no_change_eval_table <- readRDS("/restricted/projectnb/brcameta/projects/sig_recon/results/eval/tahoe/no_change_eval_table.rds")
no_change_eval_table$source <- paste0(no_change_eval_table$source, "_",toupper(no_change_eval_table$target))
combined_aggregated_df <- merge(
  no_change_eval_table %>% dplyr::select(source,gene,NES,jacc),
  all_eval_table %>% mutate(kept_alpha = kept/(kept+displaced)),
  by=c("source","gene"),
  suffixes=c("_FALSE","_TRUE")
) %>% tibble()

saveRDS(combined_aggregated_df,
        file.path(SAVE_PATH,
                  paste0("stack_",split_type,"_eval_table.rds")))
