library(tidyverse)
library(sigrecon)
library(BiocParallel)

bp <- make_bpparam(workers = 15, RNGseed = 123, type = "multicore")
BiocParallel::register(bp)

PATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/projects/sig_recon")
RECON_PATH <- file.path(PATH, "data/sigs/drugmatrix/orthos")
SAVE_PATH <- file.path(PATH, "results/eval/drugmatrix")

dir.create(SAVE_PATH, recursive = TRUE, showWarnings = FALSE)

liver_sig_up <- lapply(drugmatrix.liver, function(x) x$up)
kidney_sig_up <- lapply(drugmatrix.kidney, function(x) x$up)

liver_orthos_sigs <- readRDS(file.path(RECON_PATH, "liver_orthos_sigs.rds"))
kidney_orthos_sigs <- readRDS(file.path(RECON_PATH, "kidney_orthos_sigs.rds"))

liver_orthos_up <- lapply(liver_orthos_sigs, function(x) x$up)
kidney_orthos_up <- lapply(kidney_orthos_sigs, function(x) x$up)

no_change_file <- if (file.exists(file.path(SAVE_PATH, "nochange_table.rds"))) {
  file.path(SAVE_PATH, "nochange_table.rds")
} else {
  file.path(SAVE_PATH, "no_change_eval_table.rds")
}

all_eval_no_change_table <- readRDS(no_change_file)
all_eval_no_change_table$source <- dplyr::case_when(
  all_eval_no_change_table$source == "liver" ~ "liver_kidney",
  all_eval_no_change_table$source == "kidney" ~ "kidney_liver",
  TRUE ~ all_eval_no_change_table$source
)

liver_kidney_eval_table <- sig_eval_table(
  source_sigs = liver_sig_up,
  pred_sigs = liver_orthos_up,
  true_sigs = drugmatrix.kidney,
  source = "liver_kidney",
  target = "kidney",
  BPPARAM = bp
)

kidney_liver_eval_table <- sig_eval_table(
  source_sigs = kidney_sig_up,
  pred_sigs = kidney_orthos_up,
  true_sigs = drugmatrix.liver,
  source = "kidney_liver",
  target = "liver",
  BPPARAM = bp
)

all_eval_table <- rbind(liver_kidney_eval_table, kidney_liver_eval_table)

saveRDS(
  all_eval_table,
  file.path(SAVE_PATH, "ctrl_orthos_residual_raw_eval_table.rds")
)

combined_aggregated_df <- merge(
  all_eval_no_change_table %>% dplyr::select(source, gene, NES, jacc),
  all_eval_table %>% mutate(kept_alpha = kept / (kept + displaced)),
  by = c("source", "gene"),
  suffixes = c("_FALSE", "_TRUE")
) %>%
  tibble()

saveRDS(
  combined_aggregated_df,
  file.path(SAVE_PATH, "ctrl_orthos_residual_eval_table.rds")
)
