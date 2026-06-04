library(tidyverse)
library(sigrecon)
library(BiocParallel)

bp <- make_bpparam(workers = 15, RNGseed = 123, type = "multicore")
BiocParallel::register(bp)

PATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/projects/sig_recon")
RECON_PATH <- file.path(PATH, "data/sigs/perturb-seq/orthos")
SAVE_PATH <- file.path(PATH, "results/eval/perturb-seq")

dir.create(SAVE_PATH, recursive = TRUE, showWarnings = FALSE)

k562_sig_up <- lapply(perturbseq.k562, function(x) x$up)
rpe1_sig_up <- lapply(perturbseq.rpe1, function(x) x$up)

k562_orthos_sigs <- readRDS(file.path(RECON_PATH, "k562_orthos_sigs.rds"))
rpe1_orthos_sigs <- readRDS(file.path(RECON_PATH, "rpe1_orthos_sigs.rds"))

k562_orthos_up <- lapply(k562_orthos_sigs, function(x) x$up)
rpe1_orthos_up <- lapply(rpe1_orthos_sigs, function(x) x$up)

all_eval_no_change_table <- readRDS(file.path(SAVE_PATH, "no_change_eval_table.rds"))
all_eval_no_change_table$source <- dplyr::case_when(
  all_eval_no_change_table$source == "k562" ~ "k562_rpe1",
  all_eval_no_change_table$source == "rpe1" ~ "rpe1_k562",
  TRUE ~ all_eval_no_change_table$source
)

k562_rpe1_eval_table <- sig_eval_table(
  source_sigs = k562_sig_up,
  pred_sigs = k562_orthos_up,
  true_sigs = perturbseq.rpe1,
  source = "k562_rpe1",
  target = "rpe1",
  BPPARAM = bp
)

rpe1_k562_eval_table <- sig_eval_table(
  source_sigs = rpe1_sig_up,
  pred_sigs = rpe1_orthos_up,
  true_sigs = perturbseq.k562,
  source = "rpe1_k562",
  target = "k562",
  BPPARAM = bp
)

all_eval_table <- rbind(k562_rpe1_eval_table, rpe1_k562_eval_table)

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
