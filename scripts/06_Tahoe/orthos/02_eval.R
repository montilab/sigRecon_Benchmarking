library(tidyverse)
library(sigrecon)
library(BiocParallel)
library(doParallel)

cores <- 15
cl <- makeCluster(cores)
registerDoParallel(cl)

PATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/projects/sig_recon")
SAVE_PATH <- file.path(PATH, "results/eval/tahoe/")

dir.create(SAVE_PATH, recursive = TRUE, showWarnings = FALSE)

combined_no_recon_df <- readRDS(
  file.path(PATH, "results/eval/tahoe/no_change_eval_table.rds")
)

orthos_sigs <- readRDS(file.path(PATH, "data/sigs/tahoe/orthos/orthos_sigs.rds"))
orthos_sigs <- lapply(orthos_sigs, function(x) x$up)
nci_h23_sig <- lapply(tahoe.nci_h23, function(x) x$up)

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

celllines <- read.csv(file.path(PATH, "data/sigs/tahoe/cell_lines.csv"))$x
targets <- intersect(celllines[2:10], names(tahoe_sigs))

combined_df <- foreach(
  cellline = targets,
  .combine = bind_rows,
  .packages = c("sigrecon", "tidyverse", "BiocParallel"),
  .export = c("tahoe_sigs", "nci_h23_sig", "orthos_sigs")
) %dopar% {
  sig_eval_table(
    source_sigs = nci_h23_sig,
    pred_sigs = orthos_sigs,
    true_sigs = tahoe_sigs[[cellline]],
    source = "nci_h23",
    target = cellline,
    BPPARAM = BiocParallel::SerialParam()
  )
}

combined_aggregated_df <- merge(
  combined_no_recon_df %>% dplyr::select(source, gene, NES, jacc),
  combined_df %>% mutate(kept_alpha = kept / (kept + displaced)),
  by = c("source", "gene"),
  suffixes = c("_FALSE", "_TRUE")
) %>%
  tibble()

saveRDS(
  combined_aggregated_df,
  file.path(SAVE_PATH, "ctrl_orthos_residual_eval_table.rds")
)

stopCluster(cl)
