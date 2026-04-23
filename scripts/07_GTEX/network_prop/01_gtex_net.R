library(tidyverse)
library(sigrecon)
library(igraph)
library(SummarizedExperiment)
library(doParallel)
detectCores()
registerDoParallel(15)
DATA_PATH <- "/restricted/projectnb/gtex/montilab/CBMrepositoryData/GTEX_dbGap/processed_data"
PATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/projects/sig_recon")

filepaths <- Sys.glob(paste0(DATA_PATH, "/*.rds"))

foreach(filepath=filepaths) %dopar% {
  filename <- basename(filepath)
  tissue_name <- sub("^RNAseq_GTEx_v10_(.*)\\.rds$", "\\1", filename)


  # 1. Load Data
  se <- readRDS(filepath)

  if(ncol(se) < 30) {
    print(paste0("Not enough samples for ", tissue_name))
  } else {
    print(paste0("Processing ", tissue_name))

    se_mad_genes <- sigrecon:::rank.var.eset(se)$values[1:10000]
    se <- se[se_mad_genes,]
    se_mat <- assay(se)
    # rownames(se) <- rowData(se)$gene_id

    # 2. Learn WGCNA Networks
    se_ig <- sigrecon::wgcna.adj(t(se_mat), cor.type = "signed hybrid", igraph = TRUE, diag_zero = TRUE)
    saveRDS(se_ig, file.path(PATH, paste0("data/wgcna_networks/gtex/",tissue_name,".rds")))
  }
}
