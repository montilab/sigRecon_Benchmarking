# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

This repository benchmarks **signature recontextualization** methods — algorithms that take a gene signature from one biological context (cell line, tissue) and predict the corresponding signature in another context. It pairs an R package (`sigrecon/`) with analysis scripts for multiple datasets and a figures report.

## Repository Structure

- `sigrecon/` — R package (has its own `.git`); contains all reusable functions, benchmark metrics, and baseline method implementations
- `scripts/` — numbered pipeline scripts per dataset (02–08), each with subdirectories per method (`mean/`, `network_prop/`, `orthos/`, `projectCor/`, `scgpt/`, `stack/`)
- `results/eval/` — pre-computed `.rds` eval tables per dataset/method/split
- `results/figures.Rmd` — master figure-generation notebook; reads from `results/eval/`
- `drug_metadata.parquet` — drug annotation table used across scripts

## Dataset Numbering Convention

| Folder | Dataset |
|---|---|
| `02_PerturbSeq` | Perturb-seq (K562, RPE1) |
| `03_DrugMatrix` | DrugMatrix (kidney, liver) |
| `04_SciPlex` | SciPlex (K562, MCF7) |
| `05_Tabula_Sapiens` | Tabula Sapiens |
| `06_Tahoe` | Tahoe (multiple cell lines) |
| `07_GTEX` | GTEx (blood, brain) |
| `08_Neurips` | NeurIPS 2023 immune |

Each dataset folder has numbered scripts (`00_` = QC/EDA, `01_` = signature definition, `02_+` = method-specific runs), with `.sh` job scripts for HPC submission alongside each `.R`/`.py`.

## Key sigrecon Package Functions

- `sig_eval_table()` — main evaluation function; computes Jaccard, NES (via fgsea), and optionally ridge R² for source/predicted/true signatures
- `projectCor()` — baseline method: projects source signatures via GSVA/AUCell/eigengene scoring, then reranks genes by correlation
- `netProp()` — baseline method: one-call network propagation, learns a WGCNA co-expression network from target-context expression (`wgcna.adj()`) and propagates seed signatures across it (`network_sig()`); pass a pre-built `ig` to skip network construction and reuse a cached network across calls
- `recontextualize()` — single dispatcher across all three baseline methods (`"projectCor"`, `"networkProp"`, `"mean"`)
- `random_walk()` / `rwr_df()` — lower-level network propagation building blocks used by `network_sig()`/`netProp()`
- `paired_eval_table()` — formats eval output for paired plotting

## Evaluation Design

Methods are evaluated under three conditions stored in `results/eval/<dataset>/`:
- `ctrl_*` — none of the target perturbations are used.
- `ctrl_10th_*` — 1/10th of the target perturbations are used.
- `ctrl_90th_*` — 9/10th of the target perturbations are used.

Metrics in eval tables: `jacc_TRUE`/`jacc_FALSE` (Jaccard for predicted vs. true and source vs. true), `NES_TRUE`/`NES_FALSE` (fgsea NES). Better recontextualization = predicted closer to true than source is.

## R Package Development

```r
# Install with dependencies
BiocManager::install("montilab/sigrecon", dependencies=TRUE)

# Local dev install from repo root
devtools::install("sigrecon/")

# Run tests
devtools::test("sigrecon/")

# Run a single test file
testthat::test_file("sigrecon/tests/testthat/test-sig-eval-table.R")

# Regenerate documentation
devtools::document("sigrecon/")
```

## Rendering Figures

```r
rmarkdown::render("results/figures.Rmd")
```

Requires `PROJECT_PATH` env var pointing to the upstream data directory (set `MLAB` env var).

## Python Methods (scGPT, STACK)

Scripts in `scripts/02_PerturbSeq/scgpt/` and `scripts/*/stack/` are Python. Each has a corresponding `.sh` for HPC. Preprocessing produces AnnData `.h5ad` files; outputs are fed back into R eval scripts.
