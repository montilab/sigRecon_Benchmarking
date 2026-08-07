# sigRecon Benchmarking

This repository contains the benchmarking pipeline and figure-generation code for the sigRecon manuscript. It evaluates **signature recontextualization** methods — algorithms that take a gene signature from one biological context (cell line, tissue) and predict the corresponding signature in another context.

The reusable R package, including method implementations and evaluation utilities, lives at [montilab/sigrecon](https://github.com/montilab/sigrecon).

## Overview

Methods are benchmarked across four datasets under three data-availability regimes:

| Regime | Description |
|---|---|
| Control (null) | No target-context perturbations used |
| Low coverage (1/10) | 1/10th of target perturbations available |
| High coverage (9/10) | 9/10th of target perturbations available |

Performance is measured by Jaccard similarity and fgsea NES between predicted and ground-truth target signatures, reported as improvement over the source baseline.

## Repository Structure

```
sigrecon/          R package (methods + evaluation functions)
scripts/           Per-dataset pipeline scripts (02–08), numbered by dataset
  02_PerturbSeq/
  03_DrugMatrix/
  04_SciPlex/
  06_Tahoe/
results/
  eval/            Pre-computed evaluation tables (.rds) per dataset/method/regime
  figures.Rmd      Master figure-generation notebook (reads from results/eval/)
drug_metadata.parquet  Drug annotation table used across scripts
```

## Datasets

| Folder | Dataset |
|---|---|
| `02_PerturbSeq` | Perturb-seq (K562, RPE1) |
| `03_DrugMatrix` | DrugMatrix (kidney, liver) |
| `04_SciPlex` | SciPlex (K562, MCF7) |
| `06_Tahoe` | Tahoe (multiple cell lines) |

## Methods Benchmarked

| Method | Description |
|---|---|
| Mean | Context mean as predicted signature |
| NetProp | Network propagation via WGCNA co-expression graph (`netProp()`) |
| projCor-Eigen / projCor-GSVA | Projection-based scoring then gene reranking (`projectCor()`) |
| Orthos | Neural network–based context transfer |
| scGPT | Foundation model–based signature prediction |
| Stack | Stacked regression baseline |

## Reproducing Figures

```r
# From repo root
rmarkdown::render("results/figures.Rmd")
```

Requires the `MLAB` environment variable pointing to the upstream data directory. Pre-computed evaluation tables in `results/eval/` are sufficient to regenerate all paper figures without re-running the full pipeline.

## Data

Full-size pseudobulk expression and perturbational signatures for each dataset are archived on Zenodo:

| Dataset | Perturbational signatures | Pseudobulk expression |
|---|---|---|
| DrugMatrix | [10.5281/zenodo.21432933](https://doi.org/10.5281/zenodo.21432933) | [10.5281/zenodo.21433031](https://doi.org/10.5281/zenodo.21433031) |
| SciPlex | [10.5281/zenodo.21432935](https://doi.org/10.5281/zenodo.21432935) | [10.5281/zenodo.21433011](https://doi.org/10.5281/zenodo.21433011) |
| Perturb-seq | [10.5281/zenodo.21432937](https://doi.org/10.5281/zenodo.21432937) | [10.5281/zenodo.21433138](https://doi.org/10.5281/zenodo.21433138) |
| Tahoe | [10.5281/zenodo.21433000](https://doi.org/10.5281/zenodo.21433000) | [10.5281/zenodo.21433050](https://doi.org/10.5281/zenodo.21433050) |
