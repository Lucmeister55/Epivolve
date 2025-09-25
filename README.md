# MethylPipe

**MethylPipe** is an R package for downstream processing and analysis of cfRRBS (cell-free Reduced Representation Bisulfite Sequencing) data. It provides functions for:

- Loading and filtering coverage data
- Sample and site quality control
- Aggregating methylation data to regions
- Dimensionality reduction (PCA, UMAP, t-SNE)
- Batch effect detection and correction
- Generating publication-ready plots

---

## Installation

### From GitHub (current)

```R
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  "methylKit",
  "GenomicRanges",
  "IRanges",
  "SummarizedExperiment",
  "impute"
))

install.packages("remotes")
remotes::install_github("Lucmeister55/MethylPipe", dependencies = TRUE)

library(MethylPipe)
```