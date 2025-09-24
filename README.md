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

You will need the `remotes` package:

```R
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.21")

install.packages("remotes")
remotes::install_github("Lucmeister55/MethylPipe", dependencies = TRUE)

library(MethylPipe)
```