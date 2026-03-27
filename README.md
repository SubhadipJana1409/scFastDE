# scFastDE

[![Bioconductor devel](https://img.shields.io/badge/Bioconductor-devel-blue)](https://bioconductor.org/packages/devel/bioc/html/scFastDE.html)
[![R CMD check](https://github.com/SubhadipJana1409/scFastDE/actions/workflows/bioc-check.yml/badge.svg)](https://github.com/SubhadipJana1409/scFastDE/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

**scFastDE** provides fast, donor-weighted pseudo-bulk differential expression
for multi-donor single-cell RNA-seq experiments.

Standard pseudo-bulk DE tools (DESeq2, edgeR, muscat) loop over genes serially
and treat all donors equally regardless of cell count. scFastDE fixes this:

- **10-50x faster** — all genes tested simultaneously via vectorised matrix ops
- **Statistically principled** — samples weighted by `sqrt(n_cells)`
- **Paired design auto-detection** — same donors in multiple conditions? Automatically handled with a blocking model
- **Rare cell type aware** — sparse donor guard before aggregation
- **Bioc-native** — `SingleCellExperiment` in, `FDEResult` S4 out

## Installation

```r
BiocManager::install("scFastDE")

# Development version
BiocManager::install("SubhadipJana1409/scFastDE")
```

## Quick start

```r
library(scFastDE)

# 1. Remove donors with too few cells
sce <- filterSparseDonors(sce, donor = "donor",
                           cell_type = "cell_type", min_cells = 10)

# 2. Run DE (auto-builds pseudo-bulk + weights internally)
result <- fastDE(sce,
                  donor       = "donor",
                  cell_type   = "cell_type",
                  condition   = "disease",
                  target_type = "CD4_Tcell")

# 3. Inspect results
head(deTable(result))

# 4. Visualise
plotDEResults(result)
```

## Key functions

| Function | Description |
|---|---|
| `filterSparseDonors()` | Remove donors below min cell threshold |
| `fastPseudobulk()` | Vectorised weighted pseudo-bulk aggregation (per-donor or per-donor × condition) |
| `fastDE()` | Fast vectorised limma-voom DE with auto-detected paired/unpaired design |
| `plotDEResults()` | Volcano plot with significance colouring |

## Citation

> Subhadip Jana (2025). scFastDE: Fast Donor-Weighted Pseudo-Bulk DE for scRNA-seq.
> R package version 0.99.0.

## License

MIT © Subhadip Jana
