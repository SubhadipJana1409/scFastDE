# scFastDE 0.99.0

## New features

* Initial release submitted to Bioconductor.
* `filterSparseDonors()`: remove donors below a minimum cell count threshold
  per cell type before pseudo-bulk aggregation.
* `fastPseudobulk()`: build donor-weighted pseudo-bulk profiles from a
  SingleCellExperiment using vectorised sparse matrix operations.
* `fastDE()`: run vectorised differential expression across all genes
  simultaneously using limma-voom with donor cell-count weights.
* `plotDEResults()`: volcano plot with significance colouring and
  top-gene labelling.
* `FDEResult` S4 class for structured storage of DE results.
