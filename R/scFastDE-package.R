#' scFastDE: Fast Donor-Weighted Pseudo-Bulk DE for scRNA-seq
#'
#' @description
#' scFastDE provides fast, donor-weighted pseudo-bulk differential expression
#' analysis for multi-donor single-cell RNA-seq experiments. Unlike existing
#' tools that test genes serially in a loop, scFastDE uses vectorised sparse
#' matrix operations across all genes simultaneously, giving 10-50x speed
#' gains on large datasets (30k+ genes, 100k+ cells).
#'
#' @section Key innovations:
#' \describe{
#'   \item{Vectorised gene testing}{All genes are tested simultaneously via
#'     matrix algebra — no per-gene loop.}
#'   \item{Donor cell-count weighting}{Each sample's pseudo-bulk profile is
#'     weighted by \code{sqrt(n_cells)}, giving principled influence to
#'     well-represented samples.}
#'   \item{Paired design auto-detection}{When the same donors appear in
#'     multiple conditions, \code{fastDE} automatically aggregates per
#'     donor-condition pair and uses a blocking model
#'     (\code{~ 0 + condition + donor}) to account for inter-donor
#'     variation.}
#'   \item{Sparse pseudo-bulk guard}{Donors with fewer than
#'     \code{min_cells} cells per cell type are flagged and optionally
#'     removed before aggregation.}
#' }
#'
#' @section Main functions:
#' \describe{
#'   \item{\code{\link{filterSparseDonors}}}{Remove donors below minimum
#'     cell threshold per cell type.}
#'   \item{\code{\link{fastPseudobulk}}}{Build donor-weighted pseudo-bulk
#'     matrix from a SingleCellExperiment. Supports per-donor or
#'     per-donor-condition aggregation.}
#'   \item{\code{\link{fastDE}}}{Run vectorised DE across all genes.
#'     Auto-detects paired vs unpaired designs.}
#'   \item{\code{\link{plotDEResults}}}{Volcano plot of DE results.}
#' }
#'
#' @references
#' Crowell HL et al. (2020). muscat detects subpopulation-specific state
#' transitions from multi-sample multi-condition single-cell transcriptomics
#' data. \emph{Nature Communications}, 11, 6077.
#'
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importFrom rlang .data
#' @docType package
#' @name scFastDE-package
#' @aliases scFastDE
"_PACKAGE"
