#' @title Build Donor-Weighted Pseudo-Bulk Profiles
#'
#' @description
#' Aggregates single-cell counts into pseudo-bulk profiles (one per donor)
#' for a specified cell type, using vectorised sparse matrix operations.
#' Each donor's profile is weighted by \code{sqrt(n_cells)}, giving more
#' influence to donors with more cells while not completely discarding
#' donors with fewer cells.
#'
#' @details
#' Pseudo-bulk aggregation sums the raw counts for all cells belonging to
#' each donor within a given cell type:
#'
#' \deqn{PB_{g,d} = \sum_{c \in \text{donor}_d, \text{type}_t} X_{g,c}}
#'
#' Donor weights are then computed as:
#'
#' \deqn{w_d = \sqrt{n_{d,t}}}
#'
#' where \eqn{n_{d,t}} is the number of cells for donor \eqn{d} in cell
#' type \eqn{t}. These weights are passed to \code{\link{fastDE}} for use
#' in the weighted linear model.
#'
#' The aggregation uses sparse matrix column-sums grouped by donor,
#' avoiding a slow R-level loop over donors.
#'
#' @param sce A \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#'   object with raw counts in \code{assay(sce, "counts")}.
#' @param donor A \code{character(1)} naming the donor column in
#'   \code{colData(sce)}.
#' @param cell_type A \code{character(1)} naming the cell type column in
#'   \code{colData(sce)}.
#' @param target_type A \code{character(1)} specifying which cell type to
#'   aggregate. Must be a value present in \code{colData(sce)[[cell_type]]}.
#' @param assay_name A \code{character(1)} name of the assay to aggregate.
#'   Default: \code{"counts"}.
#'
#' @return A \code{list} with two elements:
#'   \describe{
#'     \item{\code{pseudobulk}}{A \code{matrix} of aggregated counts,
#'       rows = genes, columns = donors.}
#'     \item{\code{donor_weights}}{A named \code{numeric} vector of
#'       per-donor weights (\code{sqrt(n_cells)}).}
#'     \item{\code{donor_ncells}}{A named \code{integer} vector giving
#'       the raw cell count per donor.}
#'   }
#'
#' @examples
#' library(SingleCellExperiment)
#'
#' set.seed(42)
#' counts <- matrix(rpois(500 * 60, 5), nrow = 500, ncol = 60)
#' rownames(counts) <- paste0("Gene", seq_len(500))
#' colnames(counts) <- paste0("Cell", seq_len(60))
#' sce <- SingleCellExperiment(assays = list(counts = counts))
#' sce$donor     <- rep(paste0("D", 1:6), each = 10)
#' sce$cell_type <- "Tcell"
#' sce$condition <- rep(c("ctrl", "treat"), times = 30)
#'
#' pb <- fastPseudobulk(sce, donor = "donor",
#'                       cell_type = "cell_type",
#'                       target_type = "Tcell")
#' dim(pb$pseudobulk)
#' pb$donor_weights
#'
#' @seealso \code{\link{filterSparseDonors}}, \code{\link{fastDE}}
#'
#' @importFrom SummarizedExperiment assay colData
#' @importFrom Matrix sparseMatrix rowSums
#' @importFrom methods is
#' @importFrom stats setNames
#' @export
fastPseudobulk <- function(sce,
                            donor       = NULL,
                            cell_type   = NULL,
                            target_type = NULL,
                            assay_name  = "counts") {
    # ── Validation ────────────────────────────────────────────────────────────
    if (!is(sce, "SingleCellExperiment"))
        stop("'sce' must be a SingleCellExperiment object.")
    if (is.null(donor) || !donor %in% names(colData(sce)))
        stop("'donor' column '", donor, "' not found in colData(sce).")
    if (is.null(cell_type) || !cell_type %in% names(colData(sce)))
        stop("'cell_type' column '", cell_type, "' not found in colData(sce).")
    if (is.null(target_type))
        stop("'target_type' must be specified.")
    if (!assay_name %in% names(assays(sce)))
        stop("Assay '", assay_name, "' not found in sce.")

    ct_vals <- as.character(colData(sce)[[cell_type]])
    if (!target_type %in% ct_vals)
        stop("'target_type' = '", target_type, "' not found in ",
             "'", cell_type, "' column.")

    # ── Subset to target cell type ────────────────────────────────────────────
    sce_sub    <- sce[, ct_vals == target_type]
    donors_sub <- as.character(colData(sce_sub)[[donor]])
    donor_lvls <- unique(donors_sub)
    n_donors   <- length(donor_lvls)

    if (n_donors < 2)
        stop("At least 2 donors required for DE analysis. Found: ",
             n_donors, " donor(s) for cell type '", target_type, "'.")

    # ── Vectorised pseudo-bulk aggregation ────────────────────────────────────
    # Build a sparse indicator matrix: cells x donors
    # Then: pseudobulk = count_matrix %*% indicator
    counts_mat <- assay(sce_sub, assay_name)

    donor_idx <- match(donors_sub, donor_lvls)

    # Sparse indicator matrix (cells x donors)
    indicator <- sparseMatrix(
        i = seq_along(donor_idx),
        j = donor_idx,
        x = 1,
        dims = c(ncol(sce_sub), n_donors),
        dimnames = list(colnames(sce_sub), donor_lvls)
    )

    # Vectorised aggregation: genes x donors — one matrix multiply, no loop
    pb_matrix <- as.matrix(counts_mat %*% indicator)

    # ── Donor weights: sqrt(n_cells) ─────────────────────────────────────────
    n_cells_per_donor <- tabulate(donor_idx, nbins = n_donors)
    names(n_cells_per_donor) <- donor_lvls

    donor_weights <- sqrt(n_cells_per_donor)

    list(
        pseudobulk    = pb_matrix,
        donor_weights = donor_weights,
        donor_ncells  = n_cells_per_donor
    )
}


# ── Internal helper ───────────────────────────────────────────────────────────

#' Access assays from SCE (internal helper)
#' @keywords internal
#' @noRd
assays <- function(x) SummarizedExperiment::assays(x)
