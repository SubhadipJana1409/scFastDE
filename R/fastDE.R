#' @title Fast Vectorised Differential Expression
#'
#' @description
#' Runs differential expression analysis across all genes simultaneously
#' using a weighted limma-voom linear model. Unlike tools that loop over
#' genes serially, \code{fastDE} passes the entire pseudo-bulk matrix to
#' limma, which uses LAPACK routines to fit all gene models in one
#' vectorised call. Donor weights from \code{\link{fastPseudobulk}} are
#' incorporated directly into the linear model via \code{limma::lmFit}.
#'
#' @details
#' The analysis pipeline is:
#' \enumerate{
#'   \item Filter lowly-expressed genes (CPM > \code{min_cpm} in at least
#'     \code{min_donors} donors).
#'   \item Apply \code{limma::voom} with donor cell-count weights.
#'   \item Fit a weighted linear model: \code{~ 0 + condition}.
#'   \item Extract the contrast of interest with \code{limma::makeContrasts}.
#'   \item Apply empirical Bayes moderation with \code{limma::eBayes}.
#'   \item Return all results as a \code{\link{FDEResult}} object.
#' }
#'
#' The weighting scheme means donors with more cells have proportionally
#' more influence on the log-fold change estimate, while donors with few
#' cells still contribute rather than being discarded.
#'
#' @param sce A \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#'   object with raw counts in \code{assay(sce, "counts")}.
#' @param donor A \code{character(1)} naming the donor column in
#'   \code{colData(sce)}.
#' @param cell_type A \code{character(1)} naming the cell type column in
#'   \code{colData(sce)}.
#' @param condition A \code{character(1)} naming the condition column in
#'   \code{colData(sce)} (e.g. \code{"treatment"}, \code{"disease"}).
#'   Must have exactly two unique values.
#' @param target_type A \code{character(1)} specifying which cell type
#'   to test. Must be present in \code{colData(sce)[[cell_type]]}.
#' @param contrast A \code{character(1)} specifying the contrast in
#'   limma syntax, e.g. \code{"treat - ctrl"}. If \code{NULL}, the
#'   second level minus the first level is used. Default: \code{NULL}.
#' @param min_cpm A \code{numeric(1)} minimum CPM for a gene to be
#'   considered expressed. Default: \code{1}.
#' @param min_donors An \code{integer(1)} minimum number of donors in
#'   which a gene must exceed \code{min_cpm}. Default: \code{2}.
#' @param min_cells An \code{integer(1)} minimum cells per donor before
#'   pseudo-bulk aggregation. Passed to \code{\link{filterSparseDonors}}.
#'   Default: \code{10}.
#' @param BPPARAM A \code{\link[BiocParallel]{BiocParallelParam}} object.
#'   Default: \code{SerialParam()}.
#'
#' @return A \code{\link{FDEResult}} object containing:
#'   \itemize{
#'     \item \code{deTable}: per-gene statistics (\code{logFC},
#'       \code{AveExpr}, \code{t}, \code{P.Value}, \code{adj.P.Val},
#'       \code{B}).
#'     \item \code{pseudobulk}: the aggregated count matrix.
#'     \item \code{donorWeights}: the per-donor \code{sqrt(n_cells)}
#'       weights.
#'     \item \code{params}: the analysis parameters.
#'   }
#'
#' @examples
#' library(SingleCellExperiment)
#'
#' set.seed(42)
#' n_genes <- 200
#' counts <- matrix(rpois(n_genes * 60, 8), nrow = n_genes, ncol = 60)
#' rownames(counts) <- paste0("Gene", seq_len(n_genes))
#' colnames(counts) <- paste0("Cell", seq_len(60))
#'
#' # Inject DE signal into first 10 genes for treatment group
#' counts[1:10, 31:60] <- counts[1:10, 31:60] * 3L
#'
#' sce <- SingleCellExperiment(assays = list(counts = counts))
#' sce$donor     <- rep(paste0("D", 1:6), each = 10)
#' sce$cell_type <- "Tcell"
#' sce$condition <- rep(c("ctrl", "treat"), each = 30)
#'
#' result <- fastDE(sce,
#'                  donor       = "donor",
#'                  cell_type   = "cell_type",
#'                  condition   = "condition",
#'                  target_type = "Tcell",
#'                  min_cells   = 5)
#' result
#' head(deTable(result))
#'
#' @seealso \code{\link{fastPseudobulk}}, \code{\link{filterSparseDonors}},
#'   \code{\link{plotDEResults}}
#'
#' @importFrom limma voom lmFit makeContrasts contrasts.fit eBayes
#'   topTable
#' @importFrom SummarizedExperiment colData
#' @importFrom S4Vectors DataFrame
#' @importFrom BiocParallel SerialParam
#' @importFrom methods is
#' @importFrom stats model.matrix
#' @export
fastDE <- function(sce,
                    donor       = NULL,
                    cell_type   = NULL,
                    condition   = NULL,
                    target_type = NULL,
                    contrast    = NULL,
                    min_cpm     = 1,
                    min_donors  = 2L,
                    min_cells   = 10L,
                    BPPARAM     = SerialParam()) {
    # ── Validation ────────────────────────────────────────────────────────────
    if (!is(sce, "SingleCellExperiment"))
        stop("'sce' must be a SingleCellExperiment object.")
    for (col in c(donor, cell_type, condition)) {
        if (is.null(col) || !col %in% names(colData(sce)))
            stop("Column '", col, "' not found in colData(sce).")
    }
    if (is.null(target_type))
        stop("'target_type' must be specified.")

    # ── Step 1: filter sparse donors ─────────────────────────────────────────
    sce <- filterSparseDonors(sce,
                               donor     = donor,
                               cell_type = cell_type,
                               min_cells = min_cells,
                               action    = "remove")

    # ── Step 2: build pseudo-bulk ─────────────────────────────────────────────
    pb <- fastPseudobulk(sce,
                          donor       = donor,
                          cell_type   = cell_type,
                          target_type = target_type)

    pb_mat   <- pb$pseudobulk
    weights  <- pb$donor_weights
    n_donors <- ncol(pb_mat)

    # ── Step 3: build design matrix ───────────────────────────────────────────
    # Extract per-donor condition (take the most common condition label)
    sce_sub   <- sce[, as.character(colData(sce)[[cell_type]]) == target_type]
    donor_sub <- as.character(colData(sce_sub)[[donor]])
    cond_sub  <- as.character(colData(sce_sub)[[condition]])

    donor_condition <- tapply(cond_sub, donor_sub, function(x) {
        names(sort(table(x), decreasing = TRUE))[1]
    })
    donor_condition <- donor_condition[colnames(pb_mat)]

    cond_levels <- sort(unique(donor_condition))
    if (length(cond_levels) != 2)
        stop("'condition' must have exactly 2 levels. Found: ",
             paste(cond_levels, collapse = ", "))

    donor_condition <- factor(donor_condition, levels = cond_levels)
    design <- model.matrix(~ 0 + donor_condition)
    colnames(design) <- cond_levels

    # ── Step 4: filter lowly expressed genes ─────────────────────────────────
    lib_sizes <- colSums(pb_mat)
    cpm_mat   <- sweep(pb_mat, 2, lib_sizes / 1e6, FUN = "/")
    keep      <- rowSums(cpm_mat >= min_cpm) >= min_donors

    if (sum(keep) == 0)
        stop("No genes passed the expression filter. ",
             "Try lowering 'min_cpm' or 'min_donors'.")

    pb_filt <- pb_mat[keep, , drop = FALSE]
    message(sprintf("Testing %d / %d genes after expression filter.",
                    sum(keep), nrow(pb_mat)))

    # ── Step 5: voom + weighted linear model ─────────────────────────────────
    v <- voom(pb_filt, design,
              weights  = matrix(rep(weights, each = nrow(pb_filt)),
                                nrow = nrow(pb_filt)),
              plot     = FALSE)

    fit <- lmFit(v, design, weights = v$weights)

    # ── Step 6: contrast ──────────────────────────────────────────────────────
    if (is.null(contrast))
        contrast <- paste(cond_levels[2], "-", cond_levels[1])

    contrast_mat <- makeContrasts(contrasts = contrast, levels = design)
    fit2         <- contrasts.fit(fit, contrast_mat)
    fit2         <- eBayes(fit2)

    # ── Step 7: collect results ───────────────────────────────────────────────
    tt <- topTable(fit2, number = Inf, sort.by = "none")
    de_df <- DataFrame(
        gene      = rownames(tt),
        logFC     = tt[["logFC"]],
        AveExpr   = tt[["AveExpr"]],
        t         = tt[["t"]],
        P.Value   = tt[["P.Value"]],
        adj.P.Val = tt[["adj.P.Val"]],
        B         = tt[["B"]],
        row.names = rownames(tt)
    )

    FDEResult(
        deTable      = de_df,
        pseudobulk   = pb_filt,
        donorWeights = weights,
        params       = list(
            cell_type   = target_type,
            condition   = condition,
            donor       = donor,
            contrast    = contrast,
            min_cpm     = min_cpm,
            min_donors  = min_donors,
            min_cells   = min_cells
        )
    )
}
