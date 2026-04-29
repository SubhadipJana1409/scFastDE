# ============================================================================
# test_real_data.R — Test scFastDE with Kang et al. 2018 PBMC dataset
#
# Dataset: 8 lupus patients, stimulated (IFN-β) vs control PBMCs
# This is THE classic benchmark for pseudo-bulk DE methods.
#
# Install dependencies first (one-time):
#   BiocManager::install(c("muscData", "scFastDE", "scuttle"))
#
# Run:  source("test_real_data.R")
# ============================================================================

library(scFastDE)
library(SingleCellExperiment)

cat("========================================\n")
cat("  scFastDE Real Data Test\n")
cat("  Kang et al. 2018 — IFN-β PBMCs\n")
cat("========================================\n\n")

# ── 1. Load the Kang 2018 dataset ─────────────────────────────────────────
cat("Loading Kang et al. 2018 dataset from ExperimentHub...\n")
cat("(First run will download ~200 MB, then cached locally)\n\n")

if (!requireNamespace("muscData", quietly = TRUE)) {
    stop("Please install muscData first:\n  BiocManager::install('muscData')")
}

sce <- muscData::Kang18_8vs8()

cat(sprintf("Dataset loaded: %d genes × %d cells\n", nrow(sce), ncol(sce)))
cat(sprintf("Donors:     %s\n", paste(sort(unique(sce$ind)), collapse = ", ")))
cat(sprintf("Conditions: %s\n", paste(unique(sce$stim), collapse = ", ")))
cat(sprintf("Cell types: %s\n\n", paste(sort(unique(sce$cell)), collapse = ", ")))

# ── 2. Pick a cell type to test ───────────────────────────────────────────
# CD14+ Monocytes are the most abundant and have strong IFN-β response
target <- "CD14+ Monocytes"
cat(sprintf("Target cell type: %s\n", target))
n_target <- sum(sce$cell == target, na.rm = TRUE)
cat(sprintf("  %d cells of this type\n\n", n_target))

# ── 3. Rename columns to match scFastDE interface ────────────────────────
# muscData uses: ind (donor), stim (condition), cell (cell_type)
sce$donor     <- sce$ind
sce$cell_type <- sce$cell
sce$condition <- sce$stim

# Remove cells with NA cell type
sce <- sce[, !is.na(sce$cell_type)]

# ── 4. filterSparseDonors ────────────────────────────────────────────────
cat("Step 1: filterSparseDonors()...\n")
sce_filt <- filterSparseDonors(sce,
                                donor     = "donor",
                                cell_type = "cell_type",
                                min_cells = 10)
cat(sprintf("  Before: %d cells → After: %d cells\n\n",
            ncol(sce), ncol(sce_filt)))

# ── 5. fastPseudobulk ───────────────────────────────────────────────────
cat("Step 2: fastPseudobulk()...\n")
pb <- fastPseudobulk(sce_filt,
                      donor       = "donor",
                      cell_type   = "cell_type",
                      target_type = target)

cat(sprintf("  Pseudo-bulk: %d genes × %d donors\n",
            nrow(pb$pseudobulk), ncol(pb$pseudobulk)))
cat("  Donor weights (sqrt of cell counts):\n")
print(round(pb$donor_weights, 2))
cat(sprintf("  Cell counts per donor:\n"))
print(pb$donor_ncells)
cat("\n")

# ── 6. fastDE ────────────────────────────────────────────────────────────
cat("Step 3: fastDE()...\n")
t1 <- Sys.time()
result <- fastDE(sce_filt,
                  donor       = "donor",
                  cell_type   = "cell_type",
                  condition   = "condition",
                  target_type = target,
                  min_cells   = 10,
                  min_cpm     = 1,
                  min_donors  = 2)
t2 <- Sys.time()

cat(sprintf("  Completed in %.2f seconds\n", as.numeric(t2 - t1)))
show(result)
cat("\n")

# ── 7. DE results ────────────────────────────────────────────────────────
dt <- as.data.frame(deTable(result))
dt_sorted <- dt[order(dt$adj.P.Val), ]

n_sig_005 <- sum(dt$adj.P.Val < 0.05, na.rm = TRUE)
n_sig_001 <- sum(dt$adj.P.Val < 0.01, na.rm = TRUE)
n_up      <- sum(dt$adj.P.Val < 0.05 & dt$logFC > 1, na.rm = TRUE)
n_down    <- sum(dt$adj.P.Val < 0.05 & dt$logFC < -1, na.rm = TRUE)

cat("Step 4: DE Results Summary\n")
cat(sprintf("  Genes tested:            %d\n", nrow(dt)))
cat(sprintf("  Significant (FDR<0.05):  %d\n", n_sig_005))
cat(sprintf("  Significant (FDR<0.01):  %d\n", n_sig_001))
cat(sprintf("  Upregulated (|logFC|>1): %d\n", n_up))
cat(sprintf("  Downregulated:           %d\n\n", n_down))

# Known IFN-β response genes that SHOULD be detected
known_ifn_genes <- c("ISG15", "IFIT1", "IFIT3", "IFI6", "MX1",
                     "OAS1", "STAT1", "IRF7", "IFI44L", "ISG20")

cat("Known IFN-β response genes (should be significant):\n")
for (g in known_ifn_genes) {
    if (g %in% dt$gene) {
        row <- dt[dt$gene == g, ]
        sig <- ifelse(row$adj.P.Val < 0.05, "YES", " no")
        cat(sprintf("  %-8s  logFC=%+6.2f  FDR=%.2e  Sig=%s\n",
                    g, row$logFC, row$adj.P.Val, sig))
    } else {
        cat(sprintf("  %-8s  (filtered out / not found)\n", g))
    }
}
cat("\n")

cat("Top 20 DE genes:\n")
print(head(dt_sorted[, c("gene", "logFC", "P.Value", "adj.P.Val")], 20))
cat("\n")

# ── 8. Volcano plot ──────────────────────────────────────────────────────
cat("Step 5: plotDEResults() — Volcano plot...\n")
p <- plotDEResults(result, fdr_thresh = 0.05, lfc_thresh = 1, top_n = 15)
print(p)
cat("  Volcano plot displayed!\n\n")

# ── 9. Test a second cell type ───────────────────────────────────────────
target2 <- "CD4 T cells"
cat(sprintf("Bonus: Running fastDE on %s...\n", target2))
t1 <- Sys.time()
result2 <- fastDE(sce_filt,
                   donor       = "donor",
                   cell_type   = "cell_type",
                   condition   = "condition",
                   target_type = target2,
                   min_cells   = 10,
                   min_cpm     = 1,
                   min_donors  = 2)
t2 <- Sys.time()

dt2 <- as.data.frame(deTable(result2))
cat(sprintf("  %s: %d genes tested, %d significant (FDR<0.05) in %.2f sec\n\n",
            target2, nrow(dt2), sum(dt2$adj.P.Val < 0.05, na.rm = TRUE),
            as.numeric(t2 - t1)))

# ── Summary ──────────────────────────────────────────────────────────────
detected_ifn <- sum(known_ifn_genes %in%
                    dt$gene[dt$adj.P.Val < 0.05])

cat("========================================\n")
cat("  RESULTS SUMMARY\n")
cat("========================================\n")
cat(sprintf("  Known IFN-β genes detected: %d / %d\n",
            detected_ifn, length(known_ifn_genes)))
cat(sprintf("  Total DE genes (FDR<0.05): %d\n", n_sig_005))

if (detected_ifn >= 5) {
    cat("  ✓ PASS — scFastDE correctly detects IFN-β response!\n")
} else {
    cat("  ⚠ CHECK — Fewer known genes detected than expected.\n")
}
cat("========================================\n")
