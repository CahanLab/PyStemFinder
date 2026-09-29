# Export R stemFinder reference results for the Python parity tests.
#
# Reproduces the R vignette (https://github.com/CahanLab/stemfinder) on the
# Tabula Muris bone marrow data and writes the inputs (marker expression, kNN
# graph) and outputs (scores for each method) as CSV files. Then run
# make_bmmc_fixture.py to pack them into bmmc_r_reference.h5ad.
#
# One deliberate difference from the vignette: E2f8, which is on both the S and
# the G2M list, is used once (unique()), as in PyStemFinder. With the vignette's
# list, which counts it twice, R's scores are identical to the published
# https://cnobjects.s3.amazonaws.com/stemFinder/bmmc_sF_results.csv.
#
# Usage:
#   Rscript export_bmmc_reference.R <stemfinder R repo> <MurineBoneMarrow10X_GSE109774.rds> <out dir>
#
# The .rds is https://cnobjects.s3.amazonaws.com/stemFinder/MurineBoneMarrow10X_GSE109774.rds

args <- commandArgs(trailingOnly = TRUE)
repo <- args[1]; rds <- args[2]; out <- args[3]
dir.create(out, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})
for (f in c("run_stemFinder.R", "compute_performance_single.R", "auc_probability.R", "pct_recover.R")) {
  source(file.path(repo, "R", f))
}
for (f in c("s_genes_mouse.rda", "g2m_genes_mouse.rda")) load(file.path(repo, "data", f))

adata <- readRDS(rds)

# --- vignette steps -----------------------------------------------------------
cell_cycle_genes <- unique(c(s_genes_mouse, g2m_genes_mouse)[c(s_genes_mouse, g2m_genes_mouse) %in% rownames(adata)])
VariableFeatures(adata) <- VariableFeatures(adata)[!(VariableFeatures(adata) %in% cell_cycle_genes)]
adata <- RunPCA(adata, verbose = FALSE)
pcs <- 32
k <- round(sqrt(ncol(adata)))
adata <- FindNeighbors(adata, dims = 1:pcs, k.param = k, verbose = FALSE)
knn <- adata@graphs$RNA_nn

scores <- data.frame(row.names = colnames(adata))
for (m in c("gini", "stdev", "variance")) {
  res <- run_stemFinder(adata, k = k, nn = knn, thresh = 0, markers = cell_cycle_genes, method = m)
  scores[[paste0(m, "_stemFinder_raw")]] <- res$stemFinder_raw
  scores[[paste0(m, "_stemFinder")]] <- res$stemFinder
}
adata$stemFinder_raw <- scores$gini_stemFinder_raw
adata$stemFinder <- scores$gini_stemFinder

set.seed(1)
perf <- compute_performance_single(adata, competitor = FALSE)
pct <- pct_recover(adata)

# --- export -------------------------------------------------------------------
cells <- colnames(adata)
stopifnot(identical(colnames(knn), cells), identical(rownames(knn), cells))
meta <- data.frame(cell = cells, Phenotype = adata$Phenotype, Ground_truth = adata$Ground_truth, scores)
write.csv(meta, file.path(out, "obs.csv"), row.names = FALSE)
writeLines(cell_cycle_genes, file.path(out, "markers.txt"))
write.csv(t(as.matrix(adata@assays$RNA@scale.data[cell_cycle_genes, cells])), file.path(out, "scale_data.csv"))
write.csv(t(as.matrix(adata@assays$RNA@data[cell_cycle_genes, cells])), file.path(out, "data.csv"))
nn <- summary(as(knn, "CsparseMatrix"))
write.csv(data.frame(i = nn$i - 1, j = nn$j - 1), file.path(out, "knn.csv"), row.names = FALSE)
write.csv(
  data.frame(
    metric = c(names(perf[["stemFinder results"]]), "pct_recover"),
    value = c(unname(perf[["stemFinder results"]]), pct)
  ),
  file.path(out, "performance.csv"), row.names = FALSE
)
writeLines(
  c(paste("k", k), paste("pcs", pcs), paste("R", getRversion()),
    paste("Seurat", packageVersion("Seurat")), paste("SeuratObject", packageVersion("SeuratObject"))),
  file.path(out, "session.txt")
)
