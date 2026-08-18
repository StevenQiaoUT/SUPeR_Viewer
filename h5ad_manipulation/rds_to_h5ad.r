library(Seurat)
library(Matrix)
library(anndata)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript rds_to_h5ad.R <input.rds> <output.h5ad> [assay] [reduction] [celltype_col]")
}
rds_path      <- args[1]
h5ad_path     <- args[2]
assay_name    <- if (length(args) >= 3) args[3] else "RNA"
reduction_name<- if (length(args) >= 4) args[4] else "umap"
celltype_col  <- if (length(args) >= 5) args[5] else "CellType"

message("Loading: ", rds_path)
obj <- readRDS(rds_path)

message("Extracting raw counts from assay '", assay_name, "'")
counts <- obj[[assay_name]]$counts
counts_t <- Matrix::t(counts)
counts_t <- as(counts_t, "CsparseMatrix")

gene_ids <- rownames(counts)
cell_ids <- colnames(counts)
rm(counts); gc()

lib_sizes <- Matrix::rowSums(counts_t)
counts_t@x <- counts_t@x / lib_sizes[counts_t@i + 1L] * 1e4
counts_t@x <- log1p(counts_t@x)

if (!celltype_col %in% colnames(obj@meta.data)) {
  stop("Column '", celltype_col, "' not found in obj@meta.data. Available: ",
       paste(colnames(obj@meta.data), collapse = ", "))
}
obs_df <- data.frame(
  CellType = as.character(obj@meta.data[[celltype_col]]),
  row.names = cell_ids,
  stringsAsFactors = FALSE
)
var_df <- data.frame(row.names = gene_ids)

message("Extracting reduction '", reduction_name, "'")
umap_mat <- Embeddings(obj, reduction = reduction_name) 
umap_mat <- umap_mat[cell_ids, , drop = FALSE]        
colnames(umap_mat) <- c("UMAP1", "UMAP2")

rm(obj); gc()

message("Building AnnData (", nrow(counts_t), " cells x ", ncol(counts_t), " genes)")
ad_py <- reticulate::import("anndata")
ad <- ad_py$AnnData(
  X = counts_t,
  obs = obs_df,
  var = var_df,
  obsm = list(X_umap = umap_mat)
)

message("Writing: ", h5ad_path)
ad$write_h5ad(h5ad_path)
message("Done.")