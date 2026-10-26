# Step 2 - build the dentate Seurat object and cluster.
#
# Reads the 10x-style triplet + obs table from step 1 and runs:
#   SCTransform (regress percent.mt) -> RunPCA (50 PCs)
#   -> RunUMAP (dims 1:15) -> FindNeighbors (dims 1:10) -> FindClusters (res 0.8)
#   -> NormalizeData on RNA (LogNormalize, 10000)  # populates the RNA data slot
#
# Saves work/dentate_clustered.RDS for the annotation step.

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

set.seed(42)

DATA_DIR <- "data"
WORK_DIR <- "work"
dir.create(WORK_DIR, showWarnings = FALSE)

log <- function(msg) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg))

log("Reading counts (MTX triplet)")
counts <- ReadMtx(
  mtx      = file.path(DATA_DIR, "counts.mtx.gz"),
  cells    = file.path(DATA_DIR, "barcodes.tsv.gz"),
  features = file.path(DATA_DIR, "features.tsv.gz"),
  cell.column = 1, feature.column = 1
)

log("Reading obs.csv.gz")
obs <- read.csv(gzfile(file.path(DATA_DIR, "obs.csv.gz")),
                stringsAsFactors = FALSE, check.names = FALSE)
rownames(obs) <- obs$cell_id
obs <- obs[colnames(counts), , drop = FALSE]

log("Creating Seurat object")
obj <- CreateSeuratObject(counts = counts, meta.data = obs,
                          project = "mouse_dentate")
rm(counts, obs); gc(verbose = FALSE)

log("percent.mt")
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")

log("SCTransform (vars.to.regress = percent.mt)")
obj <- SCTransform(obj, vars.to.regress = "percent.mt", verbose = FALSE)

log("RunPCA (npcs = 50)")
obj <- RunPCA(obj, npcs = 50, verbose = FALSE)

log("RunUMAP (dims = 1:15, n.neighbors = 30, min.dist = 0.3)")
obj <- RunUMAP(obj, dims = 1:15, n.neighbors = 30, min.dist = 0.3, verbose = FALSE)

log("FindNeighbors (dims = 1:10, k.param = 20)")
obj <- FindNeighbors(obj, dims = 1:10, k.param = 20, verbose = FALSE)

log("FindClusters (resolution = 0.8)")
obj <- FindClusters(obj, resolution = 0.8, verbose = FALSE)

log("NormalizeData on RNA (LogNormalize, scale.factor = 10000)")
DefaultAssay(obj) <- "RNA"
obj <- NormalizeData(obj, normalization.method = "LogNormalize",
                     scale.factor = 10000, verbose = FALSE)

log("Saving work/dentate_clustered.RDS")
saveRDS(obj, file.path(WORK_DIR, "dentate_clustered.RDS"), compress = TRUE)
log("=== DONE ===")
print(obj)
