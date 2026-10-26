# Step 3 - annotate clusters and save the final dentate object.
#
# Clusters from step 2 were assigned to cell types by canonical marker
# expression, with hippocampal principal-cell markers taken from Hipposeq
# (Cembrowski MS et al. Hipposeq: a comprehensive RNA-seq database of gene
# expression in hippocampal principal neurons. eLife 5:e14997 (2016).
# https://doi.org/10.7554/eLife.14997), and low-quality / doublet clusters
# were dropped. To make the
# published object reproducible regardless of clustering stochasticity across
# Seurat versions, the curated per-nucleus assignment is provided as a barcode
# table (dentate_cell_annotations.csv) and joined here by barcode. Nuclei absent
# from the table are the dropped clusters and are excluded.
#
# Writes mouse_dentate.RDS (+ .yaml) for the browser / figure code.

suppressPackageStartupMessages({
  library(Seurat)
  library(yaml)
})

WORK_DIR <- "work"
ANNO_CSV <- "dentate_cell_annotations.csv"
OUT_RDS  <- "mouse_dentate.RDS"
OUT_YAML <- "mouse_dentate.yaml"

log <- function(msg) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg))

log("Loading clustered object")
obj <- readRDS(file.path(WORK_DIR, "dentate_clustered.RDS"))

log("Reading curated annotations")
anno <- read.csv(ANNO_CSV, stringsAsFactors = FALSE)
ct <- setNames(anno$cell_type, anno$barcode)

log("Restricting to annotated nuclei (drops QC/doublet clusters)")
keep <- intersect(colnames(obj), names(ct))
log(sprintf("  keeping %d / %d nuclei", length(keep), ncol(obj)))
obj <- subset(obj, cells = keep)

# preserved cell-type order of the published object
lvls <- c("CA3_pyr", "GABAergic Neurons", "Fibroblasts", "CA1_pyr",
          "Granule Cells", "Mossy_Cells", "Microglia", "Neuroblasts",
          "Oligodendrocytes", "OPC", "Astrocytes", "rNSC", "Endothelial Cells")
obj$cell_type <- factor(unname(ct[colnames(obj)]), levels = lvls)
Idents(obj) <- "cell_type"

log("Finalizing (RNA assay only)")
DefaultAssay(obj) <- "RNA"
if ("SCT" %in% names(obj@assays)) obj[["SCT"]] <- NULL

log(sprintf("Saving %s", OUT_RDS))
saveRDS(obj, OUT_RDS, compress = TRUE)

write_yaml(list(
  name     = "Mouse Brain Atlas - Dentate Gyrus",
  citation = paste("Langlieb, J., Sachdev, N.S., Balderrama, K.S. et al.",
                   "The molecular cytoarchitecture of the adult mouse brain.",
                   "Nature 624, 333-342 (2023)."),
  doi      = "https://doi.org/10.1038/s41586-023-06818-7"
), OUT_YAML)

log("=== DONE ===")
print(table(Idents(obj)))
