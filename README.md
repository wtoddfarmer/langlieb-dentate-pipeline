# Dentate gyrus single-nucleus object

Pipeline that derives the dentate gyrus object (`mouse_dentate.RDS`) from the
adult mouse brain single-nucleus atlas of Langlieb et al. 2023.

> Langlieb, J., Sachdev, N.S., Balderrama, K.S. et al. The molecular
> cytoarchitecture of the adult mouse brain. *Nature* **624**, 333–342 (2023).
> https://doi.org/10.1038/s41586-023-06818-7

## Overview

Starting from the full atlas (4,407,296 nuclei), the six dentate gyrus libraries
(`sub_region` DG1/DG2/DG3) are selected, and their nuclei are re-processed and
annotated on their own. The result is `mouse_dentate.RDS`: 47,208 nuclei,
21,899 genes, 13 annotated cell types.

```
full atlas h5ad
      │  01_subset_dentate.py        select DG1/DG2/DG3 libraries
      ▼
10x triplet + obs                    counts.mtx.gz / features / barcodes / obs
      │  02_build_dentate_seurat.R   SCTransform -> PCA -> UMAP -> clustering
      ▼
clustered Seurat object
      │  03_annotate_and_save.R      marker-based cell types, drop QC clusters
      ▼
mouse_dentate.RDS
```

## Steps

| Script | Language | Output |
|---|---|---|
| `01_subset_dentate.py` | Python (anndata, scipy, pandas) | `data/` 10x triplet + `obs.csv.gz` |
| `02_build_dentate_seurat.R` | R (Seurat v5) | `work/dentate_clustered.RDS` |
| `03_annotate_and_save.R` | R (Seurat v5) | `mouse_dentate.RDS`, `mouse_dentate.yaml` |

`dentate_cell_annotations.csv` (barcode, cell_type) is the curated per-nucleus
annotation used by step 3.

## Processing parameters

- **Subset:** atlas nuclei whose library `sub_region` is DG1, DG2, or DG3.
- **Normalization / features:** `SCTransform`, regressing out `percent.mt`.
- **Dimensionality reduction:** `RunPCA` (50 PCs); `RunUMAP` on PCs 1–15
  (`n.neighbors = 30`, `min.dist = 0.3`).
- **Clustering:** `FindNeighbors` on PCs 1–10 (`k.param = 20`);
  `FindClusters` at resolution 0.8.
- **Expression slot for plotting:** `NormalizeData` on the RNA assay
  (LogNormalize, `scale.factor = 10000`). Figures and the browser read the RNA
  `data` slot; the SCT assay is dropped from the final object.

## Annotation and reproducibility

Clusters were assigned to cell types by canonical marker expression, with
hippocampal principal-cell markers taken from Hipposeq (Cembrowski MS et al.
Hipposeq: a comprehensive RNA-seq database of gene expression in hippocampal
principal neurons. *eLife* **5**, e14997 (2016).
https://doi.org/10.7554/eLife.14997); each cluster resolved to a single type. Four clusters were
excluded as low-quality / doublets. Because Louvain clustering is stochastic and
version-sensitive, cluster numbering is not guaranteed to be identical on a
re-run, so the curated per-nucleus labels are shipped in
`dentate_cell_annotations.csv` and joined by barcode in step 3. This reproduces
the published object exactly. Barcodes are the atlas nucleus identifiers, so
they match across the full atlas, this subset, and the annotation table.

The 13 cell types: Granule Cells, Mossy_Cells, CA1_pyr, CA3_pyr,
GABAergic Neurons, Neuroblasts, rNSC, Astrocytes, Oligodendrocytes, OPC,
Microglia, Endothelial Cells, Fibroblasts.

## Running

```bash
./run_pipeline.sh /path/to/atlas_dir
```

`atlas_dir` holds the atlas h5ad and the two metadata TSVs. The wrapper creates
the Python environment in `.venv/` on first run, runs the three steps in order,
stops at the first error, and writes each step's output to `logs/`.

The steps can also be run individually:

```bash
# one-time: create the Python environment
python3 -m venv .venv
./.venv/bin/pip install -r requirements.txt

# 1. subset (ATLAS_DIR defaults to the current directory)
ATLAS_DIR=/path/to/atlas_dir ./.venv/bin/python 01_subset_dentate.py

# 2. build + cluster
Rscript 02_build_dentate_seurat.R

# 3. annotate + save  ->  mouse_dentate.RDS
Rscript 03_annotate_and_save.R
```

Requirements: Python 3 with the packages pinned in `requirements.txt`
(`anndata`, `scipy`, `pandas`, `h5py`, `numpy`); R with `Seurat` (v5), `Matrix`,
`yaml`. The full atlas h5ad, `Library_Metadata.tsv`, and `CellType_Metadata.tsv`
are available from the Langlieb et al. 2023 data release. The `.venv/` and the
large intermediate outputs (`data/`, `work/`, the `.RDS` files) are git-ignored;
`requirements.txt` reproduces the environment.

## Verification

This pipeline was run end to end and its output compared to the published
`mouse_dentate.RDS`: identical dimensions (21,899 genes x 47,208 nuclei), 100%
cell-type agreement, and a bit-identical RNA `data` slot (maximum absolute
difference 0). Non-deterministic outputs (UMAP coordinates, cluster numbering)
do not affect the result because cell types are joined by barcode and figures
read the RNA `data` slot.
