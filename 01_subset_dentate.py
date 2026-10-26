"""Step 1 - subset the Langlieb et al. 2023 atlas to the dentate gyrus.

Selects dentate gyrus libraries from the full single-nucleus atlas
(Library_Metadata.tsv sub_region in {DG1, DG2, DG3}) and exports the raw
counts of the matching nuclei as a 10x-style triplet plus a per-nucleus
metadata table, ready for the Seurat build in step 2.

Inputs, read from $ATLAS_DIR (default: current directory):
  Macosko_Mouse_Atlas_Single_Nuclei.Use_Backed.h5ad   full atlas (backed read)
  Library_Metadata.tsv                                 per-library metadata
  CellType_Metadata.tsv                                per-cluster metadata

Outputs (./data/):
  counts.mtx.gz    sparse genes x nuclei, integer counts
  features.tsv.gz  one gene symbol per line
  barcodes.tsv.gz  one nucleus barcode per line
  obs.csv.gz       per-nucleus metadata (atlas obs + library + cluster joins)
"""

import gzip
import os
import shutil
import sys
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp

# ---- paths ----
ATLAS_DIR = Path(os.environ.get("ATLAS_DIR", "."))
H5AD = ATLAS_DIR / "Macosko_Mouse_Atlas_Single_Nuclei.Use_Backed.h5ad"
LIB_TSV = ATLAS_DIR / "Library_Metadata.tsv"
CELL_TSV = ATLAS_DIR / "CellType_Metadata.tsv"
OUT = Path("data")
DG_SUB_REGIONS = {"DG1", "DG2", "DG3"}


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def main():
    OUT.mkdir(parents=True, exist_ok=True)

    log("Reading Library_Metadata.tsv")
    lib = pd.read_csv(LIB_TSV, sep="\t")
    dg_libs_df = lib.loc[lib["sub_region"].isin(DG_SUB_REGIONS)]
    dg_libs = set(dg_libs_df["library"])
    log(f"  dentate libraries: {len(dg_libs)} ({sorted(dg_libs)})")

    log("Reading CellType_Metadata.tsv")
    cell_meta = pd.read_csv(CELL_TSV, sep="\t")

    log("Opening atlas h5ad (backed)")
    adata = ad.read_h5ad(H5AD, backed="r")
    log(f"  full atlas: {adata.n_obs:,} nuclei x {adata.n_vars:,} genes")

    log("Selecting dentate nuclei (obs.derived_cell_libs in dentate libraries)")
    libs_col = adata.obs["derived_cell_libs"].astype(str).values
    mask = np.isin(libs_col, list(dg_libs))
    n_target = int(mask.sum())
    log(f"  matched {n_target:,} nuclei")
    if n_target == 0:
        sys.exit("ERROR: zero nuclei matched. Check derived_cell_libs format.")

    log("Materializing subset (loads counts for the selected nuclei)")
    subset = adata[mask, :].to_memory()

    log("Transposing to genes x nuclei")
    X = subset.X
    if not sp.issparse(X):
        X = sp.csr_matrix(X)
    counts = X.T.tocsr()

    log("Resolving gene symbols")
    # Keep every atlas feature (all 21,899). A handful of gene symbols are
    # duplicated; they are made unique downstream by Seurat's ReadMtx
    # (unique.features = TRUE), matching the published object. No collapsing.
    gene_names = subset.var["gene_name"].astype(str).values
    n_genes, n_unique = len(gene_names), len(set(gene_names))
    log(f"  genes={n_genes}, unique symbols={n_unique}, duplicates={n_genes - n_unique}")
    feature_ids = list(gene_names)

    if not np.issubdtype(counts.dtype, np.integer):
        counts.data = counts.data.astype(np.int32)

    log("Writing counts.mtx.gz")
    mtx_tmp = OUT / "counts.mtx"
    scipy.io.mmwrite(str(mtx_tmp), counts, field="integer")
    with open(mtx_tmp, "rb") as fi, gzip.open(OUT / "counts.mtx.gz", "wb", 4) as fo:
        shutil.copyfileobj(fi, fo, length=8 * 1024 * 1024)
    os.remove(mtx_tmp)

    log("Writing features.tsv.gz")
    with gzip.open(OUT / "features.tsv.gz", "wt") as f:
        f.write("\n".join(feature_ids) + "\n")

    log("Writing barcodes.tsv.gz")
    with gzip.open(OUT / "barcodes.tsv.gz", "wt") as f:
        f.write("\n".join(subset.obs.index) + "\n")

    log("Writing obs.csv.gz (atlas obs + library + cluster metadata)")
    obs = subset.obs.copy()
    obs["library"] = obs["derived_cell_libs"].astype(str)
    obs = obs.reset_index().rename(columns={obs.index.name or "index": "cell_id"})
    if "cell_id" not in obs.columns:
        obs = obs.rename(columns={obs.columns[0]: "cell_id"})
    obs = obs.merge(lib, on="library", how="left")
    obs = obs.merge(cell_meta, left_on="ClusterKey", right_on="key",
                    how="left", suffixes=("", "_clust"))
    obs.to_csv(OUT / "obs.csv.gz", compression="gzip", index=False)

    log("=== DONE ===")
    log(f"  nuclei={n_target:,}  genes={len(feature_ids):,}  -> {OUT}/")


if __name__ == "__main__":
    main()
