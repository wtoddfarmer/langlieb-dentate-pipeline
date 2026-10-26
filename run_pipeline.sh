#!/usr/bin/env bash
# Run the full pipeline: atlas -> mouse_dentate.RDS.
#
# Usage: ./run_pipeline.sh ATLAS_DIR
#
# ATLAS_DIR holds Macosko_Mouse_Atlas_Single_Nuclei.Use_Backed.h5ad,
# Library_Metadata.tsv and CellType_Metadata.tsv. The Python environment is
# created in .venv/ from requirements.txt on first run. Each step's output is
# also written to logs/.

set -euo pipefail

if [ $# -ne 1 ]; then
  echo "usage: $0 ATLAS_DIR" >&2
  exit 1
fi

ATLAS_DIR=$(cd "$1" && pwd)
export ATLAS_DIR
cd "$(dirname "$0")"

for f in Macosko_Mouse_Atlas_Single_Nuclei.Use_Backed.h5ad \
         Library_Metadata.tsv CellType_Metadata.tsv; do
  if [ ! -f "$ATLAS_DIR/$f" ]; then
    echo "missing input: $ATLAS_DIR/$f" >&2
    exit 1
  fi
done
command -v Rscript >/dev/null || { echo "Rscript not found" >&2; exit 1; }

if [ ! -x .venv/bin/python ]; then
  echo "== creating .venv"
  python3 -m venv .venv
  ./.venv/bin/pip install --quiet -r requirements.txt
fi

mkdir -p logs

echo "== step 1: subset atlas to dentate gyrus"
./.venv/bin/python 01_subset_dentate.py 2>&1 | tee logs/01_subset_dentate.log

echo "== step 2: build and cluster Seurat object"
Rscript 02_build_dentate_seurat.R 2>&1 | tee logs/02_build_dentate_seurat.log

echo "== step 3: annotate and save"
Rscript 03_annotate_and_save.R 2>&1 | tee logs/03_annotate_and_save.log

echo "== done: mouse_dentate.RDS, mouse_dentate.yaml"
