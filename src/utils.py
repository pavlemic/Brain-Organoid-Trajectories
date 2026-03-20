"""
utils.py
--------
Small helper functions used across the project:
  - directory setup
  - data loading (h5ad and 10x format)
  - basic logging/printing
"""

import os
from pathlib import Path
import scanpy as sc
import anndata as ad


# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

# Project root is one level above src/
PROJECT_ROOT = Path(__file__).resolve().parent.parent

DATA_RAW     = PROJECT_ROOT / "data" / "raw"
DATA_PROC    = PROJECT_ROOT / "data" / "processed"
RESULTS_FIG  = PROJECT_ROOT / "results" / "figures"
RESULTS_TAB  = PROJECT_ROOT / "results" / "tables"


def ensure_dirs():
    """Create all output directories if they don't exist yet."""
    for d in [DATA_RAW, DATA_PROC, RESULTS_FIG, RESULTS_TAB]:
        d.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_h5ad(filename: str) -> ad.AnnData:
    """
    Load a pre-saved AnnData file (.h5ad format).

    Parameters
    ----------
    filename : str
        File name inside data/raw/ (e.g. 'brain_organoid.h5ad')

    Returns
    -------
    AnnData object
    """
    path = DATA_RAW / filename
    print(f"Loading: {path}")
    adata = sc.read_h5ad(path)
    print(f"  -> {adata.n_obs} cells x {adata.n_vars} genes")
    return adata


def load_10x_mtx(folder: str) -> ad.AnnData:
    """
    Load a 10x Genomics output folder (matrix.mtx, barcodes.tsv, features.tsv).

    Parameters
    ----------
    folder : str
        Subfolder name inside data/raw/ (e.g. 'sample1/')

    Returns
    -------
    AnnData object
    """
    path = DATA_RAW / folder
    print(f"Loading 10x data from: {path}")
    adata = sc.read_10x_mtx(path, var_names="gene_symbols", cache=True)
    adata.var_names_make_unique()
    print(f"  -> {adata.n_obs} cells x {adata.n_vars} genes")
    return adata


def save_processed(adata: ad.AnnData, filename: str):
    """
    Save an AnnData object to data/processed/.

    Parameters
    ----------
    adata    : AnnData object to save
    filename : e.g. 'organoid_preprocessed.h5ad'
    """
    path = DATA_PROC / filename
    adata.write_h5ad(path)
    print(f"Saved: {path}")


# ---------------------------------------------------------------------------
# Convenience
# ---------------------------------------------------------------------------

def summarize(adata: ad.AnnData, label: str = "Dataset"):
    """Print a quick summary of an AnnData object."""
    print(f"\n--- {label} ---")
    print(f"  Cells : {adata.n_obs}")
    print(f"  Genes : {adata.n_vars}")
    if adata.obs.columns.tolist():
        print(f"  Obs   : {adata.obs.columns.tolist()}")
    if adata.obsm.keys():
        print(f"  Embeddings: {list(adata.obsm.keys())}")
    print()
