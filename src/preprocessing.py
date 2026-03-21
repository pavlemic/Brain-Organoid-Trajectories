"""
preprocessing.py
----------------
All preprocessing steps for scRNA-seq data:
  1. QC metrics calculation
  2. Cell and gene filtering
  3. Normalization
  4. Log-transformation
  5. Highly variable gene (HVG) selection
  6. PCA

Each function takes an AnnData object, modifies it in-place, and returns it
so you can chain calls together.
"""

import scanpy as sc
import anndata as ad


def calculate_qc_metrics(adata: ad.AnnData) -> ad.AnnData:
    """
    Add quality control metrics to adata.obs:
      - n_genes_by_counts  : number of genes detected per cell
      - total_counts       : total UMI counts per cell
      - pct_counts_mt      : % of counts from mitochondrial genes

    Mitochondrial genes are flagged by the 'MT-' prefix (human).
    Use 'mt-' for mouse data.
    """
    # Flag mitochondrial genes
    adata.var["mt"] = adata.var_names.str.startswith("MT-")

    # scanpy computes n_genes_by_counts, total_counts, pct_counts_mt etc.
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt"],
        percent_top=None,
        log1p=False,
        inplace=True,
    )

    print("QC metrics calculated.")
    print(f"  Median genes/cell : {adata.obs['n_genes_by_counts'].median():.0f}")
    print(f"  Median counts/cell: {adata.obs['total_counts'].median():.0f}")
    print(f"  Median %mito      : {adata.obs['pct_counts_mt'].median():.1f}%")

    return adata


def filter_cells_and_genes(
    adata: ad.AnnData,
    min_genes: int = 200,       # drop cells with fewer than this many genes
    max_genes: int = 6000,      # drop likely doublets (cells with too many genes)
    max_pct_mt: float = 20.0,   # drop low-quality/dying cells with high mito %
    min_cells: int = 3,         # drop genes detected in fewer than this many cells
) -> ad.AnnData:
    """
    Filter out low-quality cells and lowly expressed genes.

    Parameters
    ----------
    min_genes   : minimum genes detected per cell
    max_genes   : maximum genes detected per cell (doublet threshold)
    max_pct_mt  : maximum mitochondrial gene percentage per cell
    min_cells   : minimum cells a gene must be detected in
    """
    n_before = adata.n_obs

    # Filter cells
    sc.pp.filter_cells(adata, min_genes=min_genes)
    adata = adata[adata.obs["n_genes_by_counts"] < max_genes].copy()
    adata = adata[adata.obs["pct_counts_mt"] < max_pct_mt].copy()

    # Filter genes
    sc.pp.filter_genes(adata, min_cells=min_cells)

    n_after = adata.n_obs
    print(f"Filtering: {n_before} -> {n_after} cells retained "
          f"({n_before - n_after} removed).")
    print(f"  Genes remaining: {adata.n_vars}")

    return adata


def normalize_and_log(adata: ad.AnnData, target_sum: float = 1e4) -> ad.AnnData:
    """
    Normalize each cell to the same total count, then log-transform.

    Steps:
      1. Normalize total counts per cell to `target_sum` (e.g. 10,000)
         This corrects for differences in sequencing depth between cells.
      2. Apply log1p: log(x + 1)
         This stabilizes variance and makes the data more normally distributed.
    """
    # Store raw counts before normalization (useful for differential expression later)
    adata.raw = adata

    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)

    print(f"Normalized to {target_sum:.0f} counts/cell and log1p transformed.")
    return adata


def select_highly_variable_genes(
    adata: ad.AnnData,
    n_top_genes: int = 2000,
) -> ad.AnnData:
    """
    Identify the most variable genes across cells.

    Why: Most genes don't vary much between cell types. Focusing on the top
    highly variable genes (HVGs) reduces noise and speeds up downstream steps.

    Parameters
    ----------
    n_top_genes : number of top variable genes to select
    """
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=n_top_genes,
        flavor="seurat",   # standard method used in most scRNA-seq pipelines
        subset=False,      # keep all genes but flag HVGs in adata.var
    )

    n_hvg = adata.var["highly_variable"].sum()
    print(f"Highly variable genes selected: {n_hvg} / {adata.n_vars}")
    return adata


def run_pca(adata: ad.AnnData, n_comps: int = 50) -> ad.AnnData:
    """
    Run PCA on the highly variable genes.

    PCA reduces dimensionality from thousands of genes to a smaller number
    of principal components (PCs) that capture the main axes of variation.
    This is the input to neighbor graph construction and UMAP.

    Parameters
    ----------
    n_comps : number of principal components to compute
    """
    # use_highly_variable=True restricts PCA to the HVGs we selected above
    sc.pp.scale(adata, max_value=10)   # scale genes to unit variance
    sc.tl.pca(adata, n_comps=n_comps, use_highly_variable=True)

    print(f"PCA complete: {n_comps} components computed.")
    return adata


def run_full_preprocessing(
    adata: ad.AnnData,
    min_genes: int = 200,
    max_genes: int = 6000,
    max_pct_mt: float = 20.0,
    min_cells: int = 3,
    n_top_genes: int = 2000,
    n_pcs: int = 50,
) -> ad.AnnData:
    """
    Convenience wrapper: run the full preprocessing pipeline in one call.

    Order:
      QC metrics -> filter -> normalize -> log -> HVG -> PCA
    """
    adata = calculate_qc_metrics(adata)
    adata = filter_cells_and_genes(adata, min_genes, max_genes, max_pct_mt, min_cells)
    adata = normalize_and_log(adata)
    adata = select_highly_variable_genes(adata, n_top_genes)
    adata = run_pca(adata, n_pcs)
    return adata
