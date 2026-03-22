"""
visualization.py
----------------
Reusable plotting functions for scRNA-seq analysis.

All functions save figures to results/figures/ automatically
and also display them inline in notebooks.
"""

import matplotlib.pyplot as plt
import scanpy as sc
import anndata as ad
from pathlib import Path

# Project results directory (resolved relative to this file)
FIGURES_DIR = Path(__file__).resolve().parent.parent / "results" / "figures"


def _save(fig: plt.Figure, filename: str):
    """Save a figure to the results/figures/ directory."""
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    path = FIGURES_DIR / filename
    fig.savefig(path, dpi=150, bbox_inches="tight")
    print(f"Saved: {path}")


def plot_qc_metrics(adata: ad.AnnData, save: bool = True):
    """
    Violin plots of key QC metrics:
      - Number of genes per cell
      - Total counts per cell
      - % mitochondrial counts per cell

    These help you choose appropriate filtering thresholds.
    """
    fig, axes = plt.subplots(1, 3, figsize=(12, 4))
    fig.suptitle("QC Metrics", fontsize=14)

    metrics = [
        ("n_genes_by_counts", "Genes per cell"),
        ("total_counts",      "Total counts per cell"),
        ("pct_counts_mt",     "% Mitochondrial counts"),
    ]

    for ax, (metric, label) in zip(axes, metrics):
        ax.violinplot(adata.obs[metric], showmedians=True)
        ax.set_title(label)
        ax.set_xticks([])
        ax.set_ylabel(label)

    plt.tight_layout()

    if save:
        _save(fig, "qc_metrics.png")

    plt.show()


def plot_qc_scatter(adata: ad.AnnData, save: bool = True):
    """
    Scatter plots to spot outlier cells:
      - Total counts vs genes detected
      - Total counts vs % mitochondrial

    Cells far from the main cloud may be doublets or low-quality.
    """
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    fig.suptitle("QC Scatter Plots", fontsize=14)

    axes[0].scatter(
        adata.obs["total_counts"],
        adata.obs["n_genes_by_counts"],
        s=1, alpha=0.3
    )
    axes[0].set_xlabel("Total counts")
    axes[0].set_ylabel("Genes detected")
    axes[0].set_title("Counts vs Genes")

    axes[1].scatter(
        adata.obs["total_counts"],
        adata.obs["pct_counts_mt"],
        s=1, alpha=0.3, color="salmon"
    )
    axes[1].set_xlabel("Total counts")
    axes[1].set_ylabel("% Mitochondrial")
    axes[1].set_title("Counts vs Mito %")

    plt.tight_layout()

    if save:
        _save(fig, "qc_scatter.png")

    plt.show()


def plot_hvg(adata: ad.AnnData, save: bool = True):
    """
    Plot the highly variable genes.
    Each dot is a gene; red dots are selected as highly variable.
    """
    sc.pl.highly_variable_genes(adata, show=False)
    fig = plt.gcf()
    fig.suptitle("Highly Variable Genes", y=1.01)

    if save:
        _save(fig, "highly_variable_genes.png")

    plt.show()


def plot_pca_variance(adata: ad.AnnData, n_pcs: int = 30, save: bool = True):
    """
    Scree plot: variance explained by each principal component.
    Helps decide how many PCs to use for neighbor graph and UMAP.
    Look for the 'elbow' where explained variance levels off.
    """
    variance_ratio = adata.uns["pca"]["variance_ratio"][:n_pcs]
    cumulative     = variance_ratio.cumsum()

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    fig.suptitle("PCA Variance", fontsize=14)

    axes[0].bar(range(1, n_pcs + 1), variance_ratio)
    axes[0].set_xlabel("Principal Component")
    axes[0].set_ylabel("Variance Ratio")
    axes[0].set_title("Variance per PC")

    axes[1].plot(range(1, n_pcs + 1), cumulative, marker="o", markersize=3)
    axes[1].axhline(0.9, color="red", linestyle="--", label="90% variance")
    axes[1].set_xlabel("Principal Component")
    axes[1].set_ylabel("Cumulative Variance")
    axes[1].set_title("Cumulative Variance")
    axes[1].legend()

    plt.tight_layout()

    if save:
        _save(fig, "pca_variance.png")

    plt.show()


def plot_umap(
    adata: ad.AnnData,
    color: list = None,
    title: str = "UMAP",
    save_name: str = "umap.png",
    save: bool = True,
):
    """
    Plot UMAP embedding colored by one or more variables.

    Parameters
    ----------
    color     : list of obs column names or gene names to color by
                e.g. ['leiden', 'n_genes_by_counts', 'NEUROD2']
    title     : plot title
    save_name : filename to save to
    """
    if color is None:
        color = ["leiden"]

    sc.pl.umap(adata, color=color, title=title, show=False)
    fig = plt.gcf()

    if save:
        _save(fig, save_name)

    plt.show()
