"""Dimensionality reduction — PCA and UMAP. Mirrors DimReduce.R."""

from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from quantqc.core import QQC
from quantqc.plotting import dot_plot_style, um_plot_style, diverging_cmap


# ---------------------------------------------------------------------------
# PCA
# ---------------------------------------------------------------------------

def compute_pca(qqc: QQC, imputed: bool = True) -> QQC:
    """Compute PCA via eigenvalue decomposition of pairwise correlation matrix."""
    mat = qqc.matrices.protein_imputed if imputed else qqc.matrices.protein
    cols = qqc.matrices.protein_cols
    meta = qqc.meta_data

    # Pairwise correlation matrix (handles NaN)
    cor_mat = np.corrcoef(mat, rowvar=False)  # works on imputed data
    if np.any(np.isnan(cor_mat)):
        # Fallback: pairwise complete
        n = mat.shape[1]
        cor_mat = np.zeros((n, n))
        for i in range(n):
            for j in range(i, n):
                valid = ~np.isnan(mat[:, i]) & ~np.isnan(mat[:, j])
                if valid.sum() > 2:
                    cor_mat[i, j] = cor_mat[j, i] = np.corrcoef(mat[valid, i], mat[valid, j])[0, 1]
                else:
                    cor_mat[i, j] = cor_mat[j, i] = 0

    eigvals, eigvecs = np.linalg.eigh(cor_mat)
    # Sort descending
    order = np.argsort(-eigvals)
    eigvals = eigvals[order]
    eigvecs = eigvecs[:, order]

    pct_var = eigvals / eigvals.sum() * 100

    scx = pd.DataFrame(eigvecs, columns=[f"PC{i+1}" for i in range(eigvecs.shape[1])])
    scx["ID"] = cols
    scx = scx.merge(meta, on="ID", how="left")

    qqc.reductions["PCA"] = scx
    qqc.misc["pct_var"] = np.round(pct_var, 2)

    # Scree plot
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(range(1, len(pct_var) + 1), pct_var, "o-")
    ax.set_xlabel("PC")
    ax.set_ylabel("% of variance explained")
    plt.close(fig)

    return qqc


def plot_pca(qqc: QQC, by: str = "Condition") -> plt.Figure:
    """Plot PC1 vs PC2 colored by a metadata variable."""
    pca = qqc.reductions["PCA"]
    pv = qqc.misc["pct_var"]

    fig, ax = plt.subplots(figsize=(8, 6))

    if by == "Condition":
        for s in pca["sample"].unique():
            sub = pca[pca["sample"] == s]
            ax.scatter(sub["PC1"], sub["PC2"], label=s, s=20, alpha=0.7)
        ax.legend()
        ax.set_title("PCA by Condition")
    elif by == "Total protein":
        mask = pca["prot_total"].notna()
        sub = pca[mask]
        sc = ax.scatter(sub["PC1"], sub["PC2"], c=sub["prot_total"], cmap="RdBu_r",
                        norm=diverging_cmap(sub["prot_total"]), s=20)
        plt.colorbar(sc, ax=ax)
        ax.set_title("PCA by Total Cell Intensity")
    elif by == "Label":
        for lab in pca["label"].dropna().unique():
            sub = pca[pca["label"] == lab]
            ax.scatter(sub["PC1"], sub["PC2"], label=str(lab), s=20, alpha=0.7)
        ax.legend()
        ax.set_title("PCA by Label")
    elif by == "Run order":
        mask = pca["Order"].notna()
        sub = pca[mask]
        sc = ax.scatter(sub["PC1"], sub["PC2"], c=sub["Order"], cmap="RdBu_r",
                        norm=diverging_cmap(sub["Order"]), s=20)
        plt.colorbar(sc, ax=ax)
        ax.set_title("PCA by Run Order")

    ax.set_xlabel(f"PC1 ({pv[0]}%)")
    ax.set_ylabel(f"PC2 ({pv[1]}%)")
    dot_plot_style(ax)
    fig.tight_layout()
    return fig


def feature_pca(qqc: QQC, prot: str, imputed: bool = True) -> plt.Figure:
    """Overlay a single protein on PCA."""
    pca = qqc.reductions["PCA"]
    mat = qqc.matrices.protein_imputed if imputed else qqc.matrices.protein
    rows = qqc.matrices.protein_rows
    idx = np.where(rows == prot)[0]
    vals = mat[idx[0]] if len(idx) else np.full(len(pca), np.nan)

    fig, ax = plt.subplots(figsize=(8, 6))
    sc = ax.scatter(pca["PC1"], pca["PC2"], c=vals, cmap="RdBu_r",
                    norm=diverging_cmap(vals, midpoint=0), s=20)
    plt.colorbar(sc, ax=ax, label=prot)
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    dot_plot_style(ax)
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# UMAP
# ---------------------------------------------------------------------------

def compute_umap(qqc: QQC, n_pcs: int = 6, resolution: float = 0.5) -> QQC:
    """Compute UMAP embedding + graph-based clustering.

    Uses umap-learn for the embedding and Leiden/Louvain for clustering
    (falls back to spectral clustering if neither is installed).
    """
    try:
        import umap
    except ImportError:
        raise ImportError("umap-learn is required for UMAP. Install with: pip install umap-learn")

    from sklearn.decomposition import PCA

    mat = qqc.matrices.protein_imputed
    cols = qqc.matrices.protein_cols

    # PCA reduction (cells × genes)
    X = mat.T
    pca = PCA(n_components=n_pcs)
    X_pca = pca.fit_transform(X)

    # UMAP embedding
    reducer = umap.UMAP(n_components=2, random_state=42)
    X_umap = reducer.fit_transform(X_pca)

    # Graph-based clustering
    cluster_labels = _cluster_cells(X_pca, resolution=resolution)

    umap_df = pd.DataFrame(X_umap, columns=["umap_1", "umap_2"], index=cols)

    pca_df = qqc.reductions.get("PCA", pd.DataFrame())
    if not pca_df.empty and "ID" in pca_df.columns:
        pca_indexed = pca_df.set_index("ID")
        for src, dst in [("sample", "sample"), ("prot_total", "prot_total"),
                         ("label", "lab"), ("Order", "Order"), ("diameter", "diameter")]:
            if src in pca_indexed.columns:
                umap_df[dst] = pca_indexed.loc[cols, src].values

    umap_df["cluster"] = cluster_labels

    qqc.reductions["UMAP"] = umap_df
    return qqc


def _cluster_cells(X_pca: np.ndarray, resolution: float = 0.5,
                   n_neighbors: int = 15) -> np.ndarray:
    """Graph-based clustering: Leiden > Louvain > spectral fallback."""
    from sklearn.neighbors import kneighbors_graph

    # Build symmetric KNN connectivity graph
    knn = kneighbors_graph(X_pca, n_neighbors=n_neighbors,
                           mode="connectivity", include_self=False)
    knn = knn + knn.T
    knn.data[:] = 1  # binarize

    # Try Leiden (best match to Seurat/scanpy)
    try:
        import leidenalg
        import igraph as ig

        sources, targets = knn.nonzero()
        g = ig.Graph(n=X_pca.shape[0],
                     edges=list(zip(sources.tolist(), targets.tolist())),
                     directed=False)
        g.simplify()
        partition = leidenalg.find_partition(
            g, leidenalg.RBConfigurationVertexPartition,
            resolution_parameter=resolution,
        )
        return np.array(partition.membership, dtype=str)
    except ImportError:
        pass

    # Try python-louvain
    try:
        import community as community_louvain
        import networkx as nx

        G = nx.from_scipy_sparse_array(knn)
        partition = community_louvain.best_partition(G, resolution=resolution)
        return np.array([str(partition[i]) for i in range(X_pca.shape[0])])
    except ImportError:
        pass

    # Fallback: spectral clustering
    from sklearn.cluster import SpectralClustering
    n_clusters = max(2, min(10, X_pca.shape[0] // 50))
    sc = SpectralClustering(n_clusters=n_clusters, affinity="nearest_neighbors",
                            n_neighbors=n_neighbors, random_state=42)
    return np.array(sc.fit_predict(X_pca), dtype=str)


def plot_umap(qqc: QQC, by: str = "Cluster") -> plt.Figure:
    """Plot UMAP colored by a metadata variable."""
    umap = qqc.reductions["UMAP"]
    fig, ax = plt.subplots(figsize=(8, 6))

    if by == "Cluster":
        for cl in sorted(umap["cluster"].unique()):
            sub = umap[umap["cluster"] == cl]
            ax.scatter(sub["umap_1"], sub["umap_2"], label=str(cl), s=20, alpha=0.7)
        ax.legend()
    elif by == "Condition":
        for s in umap["sample"].unique():
            sub = umap[umap["sample"] == s]
            ax.scatter(sub["umap_1"], sub["umap_2"], label=s, s=20, alpha=0.7)
        ax.legend()
    elif by == "Total protein":
        mask = umap["prot_total"].notna()
        sub = umap[mask]
        sc_ = ax.scatter(sub["umap_1"], sub["umap_2"], c=sub["prot_total"], cmap="RdBu_r",
                         norm=diverging_cmap(sub["prot_total"]), s=20)
        plt.colorbar(sc_, ax=ax)
    elif by == "Label":
        for lab in umap["lab"].dropna().unique():
            sub = umap[umap["lab"] == lab]
            ax.scatter(sub["umap_1"], sub["umap_2"], label=str(lab), s=20, alpha=0.7)
        ax.legend()
    elif by == "Run order":
        mask = umap["Order"].notna()
        sub = umap[mask]
        sc_ = ax.scatter(sub["umap_1"], sub["umap_2"], c=sub["Order"], cmap="RdBu_r",
                         norm=diverging_cmap(sub["Order"]), s=20)
        plt.colorbar(sc_, ax=ax)
    elif by == "Diameter":
        mask = umap["diameter"].notna()
        sub = umap[mask]
        sc_ = ax.scatter(sub["umap_1"], sub["umap_2"], c=sub["diameter"], cmap="RdBu_r",
                         norm=diverging_cmap(sub["diameter"]), s=20)
        plt.colorbar(sc_, ax=ax)

    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.set_title(f"UMAP by {by}")
    um_plot_style(ax)
    fig.tight_layout()
    return fig


def feature_umap(qqc: QQC, prot: str, imputed: bool = True) -> plt.Figure:
    """Overlay a single protein on UMAP."""
    umap = qqc.reductions["UMAP"]
    mat = qqc.matrices.protein_imputed if imputed else qqc.matrices.protein
    rows = qqc.matrices.protein_rows
    idx = np.where(rows == prot)[0]
    vals = mat[idx[0]] if len(idx) else np.full(len(umap), np.nan)

    fig, ax = plt.subplots(figsize=(8, 6))
    sc_ = ax.scatter(umap["umap_1"], umap["umap_2"], c=vals, cmap="RdBu_r",
                     norm=diverging_cmap(vals, midpoint=0), s=20)
    plt.colorbar(sc_, ax=ax, label=prot)
    dot_plot_style(ax)
    fig.tight_layout()
    return fig


def feature_umap_abs(qqc: QQC, prot: str) -> plt.Figure:
    """Overlay absolute protein abundance on UMAP (log2)."""
    umap = qqc.reductions["UMAP"]
    mat = qqc.matrices.protein_abs
    rows = qqc.matrices.protein_rows
    idx = np.where(rows == prot)[0]
    vals = np.log2(mat[idx[0]]) if len(idx) else np.full(len(umap), np.nan)

    fig, ax = plt.subplots(figsize=(8, 6))
    sc_ = ax.scatter(umap["umap_1"], umap["umap_2"], c=vals, cmap="RdBu_r",
                     norm=diverging_cmap(vals), s=20)
    plt.colorbar(sc_, ax=ax, label=f"log2({prot})")
    dot_plot_style(ax)
    fig.tight_layout()
    return fig


def clust_box_plot(qqc: QQC, prot: str, imputed: bool = False) -> plt.Figure:
    """Boxplot of protein abundance by cluster."""
    umap = qqc.reductions["UMAP"]
    mat = qqc.matrices.protein_imputed if imputed else qqc.matrices.protein
    rows = qqc.matrices.protein_rows
    idx = np.where(rows == prot)[0]
    vals = mat[idx[0]] if len(idx) else np.full(len(umap), np.nan)

    df = umap[["cluster"]].copy()
    df["protein"] = vals
    clusters = sorted(df["cluster"].unique())

    fig, ax = plt.subplots(figsize=(8, 5))
    data_by_cluster = [df[df["cluster"] == c]["protein"].dropna().values for c in clusters]
    ax.boxplot(data_by_cluster, labels=[str(c) for c in clusters])
    ax.set_xlabel("Cluster")
    ax.set_ylabel(prot)
    dot_plot_style(ax)
    fig.tight_layout()
    return fig
