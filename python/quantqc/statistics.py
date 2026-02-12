"""Single-cell QC statistics and visualization — mirrors SC_statistics.R."""

from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

from quantqc.core import QQC
from quantqc.utils import (
    normalize_reference_vector, contains_missed_cleaved,
    extract_accession, create_peptide_vector, boxplot_summary,
)
from quantqc.plotting import dot_plot_style, diverging_cmap


# ---------------------------------------------------------------------------
# Peptide correlations
# ---------------------------------------------------------------------------

def shared_peptide_cor(qqc: QQC, n_obs_rec: int = 20, res: str = "sc") -> QQC:
    """Compute peptide correlations within proteins (+ null distribution)."""
    peptide_data = qqc.matrices.peptide
    protein_dat = qqc.matrices.protein
    ppm = qqc.matrices.peptide_protein_map
    protein_rows = qqc.matrices.protein_rows

    if res == "clust":
        clusters = qqc.reductions["UMAP"]["cluster"]
        unique_cl = sorted(clusters.unique())
        clust_mat = np.zeros((peptide_data.shape[0], len(unique_cl)))
        for ci, cl in enumerate(unique_cl):
            mask = (clusters == cl).values
            clust_mat[:, ci] = np.nanmean(peptide_data[:, mask], axis=1)
        peptide_data = clust_mat

    proteins = ppm["Protein"].values
    unique_prots = np.unique(proteins)

    results = []
    for prot in unique_prots:
        mask = proteins == prot
        sub = peptide_data[mask]
        if sub.shape[0] < 2:
            continue
        # Pairwise counts and correlations
        n = sub.shape[0]
        cors = []
        for a in range(n):
            for b in range(a + 1, n):
                valid = ~np.isnan(sub[a]) & ~np.isnan(sub[b])
                if valid.sum() < n_obs_rec:
                    continue
                r = np.corrcoef(sub[a, valid], sub[b, valid])[0, 1]
                cors.append(r)
        if cors:
            results.append({"Protein": prot, "Cor": float(np.nanmedian(cors)), "Obs": n})

    # Null distribution: shuffle protein labels
    shuffled = np.random.permutation(proteins)
    null_cors = []
    for prot in unique_prots:
        mask = shuffled == prot
        sub = peptide_data[mask]
        if sub.shape[0] < 2:
            continue
        n = sub.shape[0]
        for a in range(n):
            for b in range(a + 1, n):
                valid = ~np.isnan(sub[a]) & ~np.isnan(sub[b])
                if valid.sum() < 4:
                    continue
                r = np.corrcoef(sub[a, valid], sub[b, valid])[0, 1]
                null_cors.append(r)

    pep_cor_df = pd.DataFrame(results)
    if pep_cor_df.empty:
        qqc.pep_cor = [pep_cor_df, 0.0]
        return qqc

    # Match to protein matrix
    pep_cor_df = pep_cor_df[pep_cor_df["Protein"].isin(protein_rows)]

    # FC binning
    prot_abs_mean = []
    for prot in pep_cor_df["Protein"]:
        idx = np.where(protein_rows == prot)[0]
        if len(idx):
            prot_abs_mean.append(float(np.nanmean(np.abs(protein_dat[idx[0]]))))
        else:
            prot_abs_mean.append(0.0)
    pep_cor_df["FC"] = prot_abs_mean

    bins = [0, 0.4, 0.8, 1.2, 2.0, np.inf]
    labels = ["0.4", "0.8", "1.2", "2", "3"]
    pep_cor_df["FC"] = pd.cut(pep_cor_df["FC"], bins=bins, labels=labels).astype(str)

    null_median = float(np.nanmedian(null_cors)) if null_cors else 0.0
    qqc.pep_cor = [pep_cor_df, null_median]
    return qqc


def plot_pep_cor(qqc: QQC, type_: str = "box") -> plt.Figure:
    """Plot peptide correlations."""
    pep_cor = qqc.pep_cor[0]
    null_dist = qqc.pep_cor[1]

    fig, ax = plt.subplots(figsize=(8, 5))

    if type_ == "box":
        fc_vals = sorted(pep_cor["FC"].unique())
        data = [pep_cor[pep_cor["FC"] == fc]["Cor"].dropna().values for fc in fc_vals]
        bp = ax.boxplot(data, labels=fc_vals, patch_artist=True)
        for patch in bp["boxes"]:
            patch.set_facecolor("lightgray")
        ax.axhline(null_dist, color="red", lw=1)
        ax.set_xlabel("Mean abs(protein fold change)")
        ax.set_ylabel("Cor.; peptides mapping to protein")
        # Annotate counts
        for i, fc in enumerate(fc_vals):
            n = len(pep_cor[pep_cor["FC"] == fc])
            ax.text(i + 1, ax.get_ylim()[1] * 0.95, str(n), ha="center", color="blue")
    elif type_ == "hist":
        ax.hist(pep_cor["Cor"].dropna(), bins=30, color="gray", edgecolor="black")
        ax.axvline(null_dist, color="red", lw=1)
        ax.set_xlabel("Correlation; peptides mapping to protein")
        ax.set_ylabel("# Proteins")

    dot_plot_style(ax)
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Peptide / protein counts
# ---------------------------------------------------------------------------

def plot_prot_and_pep(qqc: QQC) -> plt.Figure:
    """Histogram of peptides and proteins per cell."""
    pep_counts = np.sum(~np.isnan(qqc.matrices.peptide), axis=0)
    prot_counts = np.sum(~np.isnan(qqc.matrices.protein), axis=0)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    if qqc.ms_type in ("DIA", "DIA_C"):
        # Also show unfiltered counts via mask
        pep_nf = np.sum(qqc.matrices.peptide_mask, axis=0)
        ax1.hist(pep_nf, bins=30, alpha=0.5, label="No Ch Q filter")
        prot_nf = np.sum(qqc.matrices.protein_mask, axis=0)
        ax2.hist(prot_nf, bins=30, alpha=0.5, label="No Ch Q filter")

    ax1.hist(pep_counts, bins=15, alpha=0.5, label=f"{qqc.misc.get('ChQ', '')} Ch Q")
    ax1.set_xlabel("# precursors")
    ax1.set_ylabel("# of single cells")
    ax1.set_title("# precursors per sample")
    ax1.legend()
    dot_plot_style(ax1)

    ax2.hist(prot_counts, bins=15, alpha=0.5, label=f"{qqc.misc.get('ChQ', '')} Ch Q")
    ax2.set_xlabel("# proteins")
    ax2.set_ylabel("# of single cells")
    ax2.set_title("# proteins per sample")
    ax2.legend()
    dot_plot_style(ax2)

    fig.tight_layout()
    return fig


def plot_data_complete(qqc: QQC) -> plt.Figure:
    """Histogram of data completeness per cell and per protein."""
    data = qqc.matrices.protein
    cell_complete = np.sum(~np.isnan(data), axis=0) / data.shape[0]
    prot_complete = np.sum(~np.isnan(data), axis=1) / data.shape[1]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    ax1.hist(prot_complete, bins=20, alpha=0.5)
    ax1.set_xlabel("fraction values present")
    ax1.set_ylabel("# of proteins")
    ax1.set_title(f"Protein completeness, {data.shape[0]} proteins")
    dot_plot_style(ax1)

    ax2.hist(cell_complete, bins=20, alpha=0.5)
    ax2.set_xlabel("fraction values present")
    ax2.set_ylabel("# of single cells")
    ax2.set_title("Cell completeness")
    dot_plot_style(ax2)

    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# SC to carrier ratio
# ---------------------------------------------------------------------------

def plot_sc_to_carrier_ratio(qqc: QQC) -> plt.Figure | str:
    """Plot single cell to carrier ratio."""
    if qqc.ms_type == "DDA":
        # Simplified version
        meta = qqc.meta_data
        fig, ax = plt.subplots(figsize=(7, 5))
        samples = meta["sample"].unique()
        data_by = []
        labels_ = []
        for s in samples:
            sub = meta[meta["sample"] == s]
            if "prot_total" in sub.columns:
                data_by.append(sub["prot_total"].dropna().values)
                labels_.append(s)
        ax.boxplot(data_by, labels=labels_)
        ax.set_ylabel("Single cells / Carrier")
        dot_plot_style(ax)
        fig.tight_layout()
        return fig

    if qqc.ms_type == "DIA_C":
        mat = qqc.matrices.peptide
        vals = np.nanmedian(mat, axis=0)
        meta = qqc.meta_data
        meta_sub = meta[meta["ID"].isin(qqc.matrices.peptide_cols)].copy()
        meta_sub["ratio"] = 1.0 / vals[:len(meta_sub)]

        fig, ax = plt.subplots(figsize=(7, 5))
        for s in meta_sub["sample"].unique():
            sub = meta_sub[meta_sub["sample"] == s]
            ax.boxplot(sub["ratio"].dropna().values, positions=[list(meta_sub["sample"].unique()).index(s)])
        ax.set_ylabel("# single cells in carrier")
        ax.set_title(f"Carrier is median {round(float(np.nanmedian(1/np.nanmedian(mat, axis=0))), 2)} single cells")
        dot_plot_style(ax)
        fig.tight_layout()
        return fig

    return "No carrier used"


# ---------------------------------------------------------------------------
# Digest efficiency
# ---------------------------------------------------------------------------

def plot_digest_eff(qqc: QQC) -> plt.Figure | str:
    """Plot protease digestion efficiency (regular vs missed cleaved)."""
    if qqc.ms_type not in ("DDA", "DIA_C"):
        return "No carrier used"

    mat = qqc.matrices.peptide.copy()
    peps = qqc.matrices.peptide_protein_map["seqcharge"]
    trimmed = peps.str[:-3] if qqc.ms_type == "DDA" else peps.str[:-2]

    mc = contains_missed_cleaved(trimmed)

    # Column-normalise
    for i in range(mat.shape[1]):
        med = np.nanmedian(mat[:, i])
        if med and not np.isnan(med):
            mat[:, i] /= med

    reg_vals = mat[mc.values == 0].ravel()
    mc_vals = mat[mc.values != 0].ravel()
    reg_vals = reg_vals[~np.isnan(reg_vals)]
    mc_vals = mc_vals[~np.isnan(mc_vals)]

    dig_eff = round(float(np.median(reg_vals) / np.median(mc_vals)), 2) if np.median(mc_vals) != 0 else 0

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.boxplot([np.log2(reg_vals), np.log2(mc_vals)], labels=["Regular", "Missed Cleaved"])
    ax.set_ylabel("log2(Intensity vs Carrier)")
    ax.set_title(f"{dig_eff}X digest efficiency of carrier")
    dot_plot_style(ax)
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Cell size vs intensity
# ---------------------------------------------------------------------------

def plot_cell_size_vs_intensity(qqc: QQC, type_: str = "sample") -> plt.Figure:
    """Plot cell diameter vs total protein intensity."""
    meta = qqc.meta_data.copy()
    good_cells = qqc.matrices.peptide_cols
    meta = meta[meta["ID"].isin(good_cells)]

    vol = np.log2((meta["diameter"].astype(float) / 2) ** 3)
    intens = np.log2(10 ** meta["prot_total"].astype(float))
    corr = round(float(np.corrcoef(
        meta["diameter"].astype(float).dropna(),
        meta["prot_total"].astype(float).dropna()
    )[0, 1]), 2)

    fig, ax = plt.subplots(figsize=(8, 6))

    if type_ == "sample":
        for s in meta["sample"].unique():
            sub = meta[meta["sample"] == s]
            v = np.log2((sub["diameter"].astype(float) / 2) ** 3)
            it = np.log2(10 ** sub["prot_total"].astype(float))
            ax.scatter(v, it, label=s, s=20, alpha=0.7)
        ax.legend()
    elif type_ == "Run order":
        mask = meta["Order"].notna()
        sc_ = ax.scatter(vol[mask], intens[mask], c=meta.loc[mask, "Order"], cmap="RdBu_r",
                         norm=diverging_cmap(meta.loc[mask, "Order"]), s=20)
        plt.colorbar(sc_, ax=ax)

    ax.set_xlabel("log2(Vol.) cubic uM")
    ax.set_ylabel("log2(Sum cell intensity)")
    ax.set_title(f"Correlation = {corr}")
    dot_plot_style(ax)
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Protein cluster consistency
# ---------------------------------------------------------------------------

def protein_clust_consistency(qqc: QQC, prot: str, type_: str = "line",
                              fasta_path: str | None = None) -> plt.Figure | str:
    """Plot peptide consistency across clusters for a protein."""
    from Bio import SeqIO

    if fasta_path is None:
        species = qqc.misc.get("Species", "Human")
        from quantqc.cellenone import _get_extdata
        fasta_path = _get_extdata(f"{species}.fasta")

    fasta_dict = {}
    for record in SeqIO.parse(fasta_path, "fasta"):
        acc = extract_accession(record.description)
        fasta_dict[acc] = str(record.seq)

    if prot not in fasta_dict:
        return "Protein not found in fasta"

    prot_seq = fasta_dict[prot]
    clusters = qqc.reductions["UMAP"]
    ppm = qqc.matrices.peptide_protein_map
    prot_peps = ppm[ppm["Protein"] == prot]["seqcharge"].values

    if len(prot_peps) < 2:
        return "Only 1 peptide"
    if len(prot_peps) > 4:
        prot_peps = prot_peps[:4]

    pep_mat = normalize_reference_vector(qqc.matrices.peptide, log=True)
    pep_rows = qqc.matrices.peptide_rows

    fig, axes = plt.subplots(1, len(prot_peps), figsize=(5 * len(prot_peps), 5))
    if len(prot_peps) == 1:
        axes = [axes]

    colors = ["red", "blue", "green", "purple"]
    cluster_vals = clusters["cluster"]
    unique_clusters = sorted(cluster_vals.unique())

    for pi, pep in enumerate(prot_peps):
        ax = axes[pi]
        idx = np.where(pep_rows == pep)[0]
        if not len(idx):
            continue
        vals = pep_mat[idx[0]]
        # Protein-level median per cluster
        prot_idx = np.where(pep_rows == prot_peps[0])[0]
        for cl in unique_clusters:
            cl_mask = (cluster_vals == cl).values
            pep_med = float(np.nanmedian(vals[cl_mask]))
            ax.scatter(cl, pep_med, color=colors[pi], s=40)
        ax.set_xlabel("Cluster")
        ax.set_ylabel("Log2(FoldChange)")
        ax.set_title(pep[:20])
        dot_plot_style(ax)

    fig.suptitle(prot)
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# MS1 vs MS2 correlation
# ---------------------------------------------------------------------------

def plot_ms1_vs_ms2(qqc: QQC) -> plt.Figure:
    """Plot MS1 vs MS2 peptide correlation (DIA only)."""
    raw = qqc.raw_data
    raw = raw[raw["Channel.Q.Value"] < qqc.misc.get("ChQ", 0.1)]
    pep_mat = qqc.matrices.peptide
    pep_rows = qqc.matrices.peptide_rows
    pep_cols = qqc.matrices.peptide_cols

    raw_filt = raw[raw["ID"].isin(pep_cols)]

    ms2 = raw_filt.pivot_table(index="seqcharge", columns="ID", values="Precursor.Quantity", aggfunc="first")
    sect_pep = np.intersect1d(ms2.index.values, pep_rows)
    sect_col = np.intersect1d(ms2.columns.values, pep_cols)

    row_idx = np.array([np.where(pep_rows == p)[0][0] for p in sect_pep])
    col_idx = np.array([np.where(pep_cols == c)[0][0] for c in sect_col])

    ms2_mat = ms2.loc[sect_pep, sect_col].values.astype(np.float64)
    pep_sub = pep_mat[np.ix_(row_idx, col_idx)]

    ms2_norm = normalize_reference_vector(ms2_mat, log=True)
    pep_norm = normalize_reference_vector(pep_sub, log=True)

    cors, fc_vals = [], []
    for i in range(ms2_norm.shape[0]):
        valid = ~np.isnan(ms2_norm[i]) & ~np.isnan(pep_norm[i])
        if valid.sum() > 10:
            r = np.corrcoef(ms2_norm[i, valid], pep_norm[i, valid])[0, 1]
            cors.append(r)
            fc_vals.append(float(np.nanmean(np.abs(pep_norm[i]))))

    df = pd.DataFrame({"cors": cors, "FC": fc_vals})
    bins = [0, 0.4, 0.8, 1.2, 1.6, np.inf]
    labels = ["0.4", "0.8", "1.2", "1.6", "2"]
    df["FC"] = pd.cut(df["FC"], bins=bins, labels=labels).astype(str)

    fig, ax = plt.subplots(figsize=(8, 5))
    fc_groups = sorted(df["FC"].unique())
    data = [df[df["FC"] == fc]["cors"].dropna().values for fc in fc_groups]
    ax.boxplot(data, labels=fc_groups, patch_artist=True)
    ax.set_xlabel("Average abs(log2 Fold Change)")
    ax.set_ylabel("Correlations")
    ax.set_title("Ms1, Ms2 Peptide Correlation")
    dot_plot_style(ax)
    fig.tight_layout()
    return fig
