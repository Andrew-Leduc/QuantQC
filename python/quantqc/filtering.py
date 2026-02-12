"""Cell QC and peptide filtering — mirrors FilteringData.R."""

from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from quantqc.core import QQC, mat_to_df
from quantqc.utils import normalize_reference_vector, fast_cv
from quantqc.plotting import dot_plot_style, MY_COL3


# ---------------------------------------------------------------------------
# CV computation
# ---------------------------------------------------------------------------

def _compute_cvs(qqc: QQC) -> pd.DataFrame:
    """Compute median protein CV per cell (mirrors R ``CVs``)."""
    mat_norm = normalize_reference_vector(qqc.matrices.peptide)
    ppm = qqc.matrices.peptide_protein_map

    # Build long-form DataFrame
    df = pd.DataFrame(mat_norm, columns=qqc.matrices.peptide_cols)
    df["Protein"] = ppm["Protein"].values
    df["pep"] = ppm["seqcharge"].values
    long = df.melt(id_vars=["Protein", "pep"], var_name="variable", value_name="value")

    cv_mat = fast_cv(long, protein_col="Protein")
    cv_mat = cv_mat.drop(columns=["value"], errors="ignore")

    # Count proteins with valid CV per cell
    cv_count = (
        cv_mat.dropna(subset=["cvq"])
        .groupby("variable")
        .size()
        .reset_index(name="counter")
    )
    cv_count = cv_count[cv_count["counter"] > 20]

    # Median CV per cell
    cv_med = (
        cv_mat.groupby("variable")["cvq"]
        .median()
        .reset_index()
    )
    cv_med = cv_med[cv_med["variable"].isin(cv_count["variable"])]

    # Label cells vs negatives
    pos_ids = set(qqc.meta_data[qqc.meta_data["sample"] != "neg"]["ID"])
    cv_med["value"] = cv_med["variable"].apply(lambda x: "cell" if x in pos_ids else "neg")

    return cv_med


def _count_peptides_per_cell(peptide_mat, peptide_cols, meta, good_cells=None):
    """Count peptides and total intensity per cell, split by cell vs negative control."""
    neg_ids = set(meta[meta["sample"] == "neg"]["ID"])

    # Negative controls
    neg_mask = np.isin(peptide_cols, list(neg_ids))
    neg_count = np.sum(~np.isnan(peptide_mat[:, neg_mask]), axis=0) if neg_mask.any() else np.array([0])
    neg_intense = np.nansum(peptide_mat[:, neg_mask], axis=0) if neg_mask.any() else np.array([0])
    neg_df = pd.DataFrame({"Number_precursors": neg_count, "intense": neg_intense, "type": "negative ctrl"})

    # Cells
    if good_cells is not None:
        cell_mask = np.isin(peptide_cols, good_cells)
    else:
        cell_ids = set(meta[meta["sample"] != "neg"]["ID"])
        cell_mask = np.isin(peptide_cols, list(cell_ids))
    cell_count = np.sum(~np.isnan(peptide_mat[:, cell_mask]), axis=0)
    cell_intense = np.nansum(peptide_mat[:, cell_mask], axis=0)
    pos_df = pd.DataFrame({"Number_precursors": cell_count, "intense": cell_intense, "type": "single cells"})

    result = pd.concat([pos_df, neg_df], ignore_index=True)
    result["variable"] = (
        list(peptide_cols[cell_mask]) + list(peptide_cols[neg_mask])
        if neg_mask.any()
        else list(peptide_cols[cell_mask]) + ["neg_0"]
    )
    return result


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def evaluate_negative_controls(qqc: QQC) -> QQC:
    """Evaluate negative controls — dispatches to DDA or DIA variant."""
    if qqc.meta_data.empty:
        raise ValueError("Must map sample identities to data first")

    if qqc.ms_type == "DDA":
        return _evaluate_neg_dda(qqc)
    return _evaluate_neg_dia(qqc)


def _evaluate_neg_dda(qqc: QQC) -> QQC:
    cv_mat = _compute_cvs(qqc)
    good_cells = cv_mat[(cv_mat["cvq"] < 0.9) & (cv_mat["value"] != "neg")]
    if len(good_cells) < 3:
        raise ValueError("Less than three good cells, try increasing CV filter")

    counts = _count_peptides_per_cell(
        qqc.matrices.peptide, qqc.matrices.peptide_cols,
        qqc.meta_data, good_cells=good_cells["variable"].values,
    )
    neg_meta = cv_mat.merge(counts, on="variable", how="left")
    qqc.neg_ctrl_info = neg_meta
    return qqc


def _evaluate_neg_dia(qqc: QQC) -> QQC:
    counts = _count_peptides_per_cell(
        qqc.matrices.peptide, qqc.matrices.peptide_cols, qqc.meta_data,
    )
    qqc.neg_ctrl_info = counts
    return qqc


def plot_neg_ctrl(qqc: QQC, cv_thresh: float = 0.4) -> plt.Figure:
    """Visualize negative control metrics."""
    plot_data = qqc.neg_ctrl_info

    if qqc.ms_type == "DDA":
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

        # Intensity histogram
        for t, color in zip(["single cells", "negative ctrl"], ["tab:blue", "tab:orange"]):
            sub = plot_data[plot_data["type"] == t]
            ax1.hist(np.log10(sub["intense"].clip(lower=1)), bins=30, alpha=0.5, label=t, color=color)
        ax1.set_xlabel("sum(log10(Intensity))")
        ax1.set_ylabel("# of samples")
        dot_plot_style(ax1)
        ax1.legend()

        # CV density
        for val, color in zip(["cell", "neg"], MY_COL3):
            sub = plot_data[plot_data["value"] == val]
            ax2.hist(sub["cvq"].dropna(), bins=40, alpha=0.5, density=True, label=val, color=color)
        ax2.axvline(cv_thresh, ls="--", lw=2, color="gray")
        ax2.set_xlabel("CV, peptides from same protein")
        ax2.set_ylabel("Density")
        n_good = ((plot_data["value"] == "cell") & (plot_data["cvq"] < cv_thresh)).sum()
        ax2.set_title(f"{n_good} cells pass CV < {cv_thresh}")
        dot_plot_style(ax2)

        fig.tight_layout()
        return fig

    # DIA
    fig, ax = plt.subplots(figsize=(7, 5))
    for t, color in zip(["single cells", "negative ctrl"], ["tab:blue", "tab:orange"]):
        sub = plot_data[plot_data["type"] == t]
        ax.hist(np.log10(sub["intense"].clip(lower=1)), bins=30, alpha=0.5, label=t, color=color)
    ax.set_xlabel("log10(Intensity)")
    ax.set_ylabel("# of samples")
    ax.set_title("Negative ctrl vs Single cells")
    ax.legend()
    dot_plot_style(ax)
    fig.tight_layout()
    return fig


def filter_bad_cells(qqc: QQC, cv_thresh: float = None, min_intens: float = None) -> QQC:
    """Remove low-quality cells based on CV and/or intensity thresholds."""
    neg_filter = qqc.neg_ctrl_info.copy()

    if qqc.ms_type in ("DIA", "DIA_C"):
        neg_filter = neg_filter[neg_filter["type"] != "negative ctrl"]
        if min_intens is not None:
            neg_filter = neg_filter[np.log10(neg_filter["intense"].clip(lower=1)) > min_intens]
        keep = neg_filter["variable"].values
        col_mask = np.isin(qqc.matrices.peptide_cols, keep)
        qqc.matrices.peptide = qqc.matrices.peptide[:, col_mask]
        qqc.matrices.peptide_cols = qqc.matrices.peptide_cols[col_mask]
        if hasattr(qqc.matrices, "peptide_mask") and qqc.matrices.peptide_mask.size:
            qqc.matrices.peptide_mask = qqc.matrices.peptide_mask[:, col_mask]
    else:  # DDA
        neg_filter = neg_filter[neg_filter["value"] != "neg"]
        if min_intens is not None:
            neg_filter = neg_filter[np.log10(neg_filter["intense"].clip(lower=1)) > min_intens]
        if cv_thresh is not None:
            neg_filter = neg_filter[neg_filter["cvq"] < cv_thresh]
        keep = neg_filter["variable"].values
        col_mask = np.isin(qqc.matrices.peptide_cols, keep)
        qqc.matrices.peptide = qqc.matrices.peptide[:, col_mask]
        qqc.matrices.peptide_cols = qqc.matrices.peptide_cols[col_mask]

    return qqc


def trim_extra_peptides(qqc: QQC) -> QQC:
    """Keep at most 5 best peptides per protein (by median intensity + coverage)."""
    peptide_mat = qqc.matrices.peptide
    rows = qqc.matrices.peptide_rows
    ppm = qqc.matrices.peptide_protein_map

    keep_idx = []
    for prot in ppm["Protein"].unique():
        prot_mask = (ppm["Protein"] == prot).values
        prot_indices = np.where(prot_mask)[0]
        if len(prot_indices) <= 5:
            keep_idx.extend(prot_indices.tolist())
            continue
        sub = peptide_mat[prot_indices]
        medians = np.nanmedian(sub, axis=1)
        coverage = np.sum(~np.isnan(sub), axis=1)
        top_med = np.argsort(-medians)[:5]
        top_cov = np.argsort(-coverage)[:5]
        sect = set(top_med) & set(top_cov)
        if len(sect) < 4:
            chosen = top_med[:5]
        else:
            chosen = list(sect)[:5]
        keep_idx.extend(prot_indices[chosen].tolist())

    keep_idx = sorted(set(keep_idx))
    qqc.matrices.peptide = peptide_mat[keep_idx]
    qqc.matrices.peptide_rows = rows[keep_idx] if rows is not None else None
    qqc.matrices.peptide_protein_map = ppm.iloc[keep_idx].reset_index(drop=True)

    return qqc


def filt_mat_cr(mat: np.ndarray, pct_r: float, pct_c: float) -> np.ndarray:
    """Filter matrix rows/columns by maximum fraction of missing values."""
    # Filter columns
    na_frac_c = np.sum(np.isnan(mat), axis=0) / mat.shape[0]
    mat = mat[:, na_frac_c <= pct_c]
    # Filter rows
    na_frac_r = np.sum(np.isnan(mat), axis=1) / mat.shape[1]
    mat = mat[na_frac_r <= pct_r]
    return mat
