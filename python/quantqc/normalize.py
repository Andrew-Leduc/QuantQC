"""Normalization, imputation, batch correction, and protein quantification.

Mirrors R/NormalizeAndTransform.R.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import polars as pl
from scipy.interpolate import UnivariateSpline

from quantqc.core import QQC, MatricesDDA, MatricesDIA, df_to_mat
from quantqc.utils import (
    normalize_reference_vector,
    normalize_reference_vector_log,
    normalize,
    normalize_log,
)


# ---------------------------------------------------------------------------
# Cell × Peptide matrix construction
# ---------------------------------------------------------------------------

def cell_x_peptide(qqc: QQC, tq_val: float = 1.0, ch_q_val: float = 1.0) -> QQC:
    """Build cell × peptide matrix — dispatches DDA / DIA."""
    if qqc.ms_type in ("DIA", "DIA_C"):
        return _cell_x_peptide_dia(qqc, tq_val, ch_q_val)
    return _tmt_reference_channel_norm(qqc)


def _tmt_reference_channel_norm(qqc: QQC) -> QQC:
    """TMT reference channel normalization (DDA)."""
    plex = qqc.misc["plex"]
    data = qqc.raw_data.copy()

    if plex == 14:
        ri_cols = [f"Reporter.intensity.{i}" for i in range(2, 19)]
        ref_col = ri_cols[0]
        sc_cols = [f"Reporter.intensity.{i}" for i in range(5, 19)]
    elif plex in (29, 32):
        if plex == 29:
            ri_cols = [f"Reporter.intensity.{i}" for i in range(2, 33)]
            ref_col = ri_cols[0]
        else:
            ri_cols = [f"Reporter.intensity.{i}" for i in range(1, 33)]
            ref_col = None  # no ref channel norm for 32plex
        sc_cols = [f"Reporter.intensity.{i}" for i in range(4, 33)]
    else:
        raise ValueError(f"Unsupported plex: {plex}")

    # Normalise by reference channel
    if ref_col is not None:
        for c in ri_cols:
            data[c] = data[c] / data[ref_col]
        if plex == 14:
            data = data[data[ri_cols].max(axis=1) < 2.5]

    # Build long then pivot wide
    id_cols = ["seqcharge", "Leading.razor.protein", "Raw.file", "Well", "plate"]
    keep_cols = id_cols + sc_cols
    long = data[keep_cols].melt(id_vars=id_cols, var_name="variable", value_name="value")
    long["ID"] = long["Well"].astype(str) + long["plate"].astype(str) + long["variable"].astype(str)

    wide = long.pivot_table(index=["Leading.razor.protein", "seqcharge"], columns="ID", values="value", aggfunc="first")
    wide = wide.reset_index()
    wide = wide.drop_duplicates(subset="seqcharge", keep="first")

    prot_pep_map = wide[["Leading.razor.protein", "seqcharge"]].rename(columns={"Leading.razor.protein": "Protein"})
    mat = wide.drop(columns=["Leading.razor.protein", "seqcharge"]).values.astype(np.float64)
    cols = np.array(wide.drop(columns=["Leading.razor.protein", "seqcharge"]).columns)
    rows = wide["seqcharge"].values

    matrices = MatricesDDA(
        peptide=mat, peptide_protein_map=prot_pep_map.reset_index(drop=True),
        peptide_rows=rows, peptide_cols=cols,
    )
    qqc.matrices = matrices
    return qqc


def _cell_x_peptide_dia(qqc: QQC, tq_val: float, ch_q_val: float, carrier_norm: bool = True) -> QQC:
    """mTRAQ cell × peptide matrix construction (DIA)."""
    raw = qqc.raw_data.copy()
    plex = qqc.misc["plex"]
    ms_type = qqc.ms_type

    if plex == 2 and ms_type == "DIA_C":
        plex_used = [0, 4, 8]
    elif plex == 2 and ms_type == "DIA":
        plex_used = [0, 4]
    elif plex == 3:
        plex_used = [0, 4, 8]
    else:
        raise ValueError("plex not valid")

    if not carrier_norm:
        plex_used = [0, 4]

    raw["plex"] = pd.to_numeric(raw["plex"], errors="coerce")
    raw = raw[raw["plex"].isin(plex_used)]

    filt = raw[raw["Channel.Q.Value"] < ch_q_val]
    filt = filt[filt["Translated.Q.Value"] < tq_val] if "Translated.Q.Value" in filt.columns else filt

    # Filtered pivot
    wide_filt = filt.pivot_table(
        index=["Protein.Group", "seqcharge"], columns="File.Name",
        values="Ms1.Area", aggfunc="first",
    ).reset_index()
    mat_filt = wide_filt.drop(columns=["Protein.Group", "seqcharge"]).values.astype(np.float64)
    mat_filt[mat_filt == 0] = np.nan

    # Unfiltered pivot for mask — align columns to filtered pivot
    filt_cell_cols = list(wide_filt.drop(columns=["Protein.Group", "seqcharge"]).columns)
    wide_nf = raw.pivot_table(
        index=["Protein.Group", "seqcharge"], columns="File.Name",
        values="Ms1.Area", aggfunc="first",
    ).reset_index()
    wide_nf = wide_nf[wide_nf["seqcharge"].isin(wide_filt["seqcharge"])]
    # Ensure mask has same columns as filtered matrix
    for c in filt_cell_cols:
        if c not in wide_nf.columns:
            wide_nf[c] = np.nan
    mat_nf = wide_nf[["Protein.Group", "seqcharge"] + filt_cell_cols].drop(
        columns=["Protein.Group", "seqcharge"]).values.astype(np.float64)
    mat_nf[mat_nf == 0] = np.nan
    pep_mask = ~np.isnan(mat_nf)

    # Carrier normalization
    if ms_type == "DIA_C" and carrier_norm:
        cols_filt = np.array(wide_filt.drop(columns=["Protein.Group", "seqcharge"]).columns)
        mat_filt, cols_filt, pep_mask = _dia_carrier_norm_matrix(mat_filt, cols_filt, 8, plex_used, pep_mask)
    else:
        cols_filt = np.array(wide_filt.drop(columns=["Protein.Group", "seqcharge"]).columns)

    prot_pep_map = wide_filt[["seqcharge", "Protein.Group"]].rename(columns={"Protein.Group": "Protein"})
    rows = wide_filt["seqcharge"].values

    matrices = MatricesDIA(
        peptide=mat_filt, peptide_mask=pep_mask,
        peptide_protein_map=prot_pep_map.reset_index(drop=True),
        peptide_rows=rows, peptide_cols=cols_filt,
    )
    qqc.matrices = matrices
    qqc.misc["ChQ"] = ch_q_val
    return qqc


def _dia_carrier_norm_matrix(mat, cols, carrier_ch, plex_used, mask):
    """Carrier normalization — divide cell signals by carrier, remove bad sets."""
    carrier_suffix = str(carrier_ch)
    carrier_mask = np.array([str(c).endswith(carrier_suffix) for c in cols])
    carrier_indices = np.where(carrier_mask)[0]

    remove_cols = set()
    for ci in carrier_indices:
        carrier_col = cols[ci]
        prefix = carrier_col[:-1]
        flag = False
        for p in plex_used:
            sc_col = f"{prefix}{p}"
            if sc_col in cols:
                sc_idx = np.where(cols == sc_col)[0]
                if len(sc_idx):
                    ratio = np.nanmedian(mat[:, sc_idx[0]] / mat[:, ci])
                    if ratio > 1:
                        flag = True
        if flag:
            remove_cols.add(carrier_col)
            for p in plex_used:
                remove_cols.add(f"{prefix}{p}")

    # Divide by carrier
    for ci in carrier_indices:
        carrier_col = cols[ci]
        if carrier_col in remove_cols:
            continue
        prefix = carrier_col[:-1]
        for p in plex_used:
            sc_col = f"{prefix}{p}"
            if sc_col in cols:
                sc_idx = np.where(cols == sc_col)[0]
                if len(sc_idx):
                    mat[:, sc_idx[0]] = mat[:, sc_idx[0]] / mat[:, ci]

    # Remove carrier columns and flagged columns
    all_remove = set()
    all_remove.update(remove_cols)
    for ci in carrier_indices:
        all_remove.add(cols[ci])
    keep_mask = np.array([c not in all_remove for c in cols])
    mat = mat[:, keep_mask]
    mask = mask[:, keep_mask] if mask.shape[1] == len(cols) else mask
    cols = cols[keep_mask]

    mat[np.isinf(mat)] = np.nan
    mat[mat == 0] = np.nan
    return mat, cols, mask


# ---------------------------------------------------------------------------
# Normalization (public wrappers re-exported from utils for convenience)
# ---------------------------------------------------------------------------
# normalize_reference_vector, normalize_reference_vector_log, normalize, normalize_log
# are all importable from quantqc.utils — we re-export here for API compatibility.


# ---------------------------------------------------------------------------
# Imputation
# ---------------------------------------------------------------------------

def knn_impute(qqc: QQC, k: int = 3) -> QQC:
    """K-nearest neighbors imputation on protein matrix."""
    data = qqc.matrices.protein.copy()
    imputed = data.copy()

    # Euclidean distance between columns (cells)
    from scipy.spatial.distance import pdist, squareform
    # Transpose: each column is a sample → distance between samples
    # We need to handle NaN: use pairwise complete observations
    n_cols = data.shape[1]
    dist_mat = np.full((n_cols, n_cols), np.inf)
    for i in range(n_cols):
        for j in range(i, n_cols):
            valid = ~np.isnan(data[:, i]) & ~np.isnan(data[:, j])
            if valid.sum() > 0:
                d = np.sqrt(np.sum((data[valid, i] - data[valid, j]) ** 2))
                dist_mat[i, j] = d
                dist_mat[j, i] = d
            dist_mat[i, i] = 0

    for col_idx in range(n_cols):
        vec = data[:, col_idx].copy()
        na_rows = np.where(np.isnan(vec))[0]
        order = np.argsort(dist_mat[col_idx])

        for row_idx in na_rows:
            # Find closest columns that have a value at this row
            available = [o for o in order if not np.isnan(data[row_idx, o])]
            if len(available) > k:
                available = available[:k]
            if available:
                vec[row_idx] = np.mean([data[row_idx, o] for o in available])

        imputed[:, col_idx] = vec

    imputed = normalize_reference_vector_log(imputed)
    qqc.matrices.protein_imputed = imputed
    return qqc


def min_value_impute(qqc: QQC) -> QQC:
    """Replace NaN with row minimum."""
    data = qqc.matrices.protein.copy()
    for i in range(data.shape[0]):
        row = data[i]
        min_val = np.nanmin(row)
        data[i, np.isnan(row)] = min_val
    qqc.matrices.protein_imputed = data
    return qqc


# ---------------------------------------------------------------------------
# Protein collapse
# ---------------------------------------------------------------------------

def collapse_to_protein(qqc: QQC, opt: int = 1, lc_correct: bool = False, norm: str = "ref") -> QQC:
    """Collapse peptide-level data to protein level.

    opt=1: median of normalized peptides. opt=2: MaxLFQ (requires directlfq).
    """
    pep_mat = qqc.matrices.peptide
    ppm = qqc.matrices.peptide_protein_map
    rows = qqc.matrices.peptide_rows
    cols = qqc.matrices.peptide_cols

    # --- Absolute abundances (always computed) ---
    abs_df = pd.DataFrame(pep_mat, columns=cols)
    abs_df["Protein"] = ppm["Protein"].values
    abs_long = abs_df.melt(id_vars=["Protein"], var_name="variable", value_name="value")
    abs_prot = abs_long.groupby(["Protein", "variable"])["value"].median().reset_index()
    abs_wide = abs_prot.pivot(index="Protein", columns="variable", values="value")
    qqc.matrices.protein_abs = abs_wide.values.astype(np.float64)

    if opt == 1:
        # Normalize peptide data
        if qqc.ms_type == "DDA":
            norm_pep = normalize_reference_vector(pep_mat, log=True)
            qqc.matrices.peptide = norm_pep
        else:
            if lc_correct:
                qqc = lc_batch_correct(qqc)
                norm_pep = qqc.matrices.peptide
            else:
                norm_pep = normalize_reference_vector(pep_mat, log=True)
                qqc.matrices.peptide = norm_pep

        # Filter peptides with < 3 observations
        obs_count = np.sum(~np.isnan(norm_pep), axis=1)
        keep = obs_count > 3
        norm_pep = norm_pep[keep]
        ppm_filt = ppm[keep].reset_index(drop=True)

        # Collapse: median per protein per cell
        df = pd.DataFrame(norm_pep, columns=cols)
        df["Protein"] = ppm_filt["Protein"].values
        long = df.melt(id_vars=["Protein"], var_name="variable", value_name="value")
        prot_data = long.groupby(["Protein", "variable"])["value"].median().reset_index()
        prot_wide = prot_data.pivot(index="Protein", columns="variable", values="value")

        # Re-normalise at protein level
        prot_mat = prot_wide.values.astype(np.float64)
        if norm == "std":
            prot_mat = normalize_log(prot_mat)
        else:
            prot_mat = normalize_reference_vector_log(prot_mat)

        qqc.matrices.protein = prot_mat
        qqc.matrices.protein_rows = np.array(prot_wide.index)
        qqc.matrices.protein_cols = np.array(prot_wide.columns)

        # Protein mask for DIA
        if qqc.ms_type in ("DIA", "DIA_C") and hasattr(qqc.matrices, "peptide_mask"):
            mask_df = pd.DataFrame(qqc.matrices.peptide_mask.astype(float), columns=qqc.matrices.peptide_cols)
            mask_df["Protein"] = qqc.matrices.peptide_protein_map["Protein"].values
            mask_long = mask_df.melt(id_vars=["Protein"], var_name="variable", value_name="value")
            mask_prot = mask_long.groupby(["Protein", "variable"])["value"].sum().reset_index()
            mask_prot["value"] = mask_prot["value"] > 0
            mask_wide = mask_prot.pivot(index="Protein", columns="variable", values="value")
            qqc.matrices.protein_mask = mask_wide.values

    elif opt == 2:
        raise NotImplementedError("MaxLFQ (opt=2) requires directlfq — not yet integrated")

    return qqc


# ---------------------------------------------------------------------------
# Batch correction
# ---------------------------------------------------------------------------

def batch_correct(qqc: QQC, labels: bool = True, run: bool = True,
                  batch: bool = False, norm: str = "ref") -> QQC:
    """Batch effect removal using limma-style removeBatchEffect or ComBat.

    Requires the ``combat`` Python package for ComBat, or uses a simple
    linear model approach as fallback.
    """
    if run and batch:
        raise ValueError("Cannot correct on both run and batch simultaneously")

    meta = qqc.meta_data
    prot_imp = qqc.matrices.protein_imputed.copy()
    prot_raw = qqc.matrices.protein.copy()
    cols = qqc.matrices.protein_cols

    # Align metadata to matrix columns
    meta_aligned = meta[meta["ID"].isin(cols)].copy()
    meta_aligned = meta_aligned.set_index("ID").loc[cols].reset_index()

    # Simple batch correction via mean-centering per batch
    def _remove_batch(mat, batch_var):
        corrected = mat.copy()
        for b in batch_var.unique():
            mask = (batch_var == b).values
            if mask.sum() > 0:
                batch_mean = np.nanmean(corrected[:, mask], axis=1, keepdims=True)
                global_mean = np.nanmean(corrected, axis=1, keepdims=True)
                corrected[:, mask] = corrected[:, mask] - batch_mean + global_mean
        return corrected

    corrected = prot_imp
    if labels and "label" in meta_aligned.columns:
        corrected = _remove_batch(corrected, meta_aligned["label"])
    if run and "injectWell" in meta_aligned.columns:
        corrected = _remove_batch(corrected, meta_aligned["injectWell"])
    if batch and "LCMS_Batch" in meta_aligned.columns:
        corrected = _remove_batch(corrected, meta_aligned["LCMS_Batch"])

    # Re-normalize
    if norm == "std":
        corrected = normalize_log(corrected)
    else:
        corrected = normalize_reference_vector_log(corrected)

    # Restore NAs from original
    corrected_noimp = corrected.copy()
    corrected_noimp[np.isnan(prot_raw)] = np.nan

    qqc.matrices.protein_imputed = corrected
    qqc.matrices.protein = corrected_noimp
    return qqc


def lc_batch_correct(qqc: QQC) -> QQC:
    """Run-order drift correction via spline smoothing."""
    pep_norm = qqc.matrices.peptide.copy()
    pep_norm[pep_norm == 0] = np.nan
    pep_norm = normalize_reference_vector(pep_norm, log=True)

    order_df = qqc.meta_data[qqc.meta_data["ID"].isin(qqc.matrices.peptide_cols)]
    orders = order_df.set_index("ID").loc[qqc.matrices.peptide_cols]["Order"].values.astype(float)

    rsq_save = []
    plot_residuals = []

    for i in range(pep_norm.shape[0]):
        row = pep_norm[i]
        valid = ~np.isnan(row)
        n_valid = valid.sum()

        if n_valid <= 29:
            continue

        df_spline = min(max(n_valid // 11, 2), 20)
        x = orders[valid]
        y = row[valid]

        try:
            spline = UnivariateSpline(x, y, k=3, s=len(x) * df_spline)
            predicted = spline(x)
        except Exception:
            continue

        ss_res = np.sum((y - predicted) ** 2)
        ss_tot = np.sum((y - np.mean(predicted)) ** 2)
        r_sq = 1 - ss_res / ss_tot if ss_tot > 0 else 0

        rsq_save.append(r_sq)

        if r_sq > 0.1:
            sort_idx = np.argsort(x)
            pep_norm[i, valid] = y[sort_idx] - predicted[sort_idx]

            if r_sq > 0.5 and len(plot_residuals) < 10:
                plot_residuals.append({
                    "order": x, "data": y, "model": predicted, "score": str(round(r_sq, 3)),
                })

    qqc.matrices.peptide = pep_norm
    qqc.lc_batch_deviations = [rsq_save, plot_residuals]
    return qqc
