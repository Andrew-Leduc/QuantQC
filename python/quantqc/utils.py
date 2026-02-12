"""Utility functions — mirrors misc.R helpers."""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Normalisation helpers (used by many modules)
# ---------------------------------------------------------------------------

def normalize(mat: np.ndarray, log: bool = False) -> np.ndarray:
    """Column-median then row-mean normalization in linear space."""
    mat = mat.astype(np.float64).copy()
    mat[mat == 0] = np.nan
    # Column normalisation: divide each column by its median
    col_med = np.nanmedian(mat, axis=0)
    col_med[col_med == 0] = np.nan
    mat = mat / col_med[np.newaxis, :]
    # Row normalisation: divide each row by its mean
    row_mean = np.nanmean(mat, axis=1)
    row_mean[row_mean == 0] = np.nan
    mat = mat / row_mean[:, np.newaxis]
    if log:
        mat = np.log2(mat)
    return mat


def normalize_log(mat: np.ndarray) -> np.ndarray:
    """Column-median then row-mean normalization in log space (additive)."""
    mat = mat.astype(np.float64).copy()
    mat[mat == 0] = np.nan
    col_med = np.nanmedian(mat, axis=0)
    mat = mat - col_med[np.newaxis, :]
    row_mean = np.nanmean(mat, axis=1)
    mat = mat - row_mean[:, np.newaxis]
    return mat


def normalize_reference_vector(dat: np.ndarray, log: bool = False) -> np.ndarray:
    """Reference-median normalization (Normalize_reference_vector in R)."""
    dat = dat.astype(np.float64).copy()
    dat[dat == 0] = np.nan
    ref_vec = np.nanmedian(dat, axis=1)
    for k in range(dat.shape[1]):
        col = dat[:, k]
        ratio = ref_vec / col
        dat[:, k] = col * np.nanmedian(ratio)
    for k in range(dat.shape[0]):
        row = dat[k, :]
        m = np.nanmean(row)
        if m != 0 and not np.isnan(m):
            dat[k, :] = row / m
    if log:
        with np.errstate(divide="ignore"):
            dat = np.log2(dat)
    return dat


def normalize_reference_vector_log(dat: np.ndarray) -> np.ndarray:
    """Reference-median normalization in log space."""
    dat = dat.astype(np.float64).copy()
    ref_vec = np.nanmedian(dat, axis=1)
    for k in range(dat.shape[1]):
        diff = ref_vec - dat[:, k]
        dat[:, k] = dat[:, k] + np.nanmedian(diff)
    for k in range(dat.shape[0]):
        dat[k, :] = dat[k, :] - np.nanmean(dat[k, :])
    return dat


# ---------------------------------------------------------------------------
# CV helpers
# ---------------------------------------------------------------------------

def cv(x: np.ndarray) -> float:
    """Coefficient of variation (sd / mean), ignoring NaN."""
    return float(np.nanstd(x, ddof=1) / np.nanmean(x))


def fast_cv(df: pd.DataFrame, protein_col: str = "Protein") -> pd.DataFrame:
    """Grouped CV calculation — equivalent to R fast_cv using data.table.

    Expects a long-format DataFrame with columns: ``variable``, ``Protein``, ``value``.
    Returns the same frame with an added ``cvq`` column.
    """
    valid = df.dropna(subset=["value"])
    grouped = valid.groupby(["variable", protein_col])["value"]
    stats = grouped.agg(["std", "mean"]).reset_index()
    stats["cvq"] = stats["std"] / stats["mean"]
    result = df.merge(stats[["variable", protein_col, "cvq"]], on=["variable", protein_col], how="left")
    return result


# ---------------------------------------------------------------------------
# Protein / peptide string helpers
# ---------------------------------------------------------------------------

def extract_accession(input_string: str) -> str:
    """Extract UniProt accession between pipe characters ``|ACC|``."""
    m = re.search(r"\|([^|]+)\|", input_string)
    return m.group(1) if m else input_string


def contains_missed_cleaved(pep_list: pd.Series) -> pd.Series:
    """Count missed cleavage sites (internal K or R, ignoring C-terminal)."""
    trimmed = pep_list.str[:-1]
    has_k = trimmed.str.contains("K", fixed=True).astype(int)
    has_r = trimmed.str.contains("R", fixed=True).astype(int)
    return has_k + has_r


def create_peptide_vector(protein_sequence: str, peptide: str) -> np.ndarray:
    """Return binary vector marking peptide location in protein sequence."""
    vec = np.zeros(len(protein_sequence), dtype=int)
    start = protein_sequence.find(peptide)
    if start >= 0:
        vec[start : start + len(peptide)] = 1
    return vec


# ---------------------------------------------------------------------------
# FASTA helpers
# ---------------------------------------------------------------------------

def proc_fasta(path: str | Path) -> pd.DataFrame:
    """Parse FASTA headers → DataFrame with ``split_prot`` and ``split_gene`` columns."""
    from Bio import SeqIO

    prots, genes = [], []
    for record in SeqIO.parse(str(path), "fasta"):
        header = record.description
        if "GN=" not in header:
            continue
        parts = header.split("GN=", 1)
        gene = parts[1].split()[0] if parts[1] else ""
        prot_part = parts[0]
        acc_match = re.search(r"\|([^|]+)\|", prot_part)
        acc = acc_match.group(1) if acc_match else prot_part.strip()
        prots.append(acc)
        genes.append(gene)
    return pd.DataFrame({"split_prot": prots, "split_gene": genes})


# ---------------------------------------------------------------------------
# Summary function for boxplots (mirrors R ``f``)
# ---------------------------------------------------------------------------

def boxplot_summary(y: np.ndarray) -> dict:
    """Return count and median — used for annotating boxplots."""
    return {"label": len(y), "y": float(np.nanmedian(y))}


# ---------------------------------------------------------------------------
# in-set normalization
# ---------------------------------------------------------------------------

def in_set_norm(raw_data: np.ndarray, cols: np.ndarray, cellenone_meta: pd.DataFrame) -> np.ndarray:
    """Within-injection normalization: centre each inject well separately."""
    mat = raw_data.astype(np.float64).copy()
    for well in cellenone_meta["injectWell"].unique():
        well_str = str(well)
        mask = np.array([well_str in str(c) for c in cols])
        if mask.any():
            subset = mat[:, mask]
            row_means = np.nanmean(subset, axis=1, keepdims=True)
            mat[:, mask] = subset - row_means
    return mat
