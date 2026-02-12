"""Core data structures for QuantQC — Python dataclass equivalents of the R S4 classes."""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Any

import numpy as np
import pandas as pd


@dataclass
class MatricesDDA:
    """Peptide and protein matrices for TMT (DDA) workflows."""
    peptide: np.ndarray = field(default_factory=lambda: np.empty((0, 0)))
    protein: np.ndarray = field(default_factory=lambda: np.empty((0, 0)))
    protein_abs: np.ndarray = field(default_factory=lambda: np.empty((0, 0)))
    protein_imputed: np.ndarray = field(default_factory=lambda: np.empty((0, 0)))
    peptide_protein_map: pd.DataFrame = field(default_factory=lambda: pd.DataFrame(columns=["Protein", "seqcharge"]))

    # Row / column labels stored alongside the numpy matrices
    peptide_rows: np.ndarray | None = None
    peptide_cols: np.ndarray | None = None
    protein_rows: np.ndarray | None = None
    protein_cols: np.ndarray | None = None


@dataclass
class MatricesDIA:
    """Peptide matrices with quality masks for mTRAQ (DIA) workflows."""
    peptide: np.ndarray = field(default_factory=lambda: np.empty((0, 0)))
    peptide_mask: np.ndarray = field(default_factory=lambda: np.empty((0, 0), dtype=bool))
    protein: np.ndarray = field(default_factory=lambda: np.empty((0, 0)))
    protein_mask: np.ndarray = field(default_factory=lambda: np.empty((0, 0), dtype=bool))
    protein_abs: np.ndarray = field(default_factory=lambda: np.empty((0, 0)))
    protein_imputed: np.ndarray = field(default_factory=lambda: np.empty((0, 0)))
    peptide_protein_map: pd.DataFrame = field(default_factory=lambda: pd.DataFrame(columns=["Protein", "seqcharge"]))

    peptide_rows: np.ndarray | None = None
    peptide_cols: np.ndarray | None = None
    protein_rows: np.ndarray | None = None
    protein_cols: np.ndarray | None = None


@dataclass
class MatricesMiceotopes:
    """Heavy/Light isotope quantification matrices."""
    Raw_H: np.ndarray = field(default_factory=lambda: np.empty((0, 0)))
    Raw_L: np.ndarray = field(default_factory=lambda: np.empty((0, 0)))
    HovL_pep: np.ndarray = field(default_factory=lambda: np.empty((0, 0)))
    Beta_pep: np.ndarray = field(default_factory=lambda: np.empty((0, 0)))
    Alpha_pep: np.ndarray = field(default_factory=lambda: np.empty((0, 0)))
    HovL_prot: np.ndarray = field(default_factory=lambda: np.empty((0, 0)))
    Beta_prot: np.ndarray = field(default_factory=lambda: np.empty((0, 0)))
    Alpha_prot: np.ndarray = field(default_factory=lambda: np.empty((0, 0)))
    peptide_protein_map: pd.DataFrame = field(default_factory=lambda: pd.DataFrame(columns=["Protein", "seqcharge"]))

    rows: np.ndarray | None = None
    cols: np.ndarray | None = None


@dataclass
class QQC:
    """Primary QuantQC object — holds all data for a single-cell proteomics experiment.

    Mirrors the R S4 class ``QQC`` with Python dataclasses.
    """
    ms_type: str = ""                           # 'DDA', 'DIA', or 'DIA_C'
    plex: int = 0
    raw_data: pd.DataFrame = field(default_factory=pd.DataFrame)

    matrices: MatricesDDA | MatricesDIA | None = None

    cellenone_meta: pd.DataFrame = field(default_factory=pd.DataFrame)
    meta_data: pd.DataFrame = field(default_factory=pd.DataFrame)
    run_order_statistics: list = field(default_factory=list)
    pep_cor: list = field(default_factory=list)
    neg_ctrl_info: pd.DataFrame = field(default_factory=pd.DataFrame)
    lc_batch_deviations: list = field(default_factory=list)
    reductions: dict[str, pd.DataFrame] = field(default_factory=dict)
    misc: dict[str, Any] = field(default_factory=dict)

    miceotopes: MatricesMiceotopes | None = None


# ---------------------------------------------------------------------------
# Helpers for working with labeled numpy matrices
# ---------------------------------------------------------------------------

def mat_to_df(mat: np.ndarray, rows: np.ndarray | None, cols: np.ndarray | None) -> pd.DataFrame:
    """Convert a labeled numpy matrix back to a pandas DataFrame."""
    return pd.DataFrame(mat, index=rows, columns=cols)


def df_to_mat(df: pd.DataFrame) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Convert a DataFrame to (matrix, row_labels, col_labels)."""
    return df.values.astype(np.float64), np.asarray(df.index), np.asarray(df.columns)
