"""Cell membrane permeability classification — mirrors PermeabilityClassifyer.R."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from quantqc.utils import proc_fasta


def _get_extdata(subpath: str) -> str:
    """Resolve a path inside inst/extdata of the QuantQC package."""
    here = Path(__file__).resolve().parent.parent.parent
    candidate = here / "inst" / "extdata" / subpath
    if candidate.exists():
        return str(candidate)
    raise FileNotFoundError(f"Cannot find extdata file: {subpath}")


def find_permeable_cells(
    mat: np.ndarray,
    row_names: np.ndarray,
    species: str = "Human",
) -> np.ndarray:
    """Predict cell membrane permeability using XGBoost.

    Parameters
    ----------
    mat : array (proteins × cells)
        Protein abundance matrix with UniProt accessions as row labels.
    row_names : array
        Row labels (UniProt accessions).
    species : str
        'Human' or 'Mouse'.

    Returns
    -------
    predictions : array
        Probability of permeability per cell.
    """
    try:
        import xgboost as xgb
    except ImportError:
        raise ImportError("xgboost is required. Install with: pip install xgboost")

    if species not in ("Human", "Mouse"):
        raise ValueError("Species not supported. Options: 'Human', 'Mouse'")

    # Row-normalize by SD
    row_sd = np.nanstd(mat, axis=1, ddof=1, keepdims=True)
    row_sd[row_sd == 0] = 1
    mat = mat / row_sd

    # Load training data
    training_path = _get_extdata("PermeableClassification.csv")
    training = pd.read_csv(training_path, index_col=0)
    labels = training["train"].values
    training_mat = training.drop(columns=["train"]).values
    gene_cols = training.drop(columns=["train"]).columns.values

    # FASTA mapping: accession → gene
    fasta_path = _get_extdata(f"{species}.fasta")
    fasta_map = proc_fasta(fasta_path)
    if species == "Mouse":
        fasta_map["split_gene"] = fasta_map["split_gene"].str.upper()

    # Filter to overlapping genes
    fasta_map = fasta_map[fasta_map["split_gene"].isin(gene_cols)]
    fasta_map = fasta_map[fasta_map["split_prot"].isin(row_names)]

    # Subset training data to available genes
    gene_mask = np.isin(gene_cols, fasta_map["split_gene"].values)
    training_mat = training_mat[:, gene_mask]

    # Subset and reorder input matrix
    prot_to_gene = dict(zip(fasta_map["split_prot"], fasta_map["split_gene"]))
    keep_prots = [p for p in row_names if p in prot_to_gene]
    row_mask = np.isin(row_names, keep_prots)
    mat_sub = mat[row_mask]
    gene_names = np.array([prot_to_gene[p] for p in row_names[row_mask]])

    # Reorder to match training columns
    available_genes = gene_cols[gene_mask]
    gene_order = {g: i for i, g in enumerate(gene_names)}
    reorder = [gene_order[g] for g in available_genes if g in gene_order]
    mat_sub = mat_sub[reorder].T  # cells × genes

    # Column-normalize
    col_sd = np.nanstd(mat_sub, axis=0, ddof=1)
    col_sd[col_sd == 0] = 1
    mat_sub = mat_sub / col_sd

    # Train model
    dtrain = xgb.DMatrix(data=training_mat, label=labels)
    params = {"objective": "binary:logistic", "eval_metric": "logloss"}
    model = xgb.train(params, dtrain, num_boost_round=100)

    # Predict
    dtest = xgb.DMatrix(data=mat_sub)
    predictions = model.predict(dtest)
    return predictions
