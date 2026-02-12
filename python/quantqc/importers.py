"""Data importers — MQ, DIA-NN, JMOD, FragPipe → QQC object.

Uses polars for fast file reading and initial filtering, converts to pandas
for the QQC object storage (since downstream code uses pandas/numpy).
"""

from __future__ import annotations

import os
import re
from pathlib import Path

import polars as pl
import pandas as pd
import numpy as np

from quantqc.core import QQC, MatricesDDA, MatricesDIA


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

def _read_single_file(path: Path, columns: list[str], separator: str = "\t") -> pl.DataFrame:
    """Read a single TSV or parquet file, selecting *columns*."""
    if path.suffix.lower() in (".parquet", ".pq"):
        return pl.read_parquet(str(path), columns=columns)
    return pl.read_csv(str(path), separator=separator, columns=columns, infer_schema_length=10000)


def _read_data_single_or_dir(data_path: str | Path, columns: list[str], separator: str = "\t") -> pl.DataFrame:
    """Read a single file (TSV or parquet) or all files in a directory, selecting only *columns*."""
    p = Path(data_path)
    if p.is_dir():
        frames = []
        for f in sorted(p.iterdir()):
            if f.is_file():
                frames.append(_read_single_file(f, columns, separator))
        return pl.concat(frames)
    return _read_single_file(p, columns, separator)


def _clean_protein_name(name: str) -> str:
    """Parse leading razor protein: extract accession from pipe-delimited string, strip isoform."""
    if "|" in name:
        parts = name.split("|")
        if len(parts) >= 3:
            name = parts[1]
    if "-" in name:
        name = name.split("-")[0]
    return name


# ---------------------------------------------------------------------------
# MaxQuant → QQC  (DDA / TMT)
# ---------------------------------------------------------------------------

def mq_to_qqc(data_path: str, linker: str, plex: int, PIF_in: float, PEP_in: float) -> QQC:
    """Import MaxQuant evidence file(s) into a QQC object."""
    linker_df = pd.read_csv(linker)

    ri_numb = {14: 18, 29: 32, 32: 32}[plex]

    ri_cols = [f"Reporter intensity {i}" for i in range(1, ri_numb + 1)]
    base_cols = [
        "Modified sequence", "Intensity", "Retention time", "Charge",
        "Raw file", "PEP", "PIF", "Leading razor protein",
        "Potential contaminant", "Reverse",
    ]
    columns_to_read = base_cols + ri_cols

    data = _read_data_single_or_dir(data_path, columns_to_read)

    # Normalise column names: spaces → dots
    data = data.rename({c: c.replace(" ", ".") for c in data.columns})

    data = data.to_pandas()

    # Filter to runs in the linker
    data = data[data["Raw.file"].isin(linker_df["Run"])]

    linker_df["Order"] = range(1, len(linker_df) + 1)
    data = data.merge(linker_df, left_on="Raw.file", right_on="Run", how="left")

    # Unique precursor ID
    data["seqcharge"] = data["Modified.sequence"].astype(str) + data["Charge"].astype(str)
    data["seqRun"] = data["seqcharge"] + data["Raw.file"].astype(str)
    data = data.drop_duplicates(subset="seqRun", keep="first")

    # Clean protein names
    data["Leading.razor.protein"] = data["Leading.razor.protein"].apply(_clean_protein_name)

    # Filter
    data = data[data["PEP"] < PEP_in]
    data = data[data["PIF"] > PIF_in]
    data = data[data["Reverse"] != "+"]

    return QQC(
        raw_data=data.reset_index(drop=True),
        meta_data=linker_df,
        ms_type="DDA",
        plex=plex,
        misc={"plex": plex},
    )


# ---------------------------------------------------------------------------
# DIA-NN → QQC  (DIA / mTRAQ)
# ---------------------------------------------------------------------------

def diann_to_qqc(data_path: str, linker_path: str, plex: int, carrier: bool = False) -> QQC:
    """Import DIA-NN report into a QQC object."""
    linker_df = pd.read_csv(linker_path)
    linker_df["Order"] = range(1, len(linker_df) + 1)

    columns_to_read = [
        "Genes", "Run", "Lib.PG.Q.Value", "RT", "Precursor.Id",
        "Stripped.Sequence", "Precursor.Mz", "Precursor.Charge",
        "Precursor.Quantity", "Ms1.Area", "Protein.Group",
        "Translated.Q.Value", "Channel.Q.Value",
    ]

    data = _read_data_single_or_dir(data_path, columns_to_read).to_pandas()

    data = data[data["Lib.PG.Q.Value"] < 0.01]
    data = data[data["Run"].isin(linker_df["Run"])]
    data = data.merge(linker_df, on="Run", how="left")

    data["seqcharge"] = data["Stripped.Sequence"].astype(str) + data["Precursor.Charge"].astype(str)
    data = data[data["Protein.Group"] != ""]

    # mTRAQ tag from Precursor.Id (character at position 10)
    data["plex"] = data["Precursor.Id"].str[9]

    # Unique cell ID
    data["ID"] = data["Well"].astype(str) + data["plate"].astype(str) + "." + data["plex"].astype(str)
    data["File.Name"] = data["ID"]

    # Remove duplicates
    data["uq"] = data["File.Name"].astype(str) + data["Protein.Group"].astype(str) + data["seqcharge"]
    data = data.drop_duplicates(subset="uq", keep="first")
    data = data.drop(columns=["uq"])

    ms_type = "DIA_C" if carrier else "DIA"
    return QQC(
        raw_data=data.reset_index(drop=True),
        meta_data=linker_df,
        ms_type=ms_type,
        plex=plex,
        misc={"plex": plex},
    )


# ---------------------------------------------------------------------------
# JMOD → QQC  (DIA)
# ---------------------------------------------------------------------------

def jmod_to_qqc(data_path: str, linker_path: str, plex: int, carrier: bool = False) -> QQC:
    """Import JMOD pipeline output into a QQC object."""
    linker_df = pd.read_csv(linker_path)
    linker_df["Order"] = range(1, len(linker_df) + 1)

    jmod_cols = ["file_name", "seq", "stripped_seq", "z", "plex_Area", "protein", "channel"]

    p = Path(data_path)
    if p.suffix.lower() in (".parquet", ".pq"):
        data = pl.read_parquet(str(p), columns=jmod_cols).to_pandas()
    else:
        data = pl.read_csv(str(p), separator="\t", columns=jmod_cols, infer_schema_length=10000).to_pandas()

    data.columns = [
        "Run", "Precursor.Id", "Stripped.Sequence",
        "Precursor.Charge",
        "Ms1.Area", "Protein.Group", "plex",
    ]

    data["RT"] = 0
    data["Precursor.Quantity"] = 0
    data["Precursor.Mz"] = 0
    data["Channel.Q.Value"] = 0

    data = data[data["Run"].isin(linker_df["Run"])]
    data = data.merge(linker_df, on="Run", how="left")

    data["seqcharge"] = data["Stripped.Sequence"].astype(str) + data["Precursor.Charge"].astype(str)
    data = data[data["Protein.Group"] != ""]

    data["ID"] = data["Well"].astype(str) + data["plate"].astype(str) + "." + data["plex"].astype(str)
    data["File.Name"] = data["ID"]

    ms_type = "DIA_C" if carrier else "DIA"
    return QQC(
        raw_data=data.reset_index(drop=True),
        meta_data=linker_df,
        ms_type=ms_type,
        plex=plex,
        misc={"plex": plex},
    )


# ---------------------------------------------------------------------------
# FragPipe → QQC  (DDA)
# ---------------------------------------------------------------------------

def fragpipe_to_qqc(data_path: str, linker: str, plex: int, PIF_in: float) -> QQC:
    """Import FragPipe PSM results into a QQC object."""
    p = Path(data_path)
    folders = [d for d in sorted(p.iterdir()) if d.is_dir()]

    frames: list[pd.DataFrame] = []
    tmt_grab: list[str] | None = None

    for folder in folders:
        psm_file = folder / "psm.tsv"
        if not psm_file.exists():
            continue
        df = pd.read_csv(psm_file, sep="\t")
        if tmt_grab is None:
            tmt_grab = [c for c in df.columns if c.startswith("Intensity.")]
        base = ["Peptide", "Intensity", "Apex.Retention.Time", "Charge", "Spectrum.File", "Purity", "Protein.ID"]
        df = df[base + tmt_grab]
        # Rename TMT columns to short tag names
        rename_map = {}
        for c in tmt_grab:
            tag = c.rsplit("_", 1)[-1]
            rename_map[c] = tag
        df = df.rename(columns=rename_map)
        frames.append(df)

    data = pd.concat(frames, ignore_index=True)

    # Extract raw file name from Spectrum.File path
    data["Spectrum.File"] = data["Spectrum.File"].apply(lambda x: x.replace("\\", "/").split("/")[-2] if "/" in x.replace("\\", "/") else x)

    linker_df = pd.read_csv(linker)

    if plex == 32:
        map_fp = [
            "126", "127N", "127C", "128N", "128C", "129N", "129C",
            "130N", "130C", "131N", "131C", "132N", "132C",
            "133N", "133C", "134N", "134C", "135N", "127D",
            "128ND", "128CD", "129ND", "129CD", "130ND", "130CD",
            "131ND", "131CD", "132ND", "132CD", "133ND", "133CD",
            "134ND", "134CD", "135ND", "135CD",
        ]
        map_mq = [f"Reporter.ion.{i}" for i in [
            "1", "2", "3", "5", "7", "9", "11", "13", "15", "17",
            "19", "21", "23", "25", "27", "29", "31", "33", "4",
            "6", "8", "10", "12", "14", "16", "18", "20", "22",
            "24", "26", "28", "30", "32", "34", "35",
        ]]
        columns_renamed = [
            "Modified.sequence", "Intensity", "Retention.time", "Charge",
            "Raw.file", "PIF", "Leading.razor.protein",
        ] + map_mq
        data.columns = columns_renamed

    data = data[data["Raw.file"].isin(linker_df["Run"])]
    linker_df["Order"] = range(1, len(linker_df) + 1)
    data = data.merge(linker_df, left_on="Raw.file", right_on="Run", how="left")

    data["seqcharge"] = data["Modified.sequence"].astype(str) + data["Charge"].astype(str)
    data["seqRun"] = data["seqcharge"] + data["Raw.file"].astype(str)
    data = data.drop_duplicates(subset="seqRun", keep="first")

    data = data[data["PIF"] > PIF_in]

    return QQC(
        raw_data=data.reset_index(drop=True),
        meta_data=linker_df,
        ms_type="DDA",
        plex=plex,
        misc={"plex": plex},
    )
