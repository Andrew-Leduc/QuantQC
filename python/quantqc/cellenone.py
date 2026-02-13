"""CellenONE metadata parsing — maps cell sorter output to MS data.

Mirrors R/CellenONE_metadata_map.R: analyzeCellenONE_TMT, analyzeCellenONE_mTRAQ,
link_cellenONE_Raw, PlotSlideLayout_*.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.neighbors import NearestNeighbors

from quantqc.core import QQC
from quantqc.plotting import dot_plot_style


# ---------------------------------------------------------------------------
# .fld file parser
# ---------------------------------------------------------------------------

def _parse_fld(path: str | Path, skip: int = 22) -> pd.DataFrame:
    """Parse a CellenONE .fld file into a DataFrame with position, well, volume, field."""
    import re

    rows: list[dict] = []
    with open(path, "r", encoding="latin-1") as fh:
        lines = fh.readlines()

    current_field: str | None = None
    for line in lines[skip:]:
        parts = line.rstrip("\n").split("\t")
        if not parts:
            continue
        pos = parts[0]
        # Check for field header like [0, ...]
        m = re.match(r"\[(\d+)", pos)
        if m:
            current_field = str(int(m.group(1)) + 1)
            continue
        well = parts[1] if len(parts) > 1 else ""
        vol = parts[2] if len(parts) > 2 else ""
        if well == "":
            continue
        # Parse x/y from position string "y/x"
        if "/" in pos:
            y_str, x_str = pos.split("/", 1)
        else:
            y_str, x_str = pos, ""
        rows.append({
            "position": pos,
            "well": well,
            "volume": vol,
            "field": current_field,
            "yPos": float(y_str) if y_str else np.nan,
            "xPos": float(x_str) if x_str else np.nan,
        })

    return pd.DataFrame(rows)


def _get_extdata(subpath: str) -> str:
    """Resolve a path inside inst/extdata of the QuantQC package."""
    here = Path(__file__).resolve().parent.parent.parent  # repo root
    candidate = here / "inst" / "extdata" / subpath
    if candidate.exists():
        return str(candidate)
    raise FileNotFoundError(f"Cannot find extdata file: {subpath}")


# ---------------------------------------------------------------------------
# TMT CellenONE parser
# ---------------------------------------------------------------------------

def _analyze_cellenone_tmt(cells_file: dict[str, str], plex: int) -> pd.DataFrame:
    """Parse CellenONE isolation files for TMT (DDA) experiments."""
    # Read and concat all isolation files
    frames = []
    for name, path in cells_file.items():
        df = pd.read_csv(path, sep="\t")
        df["condition"] = name
        frames.append(df)
    cells = pd.concat(frames, ignore_index=True)

    # Filter out Transmission/Blue/Green rows (R: fill(2:7, .direction="up"))
    green_df = None
    if "X" in cells.columns:
        cells.loc[cells["X"].astype(str).str.contains("Transmission", na=False), "X"] = np.nan
        # Backward fill (upward) on all columns EXCEPT X, matching R's fill(2:7, .direction="up")
        fill_cols = [c for c in cells.columns if c != "X"]
        cells[fill_cols] = cells[fill_cols].bfill()
        cells = cells.dropna(subset=["XPos"])
        if (cells["X"] == "Blue").any():
            cells.loc[cells["X"] == "Blue", "X"] = np.nan
            cells = cells.dropna(subset=["X"])
        if (cells["X"] == "Green").any():
            idx = list(range(1, len(cells), 2))
            green_df = cells.iloc[idx].copy()
            cells.loc[cells["X"] == "Green", "X"] = np.nan
            cells = cells.dropna(subset=["X"])

    # Load label and pickup .fld files
    plex_dir = {14: "14plex_files", 12: "12plex_files", 29: "29plex_files", 32: "32plex_files"}
    label = _parse_fld(_get_extdata(f"{plex_dir[plex]}/Labels.fld"))
    pickup = _parse_fld(_get_extdata(f"{plex_dir[plex]}/Pickup_mock.fld"))

    # Clean up label well names
    # Strip prefix "1G" style → numeric, and remap 1H/1P wells for 32plex
    well_clean = label["well"].str.lstrip("1").str.rstrip(",")
    # Remove alpha prefix (e.g. "G25" → "25")
    well_clean = well_clean.str.replace(r"^[A-Za-z]+", "", regex=True)
    label["well"] = pd.to_numeric(well_clean, errors="coerce")

    # Map well numbers to TMT tag indices
    unique_wells = sorted(label["well"].dropna().unique())
    well_to_tag = {w: str(i + 1) for i, w in enumerate(unique_wells)}
    label["well"] = label["well"].map(well_to_tag)

    # Normalize Field to integer strings
    cells["Field"] = pd.to_numeric(cells["Field"], errors="coerce")
    cells = cells.dropna(subset=["Field"])
    cells["Field"] = cells["Field"].astype(int).astype(str)

    # Filter to fields present in the cells data
    cell_fields = cells["Field"].unique()
    label = label[label["field"].isin(cell_fields)]
    pickup = pickup[pickup["field"].isin(cell_fields)]

    # Clean pickup well
    pickup["well"] = pickup["well"].str.lstrip("1").str.rstrip(",")
    pickup["well"] = pickup["well"].str.replace(r"^[A-Za-z]", "", regex=True)
    pickup["target"] = 1

    # Create xyf keys for merging
    cells["xyf"] = cells["XPos"].astype(str) + cells["YPos"].astype(str) + cells["Field"].astype(str)
    label["xyf"] = label["xPos"].astype(str) + label["yPos"].astype(str) + label["field"].astype(str)

    iso_lab = cells.merge(label[["xyf", "well"]], on="xyf", how="left")

    # ANN matching for pickup positions
    pickup_xy = pickup[["xPos", "yPos"]].drop_duplicates().reset_index(drop=True)
    iso_xy = iso_lab[["xPos", "yPos"]].values.astype(float)

    # Coordinate adjustments for 32plex
    if plex == 32:
        iso_xy = iso_xy.copy()
        iso_xy[iso_xy[:, 1] > 50, 1] += 10
        iso_xy[iso_xy[:, 1] > 33, 1] += 10
        iso_xy[iso_xy[:, 1] > 17, 1] += 10
        pkp = pickup_xy.values.copy()
        pkp[pkp[:, 1] > 50, 1] += 10
        pkp[pkp[:, 1] > 33, 1] += 10
        pkp[pkp[:, 1] > 17, 1] += 10
    else:
        pkp = pickup_xy.values

    nn = NearestNeighbors(n_neighbors=1).fit(pkp)
    _, indices = nn.kneighbors(iso_xy)
    iso_lab["pickupX"] = pickup_xy.iloc[indices.ravel()]["xPos"].values
    iso_lab["pickupY"] = pickup_xy.iloc[indices.ravel()]["yPos"].values

    # Merge with pickup to get inject well
    iso_lab["xyft"] = (
        iso_lab["pickupX"].astype(str)
        + iso_lab["pickupY"].astype(str)
        + iso_lab["Field"].astype(str)
        + iso_lab["Target"].astype(float).astype(int).astype(str)
    )
    pickup["xyft"] = (
        pickup["xPos"].astype(str)
        + pickup["yPos"].astype(str)
        + pickup["field"].astype(str)
        + pickup["target"].astype(str)
    )
    iso_final = iso_lab.merge(
        pickup[["xyft", "well"]].rename(columns={"well": "injectWell"}),
        on="xyft",
        how="left",
    )

    # Build result dataframe
    result = pd.DataFrame({
        "sample": iso_final["condition"],
        "diameter": iso_final.get("Diameter", np.nan),
        "elongation": iso_final.get("Elongation", np.nan),
        "field": iso_final["Field"],
        "dropXPos": iso_final["XPos"],
        "dropYPos": iso_final["YPos"],
        "label": iso_final["well_x"] if "well_x" in iso_final.columns else iso_final.get("well"),
        "pickupXPos_numb": iso_final["pickupX"].astype(float),
        "pickupYPos_numb": iso_final["pickupY"].astype(float),
        "injectWell": iso_final["injectWell"],
    })

    # Map label to Reporter.intensity.N
    labs_map = {i: f"Reporter.intensity.{i + 3}" for i in range(1, plex + 1)}
    if plex == 14:
        labs_map = {str(i): f"Reporter.intensity.{i + 4}" for i in range(1, 15)}
    elif plex in (29, 32):
        labs_map = {str(i): f"Reporter.intensity.{i + 3}" for i in range(1, plex + 1)}
    result["label"] = result["label"].astype(str).map(labs_map).fillna(result["label"])
    result = result.dropna(subset=["label"])

    result["plate"] = 1
    result["ID"] = result["injectWell"].astype(str) + result["plate"].astype(str) + result["label"].astype(str)

    # Join green fluorescence stain data back if present
    if green_df is not None:
        result["match_stain"] = (
            "F" + result["field"].astype(str)
            + "X" + result["dropXPos"].astype(str)
            + "Y" + result["dropYPos"].astype(str)
        )
        green_df["match_stain"] = (
            "F" + green_df["Field"].astype(str)
            + "X" + green_df["XPos"].astype(str)
            + "Y" + green_df["YPos"].astype(str)
        )
        green_small = green_df[["match_stain", "Intensity", "Diameter"]].copy()
        green_small = green_small.rename(columns={"Diameter": "Stain_Diameter"})
        result = result.merge(green_small, on="match_stain", how="left")
        result = result.drop(columns=["match_stain"])

    return result


# ---------------------------------------------------------------------------
# mTRAQ CellenONE parser
# ---------------------------------------------------------------------------

def _analyze_cellenone_mtraq(cells_file: dict[str, str], plex: int) -> pd.DataFrame:
    """Parse CellenONE isolation files for mTRAQ (DIA) experiments."""
    frames = []
    for name, path in cells_file.items():
        df = pd.read_csv(path, sep="\t")
        df["condition"] = name
        frames.append(df)
    cells = pd.concat(frames, ignore_index=True)

    if "X" in cells.columns:
        cells = cells[~cells["X"].isin(["Green", "0"])]
        cells.loc[cells["X"].astype(str).str.contains("Transmission", na=False), "X"] = np.nan
        # Backward fill (upward) on all columns EXCEPT X, matching R's fill(2:7, .direction="up")
        fill_cols = [c for c in cells.columns if c != "X"]
        cells[fill_cols] = cells[fill_cols].bfill()
        cells = cells.dropna(subset=["XPos"])
        cells = cells[cells["YPos"] < 69]

    # Normalize Field to integer strings (cells may have float 1.0 → "1")
    cells["Field"] = cells["Field"].dropna().astype(int).astype(str)
    cells = cells[cells["Field"].notna()]

    plex_dir = {2: "2plex_files", 3: "3plex_files"}

    label = _parse_fld(_get_extdata(f"{plex_dir[plex]}/Labels.fld"))
    pickup1 = _parse_fld(_get_extdata(f"{plex_dir[plex]}/Pickup_mock_1.fld") if plex == 3
                         else _get_extdata(f"{plex_dir[plex]}/Pickup_1_mock.fld"))
    pickup2 = _parse_fld(_get_extdata(f"{plex_dir[plex]}/Pickup_mock_2.fld") if plex == 3
                         else _get_extdata(f"{plex_dir[plex]}/Pickup_2_mock.fld"),
                         skip=27)
    pickup = pd.concat([pickup1, pickup2], ignore_index=True)

    # Clean label wells
    label["well"] = label["well"].str[2:].str.rstrip(",")
    label["well"] = pd.to_numeric(label["well"], errors="coerce")
    unique_wells = sorted(label["well"].dropna().unique())
    well_to_tag = {w: str(i + 1) for i, w in enumerate(unique_wells)}
    label["well"] = label["well"].map(well_to_tag)

    # Clean pickup wells
    pickup["well"] = pickup["well"].str[1:].str.rstrip(",")
    pickup["target"] = 1

    # Filter to matching fields
    cell_fields = cells["Field"].unique()
    label = label[label["field"].isin(cell_fields)]
    pickup = pickup[pickup["field"].isin(cell_fields)]

    # Merge label with cells via xyf
    cells["xyf"] = cells["XPos"].astype(str) + cells["YPos"].astype(str) + cells["Field"].astype(str)
    label["xyf"] = label["xPos"].astype(str) + label["yPos"].astype(str) + label["field"].astype(str)
    iso_lab = cells.merge(label[["xyf", "well"]], on="xyf", how="left")

    # ANN for pickup
    pickup_xy = pickup[["xPos", "yPos"]].drop_duplicates().reset_index(drop=True)
    nn = NearestNeighbors(n_neighbors=1).fit(pickup_xy.values.astype(float))
    _, indices = nn.kneighbors(iso_lab[["XPos", "YPos"]].values.astype(float))
    iso_lab["pickupX"] = pickup_xy.iloc[indices.ravel()]["xPos"].values
    iso_lab["pickupY"] = pickup_xy.iloc[indices.ravel()]["yPos"].values

    iso_lab["xyft"] = (
        iso_lab["pickupX"].astype(str)
        + iso_lab["pickupY"].astype(str)
        + iso_lab["Field"].astype(str)
        + iso_lab["Target"].astype(float).astype(int).astype(str)
    )
    pickup["xyft"] = (
        pickup["xPos"].astype(str)
        + pickup["yPos"].astype(str)
        + pickup["field"].astype(str)
        + pickup["target"].astype(str)
    )
    iso_final = iso_lab.merge(
        pickup[["xyft", "well"]].rename(columns={"well": "injectWell"}),
        on="xyft",
        how="left",
    )

    result = pd.DataFrame({
        "sample": iso_final["condition"],
        "diameter": iso_final.get("Diameter", np.nan),
        "elongation": iso_final.get("Elongation", np.nan),
        "field": iso_final["Field"],
        "dropXPos": iso_final["XPos"],
        "dropYPos": iso_final["YPos"],
        "label": iso_final.get("well_x", iso_final.get("well")),
        "pickupXPos_numb": iso_final["pickupX"].astype(float),
        "pickupYPos_numb": iso_final["pickupY"].astype(float),
        "injectWell": iso_final["injectWell"],
    })

    # Plate assignment: fields 1-8 → plate 1, 9-16 → plate 2
    result["plate"] = np.where(result["field"].astype(float) > 8, 2, 1)

    # mTRAQ tag assignment
    result["label"] = pd.to_numeric(result["label"], errors="coerce")
    if plex == 2:
        tag_map = {1: "0", 3: "0", 5: "0", 7: "0", 2: "4", 4: "4", 6: "4", 8: "4"}
    else:  # plex == 3
        tag_map = {
            1: "0", 4: "0", 7: "0", 10: "0",
            2: "4", 5: "4", 8: "4", 11: "4",
            3: "8", 6: "8", 9: "8", 12: "8",
        }
    result["label"] = result["label"].map(tag_map)

    result["ID"] = (
        result["injectWell"].astype(str)
        + result["plate"].astype(str)
        + "."
        + result["label"].astype(str)
    )

    return result


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def link_cellenone_raw(qqc: QQC, cells_file: dict[str, str]) -> QQC:
    """Link CellenONE cell sorter metadata to the QQC object."""
    if qqc.ms_type == "DDA":
        co_data = _analyze_cellenone_tmt(cells_file, qqc.misc["plex"])
    else:
        co_data = _analyze_cellenone_mtraq(cells_file, qqc.misc["plex"])

    peptide_mat = qqc.matrices.peptide
    peptide_cols = qqc.matrices.peptide_cols
    cell_ids = pd.DataFrame({"ID": peptide_cols})

    base_cols = ["ID", "diameter", "sample", "label", "injectWell", "plate"]
    if "Intensity" in co_data.columns:
        base_cols += ["Intensity", "Stain_Diameter"]
    co_small = co_data[[c for c in base_cols if c in co_data.columns]].drop_duplicates(subset="ID")
    cell_ids = cell_ids.merge(co_small, on="ID", how="left")
    cell_ids["sample"] = cell_ids["sample"].fillna("neg")
    cell_ids["prot_total"] = np.log2(np.nansum(peptide_mat, axis=0))

    linker = qqc.meta_data.copy()
    cell_ids["WP"] = cell_ids["plate"].astype(str) + cell_ids["injectWell"].astype(str)
    linker["WP"] = linker["plate"].astype(str) + linker["Well"].astype(str)
    cell_ids = cell_ids.merge(linker, on="WP", how="left", suffixes=("", "_linker"))
    cell_ids = cell_ids.drop(columns=["WP"], errors="ignore")

    qqc.cellenone_meta = co_data
    qqc.meta_data = cell_ids
    return qqc


def link_manual_raw(qqc: QQC, cellenone_data: pd.DataFrame) -> QQC:
    """Link manually-annotated sample metadata (when CellenONE files unavailable)."""
    peptide_mat = qqc.matrices.peptide
    peptide_cols = qqc.matrices.peptide_cols
    cell_ids = pd.DataFrame({"ID": peptide_cols})

    co_small = cellenone_data[["ID", "diameter", "sample", "label", "injectWell", "plate"]].copy()
    cell_ids = cell_ids.merge(co_small, on="ID", how="left")
    cell_ids["sample"] = cell_ids["sample"].fillna("neg")
    cell_ids["prot_total"] = np.log2(np.nansum(peptide_mat, axis=0))

    qqc.cellenone_meta = cellenone_data
    qqc.meta_data = cell_ids
    return qqc


# ---------------------------------------------------------------------------
# Visualisation
# ---------------------------------------------------------------------------

def plot_slide_layout_celltype(qqc: QQC) -> plt.Figure:
    """Spatial scatter of droplet positions colored by sample type."""
    meta = qqc.cellenone_meta
    fields = sorted(meta["field"].unique())
    ncols = 4
    nrows = max(1, (len(fields) + ncols - 1) // ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4 * nrows))
    axes = np.atleast_2d(axes)
    for idx, fld in enumerate(fields):
        ax = axes[idx // ncols, idx % ncols]
        sub = meta[meta["field"] == fld]
        for samp in sub["sample"].unique():
            s = sub[sub["sample"] == samp]
            ax.scatter(s["dropXPos"], s["dropYPos"], label=samp, s=10, alpha=0.7)
        ax.invert_yaxis()
        ax.set_title(f"Field {fld}")
        dot_plot_style(ax)
    for idx in range(len(fields), nrows * ncols):
        axes[idx // ncols, idx % ncols].set_visible(False)
    fig.tight_layout()
    return fig


def plot_slide_layout_label(qqc: QQC) -> plt.Figure:
    """Spatial scatter of droplet positions colored by label."""
    meta = qqc.cellenone_meta
    fields = sorted(meta["field"].unique())
    ncols = 4
    nrows = max(1, (len(fields) + ncols - 1) // ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4 * nrows))
    axes = np.atleast_2d(axes)
    for idx, fld in enumerate(fields):
        ax = axes[idx // ncols, idx % ncols]
        sub = meta[meta["field"] == fld]
        for lab in sub["label"].unique():
            s = sub[sub["label"] == lab]
            ax.scatter(s["dropXPos"], s["dropYPos"], label=str(lab), s=10, alpha=0.7)
        ax.invert_yaxis()
        ax.set_title(f"Field {fld}")
    for idx in range(len(fields), nrows * ncols):
        axes[idx // ncols, idx % ncols].set_visible(False)
    fig.tight_layout()
    return fig
