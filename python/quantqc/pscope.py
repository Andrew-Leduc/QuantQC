"""pSCoPE target list generation — mirrors PriotitizedSCoPE.R."""

from __future__ import annotations

import numpy as np
import pandas as pd


def gen_tlist_diann(dia: pd.DataFrame, start_time: float) -> pd.DataFrame:
    """Generate prioritized target list from DIA-NN report for pSCoPE."""
    cols = [
        "Stripped.Sequence", "Modified.Sequence", "RT", "PEP",
        "Precursor.Charge", "Precursor.Mz", "Ms1.Area", "Genes",
        "Run", "Protein.Group",
    ]
    df = dia[cols].copy()

    df["Charge_double"] = df["Precursor.Charge"].astype(float)
    df["Mass"] = df["Precursor.Mz"] * df["Charge_double"] - 1.0072826748 * df["Charge_double"]
    df["MassRound"] = df["Mass"].round(7)
    df["SeqCharge"] = df["Stripped.Sequence"].astype(str) + df["Precursor.Charge"].astype(str)

    # Best PSM per SeqCharge (lowest PEP)
    df = df.sort_values("PEP").drop_duplicates(subset="SeqCharge", keep="first")

    df["RT_sec"] = (df["RT"] * 60).round(0).astype(int)
    df["TargBool"] = True
    df["Prot"] = df["Protein.Group"].str.split(";").str[0].str.split("-").str[0]
    df["RTBool"] = True
    df["Masses"] = "376.27"
    df["Leading.razor.protein"] = df["Protein.Group"].str.split(";").str[0]

    result = df[[
        "Stripped.Sequence", "RT", "Precursor.Charge", "Ms1.Area",
        "MassRound", "TargBool", "RTBool", "Masses",
        "Leading.razor.protein", "SeqCharge",
    ]].copy()

    result.columns = [
        "Sequence", "Retention.time", "Charge", "Apex.intensity",
        "Mass", "TargBoolean", "RTC_Boolean", "Fragments.mz",
        "Leading.razor.protein", "Modified",
    ]
    result["Retention.time"] = result["Retention.time"].round(3)
    result["Apex.intensity"] = result["Apex.intensity"].round(6)
    result["Retention.time"] = result["Retention.time"] - start_time

    # Priority tiers based on intensity tertiles
    q33, q66 = result["Apex.intensity"].quantile([0.33, 0.66])
    result["Priority"] = np.where(
        result["Apex.intensity"] <= q33, 1,
        np.where(result["Apex.intensity"] <= q66, 2, 3),
    )

    # Downgrade proteins with >4 high-priority peptides
    high = result[result["Priority"] == 3]
    prot_counts = high["Leading.razor.protein"].value_counts()
    excess_prots = prot_counts[prot_counts > 4].index

    downgrade = []
    for prot in excess_prots:
        sub = high[high["Leading.razor.protein"] == prot].sort_values("Apex.intensity", ascending=False)
        downgrade.extend(sub["Modified"].iloc[4:].tolist())

    result.loc[result["Modified"].isin(downgrade), "Priority"] = 1

    return result


def reformat_pscope(final_df: pd.DataFrame, scout: pd.DataFrame) -> pd.DataFrame:
    """Update priority tiers based on scout/SCP run results."""
    scout = scout.copy()
    scout["seqcharge"] = scout["Sequence"].astype(str) + scout["Charge"].astype(str)
    scout = scout[scout["PEP"] < 0.04]
    scout = scout[scout["PIF"] > 0.75]

    result = final_df.copy()
    result.loc[result["Priority"] == 3, "Priority"] = 2
    result.loc[result["Modified"].isin(scout["seqcharge"]), "Priority"] = 3

    # Re-downgrade excess
    high = result[result["Priority"] == 3]
    prot_counts = high["Leading.razor.protein"].value_counts()
    excess_prots = prot_counts[prot_counts > 4].index

    downgrade = []
    for prot in excess_prots:
        sub = high[high["Leading.razor.protein"] == prot].sort_values("Apex.intensity", ascending=False)
        downgrade.extend(sub["Modified"].iloc[4:].tolist())

    result.loc[result["Modified"].isin(downgrade), "Priority"] = 1

    return result
