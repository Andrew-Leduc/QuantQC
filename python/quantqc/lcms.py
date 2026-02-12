"""LC/MS performance metrics — mirrors LC_MS_Performance.R."""

from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from quantqc.core import QQC
from quantqc.plotting import dot_plot_style


def calculate_run_order_statistics(qqc: QQC) -> QQC:
    """Calculate run-order statistics — dispatches DDA / DIA."""
    if qqc.ms_type == "DDA":
        return _calc_stats_dda(qqc)
    return _calc_stats_dia(qqc)


def _calc_stats_dda(qqc: QQC) -> QQC:
    raw = qqc.raw_data

    # PSM counts per run
    ids_per_run = raw.groupby("Order")["seqcharge"].nunique().reset_index()
    ids_per_run.columns = ["Run", "NumbRunIDs"]

    # MS1 intensity drift
    ms1 = raw.pivot_table(index="seqcharge", columns="Order", values="Intensity", aggfunc="first")
    first_col = ms1.columns[0]
    ms1_norm = np.log2(ms1.div(ms1[first_col], axis=0))
    ms1_means = ms1_norm.mean(skipna=True).reset_index()
    ms1_means.columns = ["Run", "MS1_means"]

    # MS2 intensity drift
    ms2 = raw.pivot_table(index="seqcharge", columns="Order", values="Reporter.intensity.1", aggfunc="first")
    ms2_norm = np.log2(ms2.div(ms2[ms2.columns[0]], axis=0))
    ms2_means = ms2_norm.mean(skipna=True).reset_index()
    ms2_means.columns = ["Run", "MS2_means"]

    # RT drift
    rt = raw.pivot_table(index="seqcharge", columns="Order", values="Retention.time", aggfunc="first")
    rt_adj = rt.sub(rt[rt.columns[0]], axis=0)
    rt_means = rt_adj.mean(skipna=True).reset_index()
    rt_means.columns = ["Run", "RT_means"]
    rt_sds = rt_adj.std(skipna=True).reset_index()
    rt_sds.columns = ["Run", "RT_sds"]
    rt_df = rt_means.merge(rt_sds, on="Run")

    qqc.run_order_statistics = [ids_per_run, ms1_means, ms2_means, rt_df]
    return qqc


def _calc_stats_dia(qqc: QQC) -> QQC:
    raw = qqc.raw_data

    ids_per_run = raw.groupby("Order").agg(NumbRunIDs=("Ms1.Area", lambda x: x.notna().sum())).reset_index()
    ids_per_run.columns = ["Run", "NumbRunIDs"]

    ms1_sum = raw.groupby("Order")["Ms1.Area"].sum().reset_index()
    ms1_sum.columns = ["Run", "MS1_means"]
    ms1_sum["MS1_means"] = ms1_sum["MS1_means"] / ms1_sum["MS1_means"].iloc[0]

    ms2_sum = raw.groupby("Order")["Precursor.Quantity"].sum().reset_index()
    ms2_sum.columns = ["Run", "MS2_means"]
    ms2_sum["MS2_means"] = ms2_sum["MS2_means"] / ms2_sum["MS2_means"].iloc[0]

    rt_stats = raw.groupby("Order")["RT"].agg(["mean", "std"]).reset_index()
    rt_stats.columns = ["Run", "RT_means", "RT_sds"]

    qqc.run_order_statistics = [ids_per_run, ms1_sum, ms2_sum, rt_stats]
    return qqc


def plot_intensity_drift(qqc: QQC) -> plt.Figure:
    """Three-panel plot: IDs, MS1, MS2 vs run order."""
    ids_df, ms1_df, ms2_df, _ = qqc.run_order_statistics

    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(8, 12))

    ax1.scatter(ids_df["Run"], ids_df["NumbRunIDs"], s=20)
    ax1.set_ylabel("# Precursor IDs")
    ax1.set_title("Run IDs as runs progress")
    dot_plot_style(ax1)

    ax2.scatter(ms1_df["Run"], ms1_df["MS1_means"], s=20)
    ax2.set_ylabel("Log2(Normalized to run 1)")
    ax2.set_title("MS1 Intensity as runs progress")
    dot_plot_style(ax2)

    ax3.scatter(ms2_df["Run"], ms2_df["MS2_means"], s=20)
    ax3.set_ylabel("Log2(Normalized to run 1)")
    ax3.set_title("MS2 Intensity as runs progress")
    ax3.set_xlabel("Run")
    dot_plot_style(ax3)

    fig.tight_layout()
    return fig


def plot_rt_drift(qqc: QQC) -> plt.Figure:
    """Two-panel plot: mean RT and RT SD vs run order."""
    rt_df = qqc.run_order_statistics[3]

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8))

    ax1.scatter(rt_df["Run"], rt_df["RT_means"], s=20)
    ax1.set_ylabel("Mean RT")
    ax1.set_title("Mean RT Drift as runs progress")
    dot_plot_style(ax1)

    ax2.scatter(rt_df["Run"], rt_df["RT_sds"], s=20)
    ax2.set_ylabel("RT SD")
    ax2.set_title("RT Standard Deviation as runs progress")
    ax2.set_xlabel("Run")
    dot_plot_style(ax2)

    fig.tight_layout()
    return fig
