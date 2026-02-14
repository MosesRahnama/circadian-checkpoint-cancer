"""
Threshold Sensitivity Sweep for AM/DC Boundary-Failure Classification
=====================================================================
Varies the percentile cutoff used to define Active Masking and Decoherence
subtypes, re-runs log-rank survival tests at each threshold, and produces
a stability plot showing whether the melanoma (SKCM) survival signal is
robust to threshold choice.

Output:
  - experiments/tcga/threshold_sensitivity_results.csv
  - results/threshold_sensitivity_sweep.png
"""
from __future__ import annotations

import os
import sys
import warnings

sys.path.insert(0, os.path.dirname(__file__))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import CensoredData, logrank

from tcga_config import (RESULTS_CSV_DIR,
    CIRCADIAN, ALL_PROJECTS, PROJECT_LABELS,
    log_transform, compute_circadian_cv,
    DATA_DIR,
)

warnings.filterwarnings("ignore", category=RuntimeWarning)

_PAPER_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS_DIR = os.path.join(_PAPER_ROOT, "results", "figures")
os.makedirs(RESULTS_DIR, exist_ok=True)

DAYS_PER_MONTH = 30.44
FOCUS_COHORTS = ["TCGA-SKCM", "TCGA-LUAD"]
PERCENTILE_GRID = list(range(25, 76, 5))  # 25, 30, 35, ..., 75


def classify_at_percentile(df, pct, project_col="project_short_name"):
    """Classify AM/DC using a given percentile threshold instead of median."""
    labels = pd.Series("Mixed", index=df.index)
    for proj in df[project_col].unique():
        mask = df[project_col] == proj
        sub = df.loc[mask]

        cd274_thresh = log_transform(sub["CD274"]).quantile(pct / 100.0)
        arntl_thresh = log_transform(sub["ARNTL"]).quantile(pct / 100.0)
        per1_thresh = log_transform(sub["PER1"]).quantile(pct / 100.0)
        b2m_thresh = log_transform(sub["B2M"]).quantile(pct / 100.0)

        cd274_log = log_transform(sub["CD274"])
        arntl_log = log_transform(sub["ARNTL"])
        per1_log = log_transform(sub["PER1"])
        b2m_log = log_transform(sub["B2M"])

        am = (cd274_log > cd274_thresh) & (arntl_log > arntl_thresh) & (per1_log < per1_thresh)
        dc = (cd274_log < cd274_thresh) & (b2m_log < b2m_thresh)

        labels.loc[mask & am] = "Active_Masking"
        labels.loc[mask & dc] = "Decoherence"
    return labels


def deduplicate_one_tumor_per_case(tumors):
    """Keep one tumor sample per case (Primary > Recurrent > Metastatic)."""
    pr = np.full(len(tumors), 4, dtype=int)
    s = tumors["sample_type_name"].fillna("").astype(str)
    pr[s.str.contains("Primary Tumor", case=False)] = 1
    pr[s.str.contains("Recurrent Tumor", case=False)] = 2
    pr[s.str.contains("Metastatic", case=False)] = 3
    out = tumors.copy()
    out["_priority"] = pr
    out = out.sort_values(["project_short_name", "case_barcode", "_priority", "sample_barcode"])
    out = out.drop_duplicates(subset=["project_short_name", "case_barcode"], keep="first")
    return out.drop(columns=["_priority"])


def main():
    print("=" * 65)
    print("  Threshold Sensitivity Sweep: AM/DC Survival Stability")
    print("=" * 65)

    expr = pd.read_csv(os.path.join(DATA_DIR, "tcga_expanded_tpm.csv"))
    clin = pd.read_csv(os.path.join(DATA_DIR, "tcga_clinical.csv"))

    tumors = expr[
        expr["sample_type_name"].str.contains("Tumor|Metastatic", case=False, na=False)
    ].copy()
    tumors = deduplicate_one_tumor_per_case(tumors)

    merged = tumors.merge(clin, on="case_barcode", how="inner", suffixes=("", "_clin"))
    if "project_short_name_clin" in merged.columns:
        merged.drop(columns=["project_short_name_clin"], inplace=True)

    merged["days_to_death"] = pd.to_numeric(merged["days_to_death"], errors="coerce")
    merged["days_to_last_follow_up"] = pd.to_numeric(merged["days_to_last_follow_up"], errors="coerce")
    merged["event"] = merged["vital_status"] == "Dead"
    merged["surv_time"] = np.where(
        merged["event"], merged["days_to_death"], merged["days_to_last_follow_up"]
    )
    merged["surv_time"] = pd.to_numeric(merged["surv_time"], errors="coerce")
    merged = merged.dropna(subset=["surv_time"])
    merged = merged[merged["surv_time"] > 0].copy()
    print(f"Patients with valid survival: {len(merged)}")

    results = []
    for pct in PERCENTILE_GRID:
        merged["bf_class"] = classify_at_percentile(merged, pct)
        for proj in FOCUS_COHORTS:
            sub = merged[merged["project_short_name"] == proj]
            am = sub[sub["bf_class"] == "Active_Masking"]
            dc = sub[sub["bf_class"] == "Decoherence"]

            n_am, n_dc = len(am), len(dc)
            if n_am < 5 or n_dc < 5:
                results.append({
                    "percentile": pct, "cohort": proj, "n_am": n_am, "n_dc": n_dc,
                    "logrank_stat": np.nan, "p_value": np.nan,
                })
                continue

            cd_am = CensoredData.right_censored(
                am["surv_time"].values.astype(float),
                ~am["event"].values,
            )
            cd_dc = CensoredData.right_censored(
                dc["surv_time"].values.astype(float),
                ~dc["event"].values,
            )
            lr = logrank(cd_am, cd_dc)
            results.append({
                "percentile": pct, "cohort": proj, "n_am": n_am, "n_dc": n_dc,
                "logrank_stat": lr.statistic, "p_value": lr.pvalue,
            })
            direction = "AM better" if lr.statistic < 0 else "DC better"
            sig = "***" if lr.pvalue < 0.001 else "**" if lr.pvalue < 0.01 else "*" if lr.pvalue < 0.05 else "ns"
            print(f"  pct={pct:2d}  {PROJECT_LABELS.get(proj, proj):12s}  "
                  f"AM={n_am:3d}  DC={n_dc:3d}  p={lr.pvalue:.4f} {sig}  ({direction})")

    res_df = pd.DataFrame(results)
    csv_path = os.path.join(RESULTS_CSV_DIR, "threshold_sensitivity_results.csv")
    res_df.to_csv(csv_path, index=False)
    print(f"\nSaved: {csv_path}")

    # --- Plot ---
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5), sharey=True)
    for idx, proj in enumerate(FOCUS_COHORTS):
        ax = axes[idx]
        sub = res_df[res_df["cohort"] == proj].copy()
        sub = sub.dropna(subset=["p_value"])

        neglog10p = -np.log10(sub["p_value"].clip(lower=1e-20))
        ax.bar(sub["percentile"], neglog10p, width=4, color="#2c7fb8", alpha=0.85)
        ax.axhline(-np.log10(0.05), color="red", linestyle="--", linewidth=1.2, label="p = 0.05")
        ax.axhline(-np.log10(0.01), color="orange", linestyle=":", linewidth=1.0, label="p = 0.01")
        ax.axvline(50, color="black", linestyle="-", linewidth=0.8, alpha=0.4, label="Median (default)")

        ax.set_xlabel("Classification Percentile Threshold")
        ax.set_title(PROJECT_LABELS.get(proj, proj), fontweight="bold")
        ax.set_xticks(PERCENTILE_GRID)
        ax.set_xticklabels([str(p) for p in PERCENTILE_GRID], fontsize=8)
        ax.legend(fontsize=8, loc="upper right")
        ax.grid(axis="y", alpha=0.2)

        # Annotate n_AM / n_DC at each bar
        for _, row in sub.iterrows():
            y_pos = -np.log10(max(row["p_value"], 1e-20))
            ax.text(
                row["percentile"], y_pos + 0.05,
                f"{int(row['n_am'])}/{int(row['n_dc'])}",
                ha="center", va="bottom", fontsize=6.5, rotation=90,
            )

    axes[0].set_ylabel(r"$-\log_{10}(p)$  [log-rank AM vs DC]")
    fig.suptitle(
        "Threshold Sensitivity: AM vs DC Survival Signal Stability",
        fontsize=13, fontweight="bold",
    )
    fig.tight_layout()
    fig_path = os.path.join(RESULTS_DIR, "threshold_sensitivity_sweep.png")
    fig.savefig(fig_path, dpi=300, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {fig_path}")

    # Summary
    skcm = res_df[res_df["cohort"] == "TCGA-SKCM"]
    sig_range = skcm[skcm["p_value"] < 0.05]
    if not sig_range.empty:
        lo, hi = sig_range["percentile"].min(), sig_range["percentile"].max()
        print(f"\nSKCM: AM-vs-DC p < 0.05 from percentile {lo} to {hi} "
              f"({len(sig_range)}/{len(skcm)} thresholds)")
    else:
        print("\nSKCM: No threshold yields p < 0.05")

    print("\nThreshold sensitivity sweep complete.")


if __name__ == "__main__":
    main()
