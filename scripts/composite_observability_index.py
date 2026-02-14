"""
Composite Observability Index
=============================
Constructs an independent "boundary observability" score using features
NOT directly reused from the AM/DC classification rule, then tests whether
this index predicts survival and correlates with boundary-failure mode.

This breaks the tautology concern: if a metric defined independently of
the subtype rule still separates AM from DC and predicts survival,
the framework has genuine predictive content.

Index components (orthogonal to AM/DC definition where possible):
  - MHC-I composite: mean log2(HLA_A, HLA_B, HLA_C)  [NOT B2M]
  - Gap junction composite: mean log2(GJA1, GJB2, GJA5, GJB6)
  - Circadian CV (orthogonal): CV computed from CLOCK, PER2, CRY1, CRY2
    (excluding ARNTL and PER1, which are used in the AM definition)

Higher index = more observable/distinguishable to the organism.

Output:
  - experiments/tcga/observability_index_results.csv
  - results/observability_index_survival.png
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
from scipy import stats
from scipy.stats import CensoredData, logrank

from tcga_config import (RESULTS_CSV_DIR,
    ALL_PROJECTS, PROJECT_LABELS,
    log_transform, classify_boundary_failure,
    DATA_DIR,
)

warnings.filterwarnings("ignore", category=RuntimeWarning)

_PAPER_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS_DIR = os.path.join(_PAPER_ROOT, "results", "figures")
os.makedirs(RESULTS_DIR, exist_ok=True)

# Orthogonal clock genes (excluding ARNTL and PER1 used in AM definition)
ORTHO_CLOCK = ["CLOCK", "PER2", "CRY1", "CRY2"]
# MHC-I components (excluding B2M used in DC definition)
MHC_I_ORTHO = ["HLA_A", "HLA_B", "HLA_C"]
# Gap junctions (not used in AM/DC definition at all)
GAP_JUNCTION = ["GJA1", "GJB2", "GJA5", "GJB6"]

FOCUS_COHORTS = ["TCGA-SKCM", "TCGA-LUAD"]
DAYS_PER_MONTH = 30.44


def compute_ortho_circadian_cv(df):
    """CV across orthogonal clock genes only."""
    cols = [c for c in ORTHO_CLOCK if c in df.columns]
    vals = df[cols].apply(log_transform)
    return vals.std(axis=1) / vals.mean(axis=1)


def compute_observability_index(df):
    """
    Composite observability index: higher = more observable.
    z-scored within cohort, then averaged.
    """
    # Component 1: MHC-I (higher = more antigen presentation)
    mhc_cols = [c for c in MHC_I_ORTHO if c in df.columns]
    mhc_mean = df[mhc_cols].apply(log_transform).mean(axis=1)

    # Component 2: Gap junction (higher = more connected)
    gj_cols = [c for c in GAP_JUNCTION if c in df.columns]
    gj_mean = df[gj_cols].apply(log_transform).mean(axis=1)

    # Component 3: Circadian coherence (lower CV = more coherent = more observable)
    # Invert so higher = more coherent
    ortho_cv = compute_ortho_circadian_cv(df)
    inv_cv = -ortho_cv  # higher = more coherent

    # Z-score each component within-cohort
    index = pd.Series(np.nan, index=df.index)
    for proj in df["project_short_name"].unique():
        mask = df["project_short_name"] == proj
        sub_mhc = mhc_mean[mask]
        sub_gj = gj_mean[mask]
        sub_cv = inv_cv[mask]

        z_mhc = (sub_mhc - sub_mhc.mean()) / sub_mhc.std()
        z_gj = (sub_gj - sub_gj.mean()) / sub_gj.std()
        z_cv = (sub_cv - sub_cv.mean()) / sub_cv.std()

        # Equal-weight average
        index[mask] = (z_mhc + z_gj + z_cv) / 3.0

    return index


def deduplicate_one_tumor_per_case(tumors):
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
    print("  Composite Observability Index: Independent of AM/DC Rule")
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

    merged["bf_class"] = classify_boundary_failure(merged)
    merged["obs_index"] = compute_observability_index(merged)
    merged = merged.dropna(subset=["obs_index"])

    results = []

    # ---- Test 1: Does the index separate AM from DC? ----
    print("\n--- Test 1: Observability Index by Boundary-Failure Subtype ---")
    for proj in FOCUS_COHORTS:
        sub = merged[merged["project_short_name"] == proj]
        label = PROJECT_LABELS.get(proj, proj)
        am = sub[sub["bf_class"] == "Active_Masking"]["obs_index"]
        dc = sub[sub["bf_class"] == "Decoherence"]["obs_index"]

        if len(am) < 5 or len(dc) < 5:
            continue

        stat, p = stats.mannwhitneyu(am, dc, alternative="two-sided")
        sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
        print(f"  {label}: AM mean={am.mean():.3f}  DC mean={dc.mean():.3f}  "
              f"MW p={p:.2e} {sig}")
        results.append({
            "cohort": proj, "test": "am_vs_dc_index",
            "am_mean": round(am.mean(), 4), "dc_mean": round(dc.mean(), 4),
            "mw_stat": stat, "p_value": p,
        })

    # ---- Test 2: Does the index predict survival (median split)? ----
    print("\n--- Test 2: Observability Index and Survival ---")
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))

    for idx, proj in enumerate(FOCUS_COHORTS):
        ax = axes[idx]
        sub = merged[merged["project_short_name"] == proj].copy()
        label = PROJECT_LABELS.get(proj, proj)

        med = sub["obs_index"].median()
        high = sub[sub["obs_index"] >= med]
        low = sub[sub["obs_index"] < med]

        if len(high) < 10 or len(low) < 10:
            continue

        cd_hi = CensoredData.right_censored(
            high["surv_time"].values.astype(float), ~high["event"].values)
        cd_lo = CensoredData.right_censored(
            low["surv_time"].values.astype(float), ~low["event"].values)

        from scipy.stats import ecdf as scipy_ecdf
        res_hi = scipy_ecdf(cd_hi)
        res_lo = scipy_ecdf(cd_lo)

        # Plot KM curves
        for res, color, lbl, n in [
            (res_hi, "#1b9e77", f"High Obs (n={len(high)})", len(high)),
            (res_lo, "#d95f02", f"Low Obs (n={len(low)})", len(low)),
        ]:
            x = np.concatenate([[0], res.sf.quantiles]) / DAYS_PER_MONTH
            y = np.concatenate([[1.0], res.sf.probabilities])
            ax.step(x, y, where="post", color=color, label=lbl, linewidth=1.5)

        lr = logrank(cd_hi, cd_lo)
        p_txt = f"p = {lr.pvalue:.2e}" if lr.pvalue < 0.01 else f"p = {lr.pvalue:.3f}"
        ax.text(0.97, 0.03, f"Log-rank\n{p_txt}",
                transform=ax.transAxes, ha="right", va="bottom", fontsize=9,
                bbox=dict(facecolor="white", alpha=0.8, edgecolor="gray", boxstyle="round"))

        ax.set_xlabel("Time (months)")
        ax.set_ylabel("Survival Probability")
        ax.set_title(f"{label}", fontweight="bold")
        ax.set_ylim(-0.02, 1.04)
        ax.legend(fontsize=9, loc="lower left")
        ax.grid(alpha=0.2)

        sig = "***" if lr.pvalue < 0.001 else "**" if lr.pvalue < 0.01 else "*" if lr.pvalue < 0.05 else "ns"
        print(f"  {label}: High-obs={len(high)}, Low-obs={len(low)}, "
              f"logrank p={lr.pvalue:.4f} {sig}")
        results.append({
            "cohort": proj, "test": "obs_index_survival",
            "am_mean": np.nan, "dc_mean": np.nan,
            "mw_stat": lr.statistic, "p_value": lr.pvalue,
        })

    fig.suptitle(
        "Survival by Composite Observability Index\n"
        "(MHC-I + Gap Junction + Orthogonal Clock; independent of AM/DC definition)",
        fontsize=12, fontweight="bold",
    )
    fig.tight_layout()
    fig_path = os.path.join(RESULTS_DIR, "observability_index_survival.png")
    fig.savefig(fig_path, dpi=300, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {fig_path}")

    # ---- Test 3: Spearman correlation between index and PD-L1 ----
    print("\n--- Test 3: Observability Index vs PD-L1 ---")
    for proj in FOCUS_COHORTS:
        sub = merged[merged["project_short_name"] == proj]
        label = PROJECT_LABELS.get(proj, proj)
        cd274_log = log_transform(sub["CD274"])
        rho, p = stats.spearmanr(sub["obs_index"], cd274_log)
        sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
        print(f"  {label}: rho={rho:.3f}  p={p:.2e}  {sig}")
        results.append({
            "cohort": proj, "test": "obs_index_vs_pdl1",
            "am_mean": np.nan, "dc_mean": np.nan,
            "mw_stat": rho, "p_value": p,
        })

    res_df = pd.DataFrame(results)
    csv_path = os.path.join(RESULTS_CSV_DIR, "observability_index_results.csv")
    res_df.to_csv(csv_path, index=False)
    print(f"\nSaved: {csv_path}")

    print("\nComposite observability index analysis complete.")


if __name__ == "__main__":
    main()
