"""
Restricted Mean Survival Time (RMST) Analysis
==============================================
Computes RMST difference (AM minus DC) at multiple time horizons
as a PH-assumption-free alternative to Cox hazard ratios.

RMST = area under the Kaplan-Meier curve up to time tau.
No proportional hazards assumption required.

Output:
  - experiments/tcga/rmst_results.csv
  - results/rmst_am_vs_dc.png
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
from scipy.stats import CensoredData, ecdf

from tcga_config import (RESULTS_CSV_DIR,
    ALL_PROJECTS, PROJECT_LABELS,
    log_transform, compute_circadian_cv, classify_boundary_failure,
    DATA_DIR,
)

warnings.filterwarnings("ignore", category=RuntimeWarning)

_PAPER_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS_DIR = os.path.join(_PAPER_ROOT, "results", "figures")
os.makedirs(RESULTS_DIR, exist_ok=True)

TAU_MONTHS = [36, 60, 84, 120]
N_BOOTSTRAP = 1000
SEED = 42


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


def compute_rmst(times, events, tau_days):
    """
    Compute RMST (area under KM curve) up to tau_days.
    Uses trapezoidal integration on the survival function.
    """
    t = np.asarray(times, dtype=float)
    e = np.asarray(events, dtype=bool)

    # Filter to observations within tau
    if np.sum(t > 0) < 2:
        return np.nan

    cd = CensoredData.right_censored(t, ~e)
    res = ecdf(cd)
    sf_x = res.sf.quantiles
    sf_y = res.sf.probabilities

    # Prepend (0, 1.0) and cap at tau
    sf_x = np.concatenate([[0], sf_x])
    sf_y = np.concatenate([[1.0], sf_y])

    # Truncate at tau
    mask = sf_x <= tau_days
    if not mask.any():
        return tau_days  # everyone survived beyond tau

    x_trunc = sf_x[mask]
    y_trunc = sf_y[mask]

    # Add tau endpoint (step function: last value extends to tau)
    if x_trunc[-1] < tau_days:
        x_trunc = np.append(x_trunc, tau_days)
        y_trunc = np.append(y_trunc, y_trunc[-1])

    # Trapezoidal integration (step function → use left-rule)
    rmst = 0.0
    for i in range(len(x_trunc) - 1):
        rmst += y_trunc[i] * (x_trunc[i + 1] - x_trunc[i])

    return rmst


def bootstrap_rmst_diff(am_times, am_events, dc_times, dc_events, tau_days, n_boot, rng):
    """Bootstrap RMST difference (AM - DC) and return 95% CI."""
    diffs = []
    for _ in range(n_boot):
        idx_am = rng.integers(0, len(am_times), size=len(am_times))
        idx_dc = rng.integers(0, len(dc_times), size=len(dc_times))
        r_am = compute_rmst(am_times[idx_am], am_events[idx_am], tau_days)
        r_dc = compute_rmst(dc_times[idx_dc], dc_events[idx_dc], tau_days)
        if np.isfinite(r_am) and np.isfinite(r_dc):
            diffs.append(r_am - r_dc)
    if len(diffs) < 50:
        return np.nan, np.nan, np.nan
    diffs = np.array(diffs)
    return np.mean(diffs), np.percentile(diffs, 2.5), np.percentile(diffs, 97.5)


def main():
    print("=" * 65)
    print("  RMST Analysis: AM vs DC (PH-free survival comparison)")
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

    rng = np.random.default_rng(SEED)

    results = []
    focus = ["TCGA-SKCM", "TCGA-LUAD"]

    for proj in focus:
        sub = merged[merged["project_short_name"] == proj]
        am = sub[sub["bf_class"] == "Active_Masking"]
        dc = sub[sub["bf_class"] == "Decoherence"]
        label = PROJECT_LABELS.get(proj, proj)

        print(f"\n--- {label} ({proj}): AM={len(am)}, DC={len(dc)} ---")

        for tau_m in TAU_MONTHS:
            tau_d = tau_m * 30.44

            rmst_am = compute_rmst(am["surv_time"].values, am["event"].values, tau_d)
            rmst_dc = compute_rmst(dc["surv_time"].values, dc["event"].values, tau_d)
            diff = rmst_am - rmst_dc if np.isfinite(rmst_am) and np.isfinite(rmst_dc) else np.nan

            mean_diff, ci_lo, ci_hi = bootstrap_rmst_diff(
                am["surv_time"].values, am["event"].values,
                dc["surv_time"].values, dc["event"].values,
                tau_d, N_BOOTSTRAP, rng,
            )

            # Convert to months for interpretability
            rmst_am_m = rmst_am / 30.44 if np.isfinite(rmst_am) else np.nan
            rmst_dc_m = rmst_dc / 30.44 if np.isfinite(rmst_dc) else np.nan
            diff_m = diff / 30.44 if np.isfinite(diff) else np.nan
            ci_lo_m = ci_lo / 30.44 if np.isfinite(ci_lo) else np.nan
            ci_hi_m = ci_hi / 30.44 if np.isfinite(ci_hi) else np.nan

            sig = "SIG" if np.isfinite(ci_lo_m) and np.isfinite(ci_hi_m) and ci_lo_m > 0 else "ns"

            results.append({
                "cohort": proj, "tau_months": tau_m,
                "rmst_am_months": round(rmst_am_m, 2) if np.isfinite(rmst_am_m) else np.nan,
                "rmst_dc_months": round(rmst_dc_m, 2) if np.isfinite(rmst_dc_m) else np.nan,
                "diff_months": round(diff_m, 2) if np.isfinite(diff_m) else np.nan,
                "boot_ci95_lo": round(ci_lo_m, 2) if np.isfinite(ci_lo_m) else np.nan,
                "boot_ci95_hi": round(ci_hi_m, 2) if np.isfinite(ci_hi_m) else np.nan,
                "n_am": len(am), "n_dc": len(dc),
                "significant": sig,
            })

            print(f"  tau={tau_m:3d}mo  RMST_AM={rmst_am_m:.1f}  RMST_DC={rmst_dc_m:.1f}  "
                  f"diff={diff_m:+.1f}mo  95%CI=[{ci_lo_m:.1f}, {ci_hi_m:.1f}]  {sig}")

    res_df = pd.DataFrame(results)
    csv_path = os.path.join(RESULTS_CSV_DIR, "rmst_results.csv")
    res_df.to_csv(csv_path, index=False)
    print(f"\nSaved: {csv_path}")

    # --- Plot ---
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    for idx, proj in enumerate(focus):
        ax = axes[idx]
        sub = res_df[res_df["cohort"] == proj].copy()
        sub = sub.dropna(subset=["diff_months"])

        x = sub["tau_months"].values
        y = sub["diff_months"].values
        lo = sub["boot_ci95_lo"].values
        hi = sub["boot_ci95_hi"].values

        colors = ["#1b9e77" if l > 0 else "#7f7f7f" for l in lo]

        ax.bar(x, y, width=12, color=colors, alpha=0.85, edgecolor="black", linewidth=0.5)
        ax.errorbar(x, y, yerr=[y - lo, hi - y], fmt="none", ecolor="black",
                     capsize=5, linewidth=1.2)
        ax.axhline(0, color="black", linewidth=0.8)
        ax.set_xlabel("Time Horizon (months)")
        ax.set_ylabel("RMST Difference (AM - DC, months)")
        ax.set_title(PROJECT_LABELS.get(proj, proj), fontweight="bold")
        ax.set_xticks(x)
        ax.grid(axis="y", alpha=0.2)

        # Annotate
        for i in range(len(sub)):
            sig_txt = "*" if lo[i] > 0 else ""
            ax.text(x[i], hi[i] + 0.5, sig_txt, ha="center", fontsize=14, fontweight="bold")

    fig.suptitle(
        "RMST Difference: Active Masking vs Decoherence (Bootstrap 95% CI)",
        fontsize=13, fontweight="bold",
    )
    fig.tight_layout()
    fig_path = os.path.join(RESULTS_DIR, "rmst_am_vs_dc.png")
    fig.savefig(fig_path, dpi=300, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {fig_path}")

    print("\nRMST analysis complete.")


if __name__ == "__main__":
    main()
