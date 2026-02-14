"""
Immune-Fraction Residualization and Purity Stratification
=========================================================
Tests whether the circadian CV-PD-L1 coupling survives after removing
immune infiltration effects, addressing the stromal contamination concern.

Two complementary approaches:
  1. Continuous residualization: regress circadian CV on immune covariates,
     then correlate residual CV with PD-L1
  2. Purity-stratified: split into purity tertiles and recompute correlations
     within each stratum

Output:
  - experiments/tcga/immune_residualization_results.csv
  - results/immune_residualization_scatter.png
  - results/purity_stratified_correlations.png
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

from tcga_config import (RESULTS_CSV_DIR,
    ALL_PROJECTS, PROJECT_LABELS,
    log_transform, compute_circadian_cv,
    DATA_DIR,
)

warnings.filterwarnings("ignore", category=RuntimeWarning)

_PAPER_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS_DIR = os.path.join(_PAPER_ROOT, "results", "figures")
os.makedirs(RESULTS_DIR, exist_ok=True)

FOCUS_COHORTS = ["TCGA-SKCM", "TCGA-LUAD", "TCGA-BRCA", "TCGA-HNSC"]
IMMUNE_COVARIATES = [
    "leukocyte_fraction", "lymphocyte_infiltration",
    "ifn_gamma_response", "purity", "stromal_fraction",
]


def main():
    print("=" * 65)
    print("  Immune Residualization & Purity Stratification")
    print("=" * 65)

    expr = pd.read_csv(os.path.join(DATA_DIR, "tcga_expanded_tpm.csv"))
    immune = pd.read_csv(os.path.join(DATA_DIR, "tcga_purity_immune_covariates.csv"))

    tumors = expr[
        expr["sample_type_name"].str.contains("Tumor|Metastatic", case=False, na=False)
    ].copy()

    tumors["circ_cv"] = compute_circadian_cv(tumors)
    tumors["cd274_log"] = log_transform(tumors["CD274"])

    merged = tumors.merge(immune, on="case_barcode", how="inner", suffixes=("", "_imm"))
    if "project_short_name_imm" in merged.columns:
        merged.drop(columns=["project_short_name_imm"], inplace=True)
    print(f"Merged tumor + immune rows: {len(merged)}")

    # Ensure numeric
    for col in IMMUNE_COVARIATES:
        if col in merged.columns:
            merged[col] = pd.to_numeric(merged[col], errors="coerce")

    results = []

    # ================================================================
    #  PART 1: Continuous Residualization
    # ================================================================
    print("\n--- Part 1: Continuous Residualization ---")
    fig1, axes1 = plt.subplots(2, 2, figsize=(13, 11))
    axes1_flat = axes1.flatten()

    for idx, proj in enumerate(FOCUS_COHORTS):
        sub = merged[merged["project_short_name"] == proj].copy()
        label = PROJECT_LABELS.get(proj, proj)

        # Available covariates for this cohort
        avail_covs = [c for c in IMMUNE_COVARIATES if c in sub.columns and sub[c].notna().sum() > 10]
        valid = sub.dropna(subset=["circ_cv", "cd274_log"] + avail_covs).copy()
        n_valid = len(valid)

        if n_valid < 20:
            print(f"  {label}: insufficient data ({n_valid})")
            continue

        # Raw correlation
        rho_raw, p_raw = stats.spearmanr(valid["circ_cv"], valid["cd274_log"])

        # OLS residualization of circadian CV on immune covariates
        X = valid[avail_covs].values
        y_cv = valid["circ_cv"].values

        # Add intercept
        X_design = np.column_stack([np.ones(len(X)), X])
        try:
            beta, _, _, _ = np.linalg.lstsq(X_design, y_cv, rcond=None)
            y_pred = X_design @ beta
            residual_cv = y_cv - y_pred
            r2 = 1.0 - np.var(residual_cv) / np.var(y_cv)
        except np.linalg.LinAlgError:
            print(f"  {label}: OLS failed")
            continue

        # Correlation of residual CV with PD-L1
        rho_resid, p_resid = stats.spearmanr(residual_cv, valid["cd274_log"].values)

        results.append({
            "cohort": proj, "analysis": "residualization",
            "n": n_valid,
            "rho_raw": round(rho_raw, 4), "p_raw": p_raw,
            "immune_r2_on_cv": round(r2, 4),
            "rho_residual": round(rho_resid, 4), "p_residual": p_resid,
            "covariates_used": ", ".join(avail_covs),
        })

        sig_raw = "***" if p_raw < 0.001 else "**" if p_raw < 0.01 else "*" if p_raw < 0.05 else "ns"
        sig_res = "***" if p_resid < 0.001 else "**" if p_resid < 0.01 else "*" if p_resid < 0.05 else "ns"
        print(f"  {label} (n={n_valid}):  raw rho={rho_raw:.3f} {sig_raw}  |  "
              f"R2(immune->CV)={r2:.3f}  |  residual rho={rho_resid:.3f} {sig_res}")

        # Scatter plot
        ax = axes1_flat[idx]
        ax.scatter(valid["cd274_log"].values, residual_cv, s=8, alpha=0.35, color="#2c7fb8")
        # Trend line
        z = np.polyfit(valid["cd274_log"].values, residual_cv, 1)
        x_line = np.linspace(valid["cd274_log"].min(), valid["cd274_log"].max(), 100)
        ax.plot(x_line, np.polyval(z, x_line), color="red", linewidth=1.5, linestyle="--")
        ax.set_xlabel("log2(PD-L1 + 1)")
        ax.set_ylabel("Residual Circadian CV\n(immune effects removed)")
        ax.set_title(f"{label} (n={n_valid})\nresid rho={rho_resid:.3f}, p={p_resid:.2e}")
        ax.grid(alpha=0.2)
        ax.text(0.03, 0.97, f"Immune R2={r2:.3f}", transform=ax.transAxes,
                va="top", fontsize=9, bbox=dict(facecolor="white", alpha=0.8, edgecolor="none"))

    fig1.suptitle(
        "Circadian CV vs PD-L1 After Immune-Fraction Residualization",
        fontsize=13, fontweight="bold",
    )
    fig1.tight_layout()
    fig1_path = os.path.join(RESULTS_DIR, "immune_residualization_scatter.png")
    fig1.savefig(fig1_path, dpi=300, facecolor="white", bbox_inches="tight")
    plt.close(fig1)
    print(f"  Saved: {fig1_path}")

    # ================================================================
    #  PART 2: Purity-Stratified Correlations
    # ================================================================
    print("\n--- Part 2: Purity-Stratified Correlations ---")
    fig2, axes2 = plt.subplots(2, 2, figsize=(13, 10))
    axes2_flat = axes2.flatten()

    for idx, proj in enumerate(FOCUS_COHORTS):
        sub = merged[merged["project_short_name"] == proj].copy()
        label = PROJECT_LABELS.get(proj, proj)

        if "purity" not in sub.columns or sub["purity"].notna().sum() < 30:
            print(f"  {label}: insufficient purity data")
            axes2_flat[idx].text(0.5, 0.5, "Insufficient purity data",
                                 transform=axes2_flat[idx].transAxes, ha="center")
            continue

        valid = sub.dropna(subset=["circ_cv", "cd274_log", "purity"]).copy()

        # Tertile split
        valid["purity_tertile"] = pd.qcut(valid["purity"], 3, labels=["Low", "Mid", "High"])

        ax = axes2_flat[idx]
        colors_t = {"Low": "#d95f02", "Mid": "#7570b3", "High": "#1b9e77"}
        tertile_results = []

        for tert in ["Low", "Mid", "High"]:
            tsub = valid[valid["purity_tertile"] == tert]
            if len(tsub) < 10:
                continue
            rho_t, p_t = stats.spearmanr(tsub["circ_cv"], tsub["cd274_log"])
            tertile_results.append({
                "tertile": tert, "n": len(tsub), "rho": rho_t, "p": p_t,
            })
            results.append({
                "cohort": proj, "analysis": f"purity_tertile_{tert}",
                "n": len(tsub),
                "rho_raw": round(rho_t, 4), "p_raw": p_t,
                "immune_r2_on_cv": np.nan, "rho_residual": np.nan, "p_residual": np.nan,
                "covariates_used": f"purity_tertile={tert}",
            })

            ax.scatter(tsub["cd274_log"], tsub["circ_cv"], s=10, alpha=0.3,
                       color=colors_t[tert], label=f"{tert} (n={len(tsub)}, rho={rho_t:.2f})")

            sig_t = "***" if p_t < 0.001 else "**" if p_t < 0.01 else "*" if p_t < 0.05 else "ns"
            print(f"  {label} purity={tert:4s} (n={len(tsub)}): rho={rho_t:.3f} {sig_t}")

        ax.set_xlabel("log2(PD-L1 + 1)")
        ax.set_ylabel("Circadian CV")
        ax.set_title(f"{label}: Purity-Stratified", fontweight="bold")
        ax.legend(fontsize=8, loc="upper right")
        ax.grid(alpha=0.2)

    fig2.suptitle(
        "Circadian CV vs PD-L1 Within Tumor-Purity Tertiles",
        fontsize=13, fontweight="bold",
    )
    fig2.tight_layout()
    fig2_path = os.path.join(RESULTS_DIR, "purity_stratified_correlations.png")
    fig2.savefig(fig2_path, dpi=300, facecolor="white", bbox_inches="tight")
    plt.close(fig2)
    print(f"  Saved: {fig2_path}")

    # Save all results
    res_df = pd.DataFrame(results)
    csv_path = os.path.join(RESULTS_CSV_DIR, "immune_residualization_results.csv")
    res_df.to_csv(csv_path, index=False)
    print(f"\nSaved: {csv_path}")

    print("\nImmune residualization + purity stratification complete.")


if __name__ == "__main__":
    main()
