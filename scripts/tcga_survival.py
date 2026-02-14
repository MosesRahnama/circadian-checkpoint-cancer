"""
Kaplan-Meier Survival Analysis by Circadian Coherence
======================================================
Stratifies patients by circadian CV quartiles and boundary-failure
classification, then computes KM survival curves and log-rank tests
using scipy.stats (CensoredData, ecdf, logrank).

Data: TCGA expanded TPM + clinical data across 6 cancer types.
"""
import os
import sys
import warnings

# ---- path setup & plotting backend BEFORE pyplot import ----
sys.path.insert(0, os.path.dirname(__file__))

import matplotlib
matplotlib.use("Agg")

import numpy as np
import pandas as pd
from scipy.stats import CensoredData, ecdf, logrank, false_discovery_control

from tcga_config import (RESULTS_CSV_DIR,
    CIRCADIAN, ALL_PROJECTS, PROJECT_LABELS,
    log_transform, compute_circadian_cv, classify_boundary_failure,
    DATA_DIR, FIGURE_DIR, setup_plotting, upload_to_gcs,
)

plt = setup_plotting()
warnings.filterwarnings("ignore", category=RuntimeWarning)

# ── Constants ─────────────────────────────────────────────────────────────
DAYS_PER_MONTH = 30.44
Q_COLORS = {1: "#2166ac", 2: "#4dac26", 3: "#f98e09", 4: "#b2182b"}
Q_LABELS = {1: "Q1 (most coherent)", 2: "Q2", 3: "Q3", 4: "Q4 (most decoherent)"}
MIN_PATIENTS = 50


# ── Helpers ───────────────────────────────────────────────────────────────
def deduplicate_one_tumor_sample_per_case(tumors):
    """
    Keep one tumor sample per (project_short_name, case_barcode).

    Priority:
      1) Primary Tumor
      2) Recurrent Tumor
      3) Metastatic / Additional Metastatic
      4) Any other tumor-like label

    This avoids violating patient-level independence in survival analysis.
    """
    pr = np.full(len(tumors), 4, dtype=int)
    s = tumors["sample_type_name"].fillna("").astype(str)
    pr[s.str.contains("Primary Tumor", case=False)] = 1
    pr[s.str.contains("Recurrent Tumor", case=False)] = 2
    pr[s.str.contains("Metastatic", case=False)] = 3

    out = tumors.copy()
    out["_sample_priority"] = pr
    out = out.sort_values(
        ["project_short_name", "case_barcode", "_sample_priority", "sample_barcode"]
    )
    out = out.drop_duplicates(
        subset=["project_short_name", "case_barcode"], keep="first"
    )
    out = out.drop(columns=["_sample_priority"])
    return out


def build_censored(times, events):
    """Return a CensoredData object from arrays of times and boolean events.

    CensoredData.right_censored(x, censored) takes:
      x: all observation times
      censored: boolean array where True = right-censored (alive/no event)
    """
    t = times.values.astype(float)
    is_censored = ~events.values  # True where NOT an event (alive = censored)
    return CensoredData.right_censored(t, is_censored)


def km_curve(cd):
    """Return (x_months, survival_prob) step arrays from a CensoredData object."""
    res = ecdf(cd)
    x = res.sf.quantiles
    y = res.sf.probabilities
    # prepend time=0, surv=1 for plotting
    x = np.concatenate([[0], x]) / DAYS_PER_MONTH
    y = np.concatenate([[1.0], y])
    return x, y


# ── Main ──────────────────────────────────────────────────────────────────
def main():
    print("=" * 65)
    print("  TCGA Survival Analysis -- Circadian Coherence Quartiles")
    print("=" * 65)

    # ---- load & merge ----
    expr = pd.read_csv(os.path.join(DATA_DIR, "tcga_expanded_tpm.csv"))
    clin = pd.read_csv(os.path.join(DATA_DIR, "tcga_clinical.csv"))
    print(f"Expression rows: {len(expr)}")
    print(f"Clinical  rows : {len(clin)}")

    # keep tumor samples only
    tumors = expr[
        expr["sample_type_name"].str.contains("Tumor|Metastatic", case=False, na=False)
    ].copy()
    print(f"Tumor rows (pre-dedup): {len(tumors)}")
    print(f"Unique cases (pre-dedup): {tumors['case_barcode'].nunique()}")
    tumors = deduplicate_one_tumor_sample_per_case(tumors)
    print(f"Tumor rows (deduped): {len(tumors)}")
    print(f"Unique cases (deduped): {tumors['case_barcode'].nunique()}")

    # merge (one clinical row per patient, many expression rows possible)
    merged = tumors.merge(clin, on="case_barcode", how="inner", suffixes=("", "_clin"))
    # resolve duplicated project column
    if "project_short_name_clin" in merged.columns:
        merged.drop(columns=["project_short_name_clin"], inplace=True)
    print(f"Merged rows    : {len(merged)}")

    # ---- survival time ----
    merged["days_to_death"] = pd.to_numeric(merged["days_to_death"], errors="coerce")
    merged["days_to_last_follow_up"] = pd.to_numeric(
        merged["days_to_last_follow_up"], errors="coerce"
    )

    merged["event"] = merged["vital_status"] == "Dead"
    merged["surv_time"] = np.where(
        merged["event"],
        merged["days_to_death"],
        merged["days_to_last_follow_up"],
    )
    merged["surv_time"] = pd.to_numeric(merged["surv_time"], errors="coerce")
    merged = merged.dropna(subset=["surv_time"])
    merged = merged[merged["surv_time"] > 0].copy()
    print(f"With valid survival: {len(merged)}")

    # ---- circadian CV & quartiles ----
    merged["circ_cv"] = compute_circadian_cv(merged)
    merged = merged.dropna(subset=["circ_cv"]).copy()

    merged["circ_q"] = merged.groupby("project_short_name")["circ_cv"].transform(
        lambda s: pd.qcut(s, 4, labels=[1, 2, 3, 4])
    )
    merged["circ_q"] = merged["circ_q"].astype(int)

    # ---- boundary failure classification ----
    merged["bf_class"] = classify_boundary_failure(merged)

    # ==================================================================
    #  FIGURE 1 -- KM curves by circadian-CV quartile, per cancer type
    # ==================================================================
    eligible = []
    for proj in ALL_PROJECTS:
        sub = merged[merged["project_short_name"] == proj]
        if len(sub) >= MIN_PATIENTS:
            eligible.append(proj)
    print(f"\nCancer types with >= {MIN_PATIENTS} patients: {eligible}")

    n_panels = min(len(eligible), 6)
    if n_panels == 0:
        print("No cancer types with enough patients. Exiting.")
        return

    nrows, ncols = 2, 3
    fig1, axes1 = plt.subplots(nrows, ncols, figsize=(17, 10))
    axes1 = axes1.flatten()

    surv_results = []

    for idx, proj in enumerate(eligible[:6]):
        ax = axes1[idx]
        sub = merged[merged["project_short_name"] == proj].copy()
        label = PROJECT_LABELS.get(proj, proj)

        for q in [1, 2, 3, 4]:
            qsub = sub[sub["circ_q"] == q]
            if len(qsub) < 5:
                continue
            cd = build_censored(qsub["surv_time"], qsub["event"])
            x, y = km_curve(cd)
            ax.step(x, y, where="post", color=Q_COLORS[q],
                    label=f"{Q_LABELS[q]} (n={len(qsub)})", linewidth=1.4)

        # log-rank Q1 vs Q4
        q1 = sub[sub["circ_q"] == 1]
        q4 = sub[sub["circ_q"] == 4]
        if len(q1) >= 5 and len(q4) >= 5:
            cd1 = build_censored(q1["surv_time"], q1["event"])
            cd4 = build_censored(q4["surv_time"], q4["event"])
            lr = logrank(cd1, cd4)
            p_txt = f"p = {lr.pvalue:.2e}" if lr.pvalue < 0.01 else f"p = {lr.pvalue:.3f}"
            ax.text(0.97, 0.03, f"Log-rank Q1 vs Q4\n{p_txt}",
                    transform=ax.transAxes, ha="right", va="bottom", fontsize=8,
                    bbox=dict(facecolor="white", alpha=0.8, edgecolor="gray", boxstyle="round"))
            surv_results.append({
                "cancer_type": proj,
                "test": "logrank_Q1_vs_Q4",
                "statistic": lr.statistic,
                "p_value": lr.pvalue,
                "n_Q1": len(q1),
                "n_Q4": len(q4),
            })
        else:
            surv_results.append({
                "cancer_type": proj,
                "test": "logrank_Q1_vs_Q4",
                "statistic": np.nan,
                "p_value": np.nan,
                "n_Q1": len(q1),
                "n_Q4": len(q4),
            })

        ax.set_title(label, fontweight="bold")
        ax.set_xlabel("Time (months)")
        ax.set_ylabel("Survival probability")
        ax.set_ylim(-0.02, 1.04)
        ax.legend(fontsize=7, loc="lower left")

    # blank unused panels
    for idx in range(len(eligible[:6]), nrows * ncols):
        axes1[idx].set_visible(False)

    fig1.suptitle("Kaplan-Meier Survival by Circadian Coherence Quartile",
                  fontsize=14, fontweight="bold", y=1.01)
    fig1.tight_layout()
    fig1_path = os.path.join(FIGURE_DIR, "survival_circadian_quartile.png")
    fig1.savefig(fig1_path, dpi=300, facecolor="white", bbox_inches="tight")
    plt.close(fig1)
    print(f"  Saved: {fig1_path}")

    # ==================================================================
    #  FIGURE 2 -- KM curves: Active Masking vs Decoherence
    # ==================================================================
    bf_colors = {"Active_Masking": "#7570b3", "Decoherence": "#d95f02"}

    fig2, axes2 = plt.subplots(nrows, ncols, figsize=(17, 10))
    axes2 = axes2.flatten()

    for idx, proj in enumerate(eligible[:6]):
        ax = axes2[idx]
        sub = merged[merged["project_short_name"] == proj].copy()
        label = PROJECT_LABELS.get(proj, proj)

        for bf_label, color in bf_colors.items():
            bfsub = sub[sub["bf_class"] == bf_label]
            if len(bfsub) < 5:
                continue
            cd = build_censored(bfsub["surv_time"], bfsub["event"])
            x, y = km_curve(cd)
            ax.step(x, y, where="post", color=color,
                    label=f"{bf_label.replace('_', ' ')} (n={len(bfsub)})",
                    linewidth=1.5)

        # log-rank AM vs DC
        am = sub[sub["bf_class"] == "Active_Masking"]
        dc = sub[sub["bf_class"] == "Decoherence"]
        if len(am) >= 5 and len(dc) >= 5:
            cd_am = build_censored(am["surv_time"], am["event"])
            cd_dc = build_censored(dc["surv_time"], dc["event"])
            lr = logrank(cd_am, cd_dc)
            p_txt = f"p = {lr.pvalue:.2e}" if lr.pvalue < 0.01 else f"p = {lr.pvalue:.3f}"
            ax.text(0.97, 0.03,
                    f"Log-rank AM vs DC\n{p_txt}",
                    transform=ax.transAxes, ha="right", va="bottom", fontsize=8,
                    bbox=dict(facecolor="white", alpha=0.8, edgecolor="gray", boxstyle="round"))
            surv_results.append({
                "cancer_type": proj,
                "test": "logrank_AM_vs_DC",
                "statistic": lr.statistic,
                "p_value": lr.pvalue,
                "n_Q1": len(am),
                "n_Q4": len(dc),
            })

        ax.set_title(label, fontweight="bold")
        ax.set_xlabel("Time (months)")
        ax.set_ylabel("Survival probability")
        ax.set_ylim(-0.02, 1.04)
        ax.legend(fontsize=8, loc="lower left")

    for idx in range(len(eligible[:6]), nrows * ncols):
        axes2[idx].set_visible(False)

    fig2.suptitle("Kaplan-Meier Survival: Active Masking vs Decoherence",
                  fontsize=14, fontweight="bold", y=1.01)
    fig2.tight_layout()
    fig2_path = os.path.join(FIGURE_DIR, "survival_boundary_failure.png")
    fig2.savefig(fig2_path, dpi=300, facecolor="white", bbox_inches="tight")
    plt.close(fig2)
    print(f"  Saved: {fig2_path}")

    # ---- save results CSV ----
    res_df = pd.DataFrame(surv_results)
    if not res_df.empty:
        # BH-FDR within each test family (Q1 vs Q4, AM vs DC) across 6 cohorts.
        res_df["q_value_test_family"] = np.nan
        for _, g in res_df.groupby("test", dropna=False):
            pvals = g["p_value"].astype(float)
            valid = pvals.notna().values
            if valid.sum() > 0:
                qvals = np.full(len(g), np.nan, dtype=float)
                qvals[valid] = false_discovery_control(pvals[valid], method="bh")
                res_df.loc[g.index, "q_value_test_family"] = qvals
        # BH-FDR across all 12 survival tests
        p_all = res_df["p_value"].astype(float)
        valid_all = p_all.notna().values
        q_all = np.full(len(res_df), np.nan, dtype=float)
        if valid_all.sum() > 0:
            q_all[valid_all] = false_discovery_control(p_all[valid_all], method="bh")
        res_df["q_value_all_survival_tests"] = q_all
        res_df["significant_fdr05_family"] = res_df["q_value_test_family"] < 0.05
        res_df["significant_fdr05_all"] = res_df["q_value_all_survival_tests"] < 0.05

    csv_path = os.path.join(RESULTS_CSV_DIR, "survival_logrank_results.csv")
    res_df.to_csv(csv_path, index=False)
    print(f"  Saved: {csv_path}")
    print(res_df.to_string(index=False))

    # ---- upload to GCS ----
    try:
        upload_to_gcs(csv_path, "analysis/survival_logrank_results.csv")
        upload_to_gcs(fig1_path, "figures/survival_circadian_quartile.png")
        upload_to_gcs(fig2_path, "figures/survival_boundary_failure.png")
    except Exception as e:
        print(f"  GCS upload skipped: {e}")

    print("\nSurvival analysis complete.")


if __name__ == "__main__":
    main()
