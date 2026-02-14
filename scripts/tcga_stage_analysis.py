"""
Stage-Stratified Circadian Coherence Analysis + Comprehensive FDR
=================================================================
1. Analyzes circadian coherence (CV) across AJCC stages (I-IV).
2. Spearman trend + Kruskal-Wallis for CV, CD274, B2M, gap junction mean.
3. Collects ALL p-values from every analysis script and applies BH-FDR.
4. Saves a master results table.

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
from scipy import stats

from tcga_config import (RESULTS_CSV_DIR,
    CIRCADIAN, GAP_JUNCTION, ALL_PROJECTS, PROJECT_LABELS,
    log_transform, compute_circadian_cv, classify_boundary_failure,
    normalize_stage,
    DATA_DIR, FIGURE_DIR, setup_plotting, upload_to_gcs,
)

plt = setup_plotting()
warnings.filterwarnings("ignore", category=RuntimeWarning)

# ── Constants ─────────────────────────────────────────────────────────────
MIN_STAGED = 30
STAGE_ORDER = ["I", "II", "III", "IV"]
STAGE_NUM = {"I": 1, "II": 2, "III": 3, "IV": 4}


# ── Benjamini-Hochberg FDR ────────────────────────────────────────────────
def bh_fdr(pvals):
    """Return BH-adjusted p-values (q-values) from a 1-D array of raw p-values."""
    pvals = np.asarray(pvals, dtype=float)
    n = len(pvals)
    if n == 0:
        return np.array([])
    order = np.argsort(pvals)
    ranked = np.empty(n)
    ranked[order] = np.arange(1, n + 1)
    q = pvals * n / ranked
    # enforce monotonicity: walk backwards through sorted order
    q_sorted_idx = np.argsort(pvals)[::-1]
    cummin = np.inf
    for i in q_sorted_idx:
        cummin = min(cummin, q[i])
        q[i] = cummin
    return np.clip(q, 0, 1)


# ── Main ──────────────────────────────────────────────────────────────────
def main():
    print("=" * 65)
    print("  TCGA Stage-Stratified Circadian Analysis + FDR Correction")
    print("=" * 65)

    # ---- load & merge ----
    expr = pd.read_csv(os.path.join(DATA_DIR, "tcga_expanded_tpm.csv"))
    clin = pd.read_csv(os.path.join(DATA_DIR, "tcga_clinical.csv"))
    print(f"Expression rows: {len(expr)}")
    print(f"Clinical  rows : {len(clin)}")

    tumors = expr[
        expr["sample_type_name"].str.contains("Tumor|Metastatic", case=False, na=False)
    ].copy()

    merged = tumors.merge(clin, on="case_barcode", how="inner", suffixes=("", "_clin"))
    if "project_short_name_clin" in merged.columns:
        merged.drop(columns=["project_short_name_clin"], inplace=True)
    print(f"Merged tumor rows: {len(merged)}")

    # ---- normalize stage ----
    merged["stage_normalized"] = merged["ajcc_stage"].apply(normalize_stage)
    staged = merged.dropna(subset=["stage_normalized"]).copy()
    staged["stage_num"] = staged["stage_normalized"].map(STAGE_NUM)
    print(f"With valid stage : {len(staged)}")

    # ---- circadian CV ----
    staged["circ_cv"] = compute_circadian_cv(staged)

    # ---- gap junction mean (log2) ----
    gj_cols = [c for c in GAP_JUNCTION if c in staged.columns]
    staged["gj_mean"] = staged[gj_cols].apply(log_transform).mean(axis=1)

    # ---- CD274 and B2M log ----
    staged["cd274_log"] = log_transform(staged["CD274"])
    staged["b2m_log"] = log_transform(staged["B2M"])

    # ---- boundary failure ----
    staged["bf_class"] = classify_boundary_failure(staged)

    # ==================================================================
    #  Per-cancer-type stage analysis
    # ==================================================================
    stage_results = []
    eligible = []

    for proj in ALL_PROJECTS:
        sub = staged[staged["project_short_name"] == proj]
        if len(sub) >= MIN_STAGED:
            eligible.append(proj)
    print(f"\nCancer types with >= {MIN_STAGED} staged tumors: {eligible}")

    metrics = [
        ("circ_cv",    "Circadian CV"),
        ("cd274_log",  "log2(CD274+1)"),
        ("b2m_log",    "log2(B2M+1)"),
        ("gj_mean",    "Gap Junction Mean"),
    ]

    for proj in eligible:
        sub = staged[staged["project_short_name"] == proj].copy()
        label = PROJECT_LABELS.get(proj, proj)
        print(f"\n--- {label} ({proj}): n={len(sub)} ---")

        for metric_col, metric_name in metrics:
            valid = sub.dropna(subset=[metric_col, "stage_num"])
            if len(valid) < 10:
                continue

            # Spearman trend (ordinal stage vs metric)
            rho, p_sp = stats.spearmanr(valid["stage_num"], valid[metric_col])
            stage_results.append({
                "cancer_type": proj,
                "metric": metric_name,
                "test": "spearman_stage_trend",
                "statistic": rho,
                "p_value": p_sp,
                "n": len(valid),
            })

            # Kruskal-Wallis across stage groups
            groups = [
                valid.loc[valid["stage_normalized"] == s, metric_col].dropna().values
                for s in STAGE_ORDER
            ]
            groups = [g for g in groups if len(g) >= 3]
            if len(groups) >= 2:
                kw_stat, p_kw = stats.kruskal(*groups)
                stage_results.append({
                    "cancer_type": proj,
                    "metric": metric_name,
                    "test": "kruskal_wallis_stage",
                    "statistic": kw_stat,
                    "p_value": p_kw,
                    "n": len(valid),
                })

            sig = "***" if p_sp < 0.001 else "**" if p_sp < 0.01 else "*" if p_sp < 0.05 else "ns"
            print(f"  {metric_name:25s}  Spearman rho={rho:+.3f}  p={p_sp:.2e}  {sig}")

    # ==================================================================
    #  FIGURE -- box plots of circadian CV by stage
    # ==================================================================
    n_panels = min(len(eligible), 6)
    nrows, ncols = 2, 3
    fig, axes = plt.subplots(nrows, ncols, figsize=(17, 10))
    axes = axes.flatten()

    box_colors = {"I": "#4575b4", "II": "#74add1", "III": "#fdae61", "IV": "#d73027"}

    for idx, proj in enumerate(eligible[:6]):
        ax = axes[idx]
        sub = staged[staged["project_short_name"] == proj].copy()
        label = PROJECT_LABELS.get(proj, proj)

        data_by_stage = []
        stage_labels = []
        for s in STAGE_ORDER:
            vals = sub.loc[sub["stage_normalized"] == s, "circ_cv"].dropna()
            if len(vals) > 0:
                data_by_stage.append(vals.values)
                stage_labels.append(f"{s}\n(n={len(vals)})")

        if len(data_by_stage) < 2:
            ax.text(0.5, 0.5, "Insufficient data", transform=ax.transAxes,
                    ha="center", va="center")
            ax.set_title(label)
            continue

        bp = ax.boxplot(data_by_stage, tick_labels=stage_labels, patch_artist=True,
                        widths=0.55, showfliers=False)
        for i, (patch, s) in enumerate(zip(bp["boxes"], STAGE_ORDER[:len(data_by_stage)])):
            patch.set_facecolor(box_colors.get(s, "#cccccc"))
            patch.set_alpha(0.7)

        # overlay jitter
        for i, vals in enumerate(data_by_stage):
            jitter = np.random.normal(0, 0.06, size=len(vals))
            ax.scatter(np.full(len(vals), i + 1) + jitter, vals,
                       alpha=0.25, s=8, color="black", zorder=3)

        # annotate Spearman trend
        valid = sub.dropna(subset=["circ_cv", "stage_num"])
        if len(valid) >= 10:
            rho, p_sp = stats.spearmanr(valid["stage_num"], valid["circ_cv"])
            sig = "***" if p_sp < 0.001 else "**" if p_sp < 0.01 else "*" if p_sp < 0.05 else "ns"
            ax.text(0.97, 0.97,
                    f"Spearman rho = {rho:+.3f}\np = {p_sp:.2e} {sig}",
                    transform=ax.transAxes, ha="right", va="top", fontsize=8,
                    bbox=dict(facecolor="white", alpha=0.85, edgecolor="gray",
                              boxstyle="round"))

        ax.set_title(label, fontweight="bold")
        ax.set_xlabel("AJCC Stage")
        ax.set_ylabel("Circadian CV (lower = more coherent)")

    for idx in range(len(eligible[:6]), nrows * ncols):
        axes[idx].set_visible(False)

    fig.suptitle("Circadian Coherence (CV) by AJCC Stage",
                 fontsize=14, fontweight="bold", y=1.01)
    fig.tight_layout()
    fig_path = os.path.join(FIGURE_DIR, "stage_circadian_cv_boxplot.png")
    fig.savefig(fig_path, dpi=300, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    print(f"\n  Saved: {fig_path}")

    # ---- save stage results CSV ----
    stage_df = pd.DataFrame(stage_results)
    stage_csv = os.path.join(RESULTS_CSV_DIR, "stage_analysis_results.csv")
    stage_df.to_csv(stage_csv, index=False)
    print(f"  Saved: {stage_csv}")

    # ==================================================================
    #  COMPREHENSIVE FDR -- collect ALL p-values across analyses
    # ==================================================================
    print("\n" + "=" * 65)
    print("  Comprehensive BH-FDR correction across all analyses")
    print("=" * 65)

    all_pvals = []

    def append_entry(
        source, cancer_type, test_description, test_type, statistic, p_value, n
    ):
        all_pvals.append({
            "source": source,
            "cancer_type": cancer_type,
            "test_description": test_description,
            "test_type": test_type,
            "statistic": statistic,
            "p_value": p_value,
            "n": n,
        })

    # 1) Multi-cancer correlation results (preferred expanded file)
    mc_path = os.path.join(DATA_DIR, "tcga_multicancer_correlations.csv")
    loaded_multicancer = False
    if os.path.exists(mc_path):
        mc = pd.read_csv(mc_path)
        for _, row in mc.iterrows():
            append_entry(
                source="multi_cancer",
                cancer_type=row.get("project", ""),
                test_description=row.get("comparison", ""),
                test_type="spearman",
                statistic=row.get("rho", np.nan),
                p_value=row.get("p_value", np.nan),
                n=row.get("n", np.nan),
            )
        loaded_multicancer = True
        print(f"  Loaded {len(mc)} p-values from tcga_multicancer_correlations.csv")

    # 2) Legacy hypothesis-1 file only if expanded multi-cancer file is absent
    h1_path = os.path.join(DATA_DIR, "hypothesis1_correlations.csv")
    if (not loaded_multicancer) and os.path.exists(h1_path):
        h1 = pd.read_csv(h1_path)
        for _, row in h1.iterrows():
            append_entry(
                source="hypothesis1_correlations_legacy",
                cancer_type=row.get("project", ""),
                test_description=row.get("comparison", ""),
                test_type="spearman",
                statistic=row.get("rho", np.nan),
                p_value=row.get("p_value", np.nan),
                n=row.get("n", np.nan),
            )
        print(f"  Loaded {len(h1)} p-values from hypothesis1_correlations.csv")

    # 3) Tumor-normal paired/unpaired tests
    tn_path = os.path.join(DATA_DIR, "tumor_normal_comparison.csv")
    if os.path.exists(tn_path):
        tn = pd.read_csv(tn_path)
        for _, row in tn.iterrows():
            ctype = row.get("project", "")
            analyte = row.get("analyte", "")
            n_pairs = row.get("n_pairs", np.nan)
            n_t_unpaired = row.get("n_tumor_unpaired", np.nan)
            n_n_unpaired = row.get("n_normal_unpaired", np.nan)
            n_unpaired = (
                n_t_unpaired + n_n_unpaired
                if np.isfinite(n_t_unpaired) and np.isfinite(n_n_unpaired)
                else n_pairs
            )
            append_entry(
                source="tumor_normal",
                cancer_type=ctype,
                test_description=f"{analyte}_wilcoxon",
                test_type="wilcoxon",
                statistic=row.get("wilcoxon_stat", np.nan),
                p_value=row.get("wilcoxon_p", np.nan),
                n=n_pairs,
            )
            append_entry(
                source="tumor_normal",
                cancer_type=ctype,
                test_description=f"{analyte}_mannwhitney",
                test_type="mannwhitney",
                statistic=row.get("mannwhitney_stat", np.nan),
                p_value=row.get("mannwhitney_p", np.nan),
                n=n_unpaired,
            )
        print(f"  Loaded {2 * len(tn)} p-values from tumor_normal_comparison.csv")

    # 4) Immune subtype tests
    im_path = os.path.join(DATA_DIR, "immune_subtype_comparison.csv")
    if os.path.exists(im_path):
        im = pd.read_csv(im_path)
        for _, row in im.iterrows():
            ctype = row.get("project", "")
            append_entry(
                source="immune_subtype",
                cancer_type=ctype,
                test_description="circadian_cv_subtypes_kruskal",
                test_type="kruskal",
                statistic=row.get("kruskal_stat", np.nan),
                p_value=row.get("kruskal_p", np.nan),
                n=row.get("n_total", np.nan),
            )
            append_entry(
                source="immune_subtype",
                cancer_type=ctype,
                test_description="active_masking_vs_decoherence_mannwhitney",
                test_type="mannwhitney",
                statistic=row.get("mw_am_vs_dc_stat", np.nan),
                p_value=row.get("mw_am_vs_dc_p", np.nan),
                n=row.get("n_total", np.nan),
            )
        print(f"  Loaded {2 * len(im)} p-values from immune_subtype_comparison.csv")

    # 5) Survival log-rank results
    sv_path = os.path.join(DATA_DIR, "survival_logrank_results.csv")
    if os.path.exists(sv_path):
        sv = pd.read_csv(sv_path)
        for _, row in sv.iterrows():
            append_entry(
                source="survival",
                cancer_type=row.get("cancer_type", ""),
                test_description=row.get("test", ""),
                test_type="logrank",
                statistic=row.get("statistic", np.nan),
                p_value=row.get("p_value", np.nan),
                n=row.get("n_Q1", np.nan),
            )
        print(f"  Loaded {len(sv)} p-values from survival_logrank_results.csv")

    # 6) Stage analysis results (this script)
    for r in stage_results:
        all_pvals.append({
            "source": "stage_analysis",
            "cancer_type": r["cancer_type"],
            "test_description": f"{r['metric']}_{r['test']}",
            "test_type": r["test"],
            "statistic": r["statistic"],
            "p_value": r["p_value"],
            "n": r["n"],
        })

    master = pd.DataFrame(all_pvals)
    # drop rows without valid p-values
    master = master.dropna(subset=["p_value"]).copy()
    master = master[np.isfinite(master["p_value"])].copy()
    print(f"\n  Total p-values for FDR: {len(master)}")

    if len(master) > 0:
        master["q_value"] = bh_fdr(master["p_value"].values)
        master["significant_fdr05"] = master["q_value"] < 0.05

        n_sig = master["significant_fdr05"].sum()
        n_total = len(master)
        print(f"  Significant after BH-FDR (q < 0.05): {n_sig} / {n_total}")
        print(f"  ({n_sig / n_total * 100:.1f}%)")

        # sort by q-value
        master = master.sort_values("q_value").reset_index(drop=True)

        # print top results
        print("\n  Top 20 results by q-value:")
        top_cols = ["source", "cancer_type", "test_description", "statistic",
                    "p_value", "q_value", "significant_fdr05"]
        print(master[top_cols].head(20).to_string(index=False))

        # save master
        master_path = os.path.join(RESULTS_CSV_DIR, "master_fdr_results.csv")
        master.to_csv(master_path, index=False)
        print(f"\n  Saved: {master_path}")

        # summary by source
        print("\n  Summary by source:")
        summary = master.groupby("source").agg(
            total=("p_value", "count"),
            sig_fdr05=("significant_fdr05", "sum"),
            min_q=("q_value", "min"),
        ).reset_index()
        summary["pct_sig"] = (summary["sig_fdr05"] / summary["total"] * 100).round(1)
        print(summary.to_string(index=False))

    # ---- upload to GCS ----
    try:
        upload_to_gcs(stage_csv, "analysis/stage_analysis_results.csv")
        upload_to_gcs(fig_path, "figures/stage_circadian_cv_boxplot.png")
        if len(master) > 0:
            upload_to_gcs(master_path, "analysis/master_fdr_results.csv")
    except Exception as e:
        print(f"  GCS upload skipped: {e}")

    print("\nStage analysis + FDR correction complete.")


if __name__ == "__main__":
    main()
