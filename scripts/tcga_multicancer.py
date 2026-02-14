"""
Multi-Cancer Correlation Analysis (6 TCGA Cohorts)
====================================================
Extends Hypothesis 1 across SKCM, LUAD, BRCA, COAD, HNSC, LUSC.

For each cancer type, computes Spearman correlations between:
  - CD274 (PD-L1) vs gap junction genes
  - CD274 vs circadian genes
  - CD274 vs circadian CV (coherence proxy)
  - CD274 vs MHC-I genes
  - Gap junction mean vs MHC-I mean
  - CDH1 vs VIM (EMT axis)

Applies Benjamini-Hochberg FDR correction across all tests.

Outputs:
  Figure 1 -- Correlation heatmap (rows=comparisons, cols=cancer types)
  Figure 2 -- Forest plot of circadian CV vs CD274 rho with 95% CI
  CSV     -- tcga_multicancer_correlations.csv

Data: tcga_expanded_tpm.csv (21 gene columns, 6 cancer types)
"""
import sys
import os

sys.path.insert(0, os.path.dirname(__file__))

import numpy as np
import pandas as pd
from scipy import stats

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from tcga_config import (RESULTS_CSV_DIR,
    CHECKPOINT, MHC_I, GAP_JUNCTION, CIRCADIAN, DIFFERENTIATION,
    ALL_PROJECTS, PROJECT_LABELS,
    log_transform, compute_circadian_cv, safe_alias,
    DATA_DIR, FIGURE_DIR, setup_plotting, upload_to_gcs,
)

# ── Plotting defaults ─────────────────────────────────────────────────────
plt = setup_plotting()


# ── Load data ──────────────────────────────────────────────────────────────
DATA_PATH = os.path.join(DATA_DIR, "tcga_expanded_tpm.csv")
df_all = pd.read_csv(DATA_PATH)
print(f"Loaded {len(df_all)} samples from {DATA_PATH}")
print(f"Projects present: {df_all['project_short_name'].value_counts().to_dict()}")


# ── Correlation engine ─────────────────────────────────────────────────────
def run_correlations(df_subset, project_name):
    """
    Run Spearman correlations for one cancer type.
    Returns a DataFrame with columns: project, comparison, category, rho, p_value, n.
    """
    results = []
    pdl1 = log_transform(df_subset["CD274"])

    # CD274 vs each gap junction gene (4 tests)
    for gj in GAP_JUNCTION:
        if gj not in df_subset.columns:
            continue
        gj_expr = log_transform(df_subset[gj])
        mask = pdl1.notna() & gj_expr.notna()
        r, p = stats.spearmanr(pdl1[mask], gj_expr[mask])
        results.append({
            "project": project_name,
            "comparison": f"CD274 vs {gj}",
            "category": "PD-L1 vs Gap Junction",
            "rho": r, "p_value": p, "n": int(mask.sum()),
        })

    # CD274 vs each circadian gene (6 tests)
    for circ in CIRCADIAN:
        if circ not in df_subset.columns:
            continue
        circ_expr = log_transform(df_subset[circ])
        mask = pdl1.notna() & circ_expr.notna()
        r, p = stats.spearmanr(pdl1[mask], circ_expr[mask])
        results.append({
            "project": project_name,
            "comparison": f"CD274 vs {circ}",
            "category": "PD-L1 vs Circadian",
            "rho": r, "p_value": p, "n": int(mask.sum()),
        })

    # CD274 vs circadian CV (1 test)
    circ_cv = compute_circadian_cv(df_subset)
    mask = pdl1.notna() & circ_cv.notna() & np.isfinite(circ_cv)
    r, p = stats.spearmanr(pdl1[mask], circ_cv[mask])
    results.append({
        "project": project_name,
        "comparison": "CD274 vs Circadian_CV",
        "category": "PD-L1 vs Circadian Coherence",
        "rho": r, "p_value": p, "n": int(mask.sum()),
    })

    # CD274 vs each MHC-I gene (4 tests)
    for mhc in MHC_I:
        if mhc not in df_subset.columns:
            continue
        mhc_expr = log_transform(df_subset[mhc])
        mask = pdl1.notna() & mhc_expr.notna()
        r, p = stats.spearmanr(pdl1[mask], mhc_expr[mask])
        results.append({
            "project": project_name,
            "comparison": f"CD274 vs {mhc}",
            "category": "PD-L1 vs MHC-I",
            "rho": r, "p_value": p, "n": int(mask.sum()),
        })

    # Gap junction mean vs MHC-I mean (1 test)
    gj_cols = [c for c in GAP_JUNCTION if c in df_subset.columns]
    mhc_cols = [c for c in MHC_I if c in df_subset.columns]
    gj_mean = df_subset[gj_cols].apply(log_transform).mean(axis=1)
    mhc_mean = df_subset[mhc_cols].apply(log_transform).mean(axis=1)
    mask = gj_mean.notna() & mhc_mean.notna()
    r, p = stats.spearmanr(gj_mean[mask], mhc_mean[mask])
    results.append({
        "project": project_name,
        "comparison": "GapJunction_mean vs MHC_I_mean",
        "category": "Gap Junction vs MHC-I",
        "rho": r, "p_value": p, "n": int(mask.sum()),
    })

    # CDH1 vs VIM -- EMT axis (1 test)
    if "CDH1" in df_subset.columns and "VIM" in df_subset.columns:
        cdh1 = log_transform(df_subset["CDH1"])
        vim = log_transform(df_subset["VIM"])
        mask = cdh1.notna() & vim.notna()
        r, p = stats.spearmanr(cdh1[mask], vim[mask])
        results.append({
            "project": project_name,
            "comparison": "CDH1 vs VIM",
            "category": "EMT axis",
            "rho": r, "p_value": p, "n": int(mask.sum()),
        })

    return pd.DataFrame(results)


# ── Run across all 6 cancer types ─────────────────────────────────────────
def run_all_cancers():
    """Run correlations for every project in ALL_PROJECTS."""
    frames = []
    for proj in ALL_PROJECTS:
        subset = df_all[df_all["project_short_name"] == proj].copy()
        tumor = subset[
            subset["sample_type_name"].str.contains(
                "Tumor|Metastatic", case=False, na=False
            )
        ]
        n_tumor = len(tumor)
        if n_tumor == 0:
            print(f"  WARNING: {proj} has 0 tumor samples -- skipping")
            continue

        print(f"\n{'=' * 60}")
        print(f"  {proj} ({PROJECT_LABELS.get(proj, proj)}): {n_tumor} tumor samples")
        print(f"{'=' * 60}")

        res = run_correlations(tumor, proj)
        frames.append(res)

        for _, row in res.iterrows():
            sig = (
                "***" if row["p_value"] < 0.001
                else "**" if row["p_value"] < 0.01
                else "*" if row["p_value"] < 0.05
                else "ns"
            )
            print(
                f"  {row['comparison']:35s}  rho={row['rho']:+.3f}"
                f"  p={row['p_value']:.2e}  n={row['n']:4d}  {sig}"
            )

    results = pd.concat(frames, ignore_index=True)
    return results


# ── FDR correction ─────────────────────────────────────────────────────────
def apply_fdr(results_df):
    """
    Apply Benjamini-Hochberg FDR correction across ALL p-values
    from every cancer type. Adds q_value and significant columns.
    """
    p_values = results_df["p_value"].values
    q_values = stats.false_discovery_control(p_values, method="bh")
    results_df = results_df.copy()
    results_df["q_value"] = q_values
    results_df["significant"] = q_values < 0.05
    n_sig = results_df["significant"].sum()
    print(f"\nFDR correction: {n_sig}/{len(results_df)} tests significant (q < 0.05)")
    return results_df


# ── Figure 1: Correlation heatmap ──────────────────────────────────────────
def plot_heatmap(results_df):
    """
    Heatmap with rows = comparisons, cols = cancer types.
    Color = Spearman rho (diverging RdBu_r), asterisks for q < 0.05.
    """
    comparisons = results_df["comparison"].unique()
    projects = ALL_PROJECTS

    n_rows = len(comparisons)
    n_cols = len(projects)

    rho_matrix = np.full((n_rows, n_cols), np.nan)
    sig_matrix = np.full((n_rows, n_cols), False)

    for i, comp in enumerate(comparisons):
        for j, proj in enumerate(projects):
            match = results_df[
                (results_df["comparison"] == comp) &
                (results_df["project"] == proj)
            ]
            if len(match) == 1:
                rho_matrix[i, j] = match["rho"].values[0]
                sig_matrix[i, j] = match["significant"].values[0]

    # Symmetric color limits
    vmax = np.nanmax(np.abs(rho_matrix))
    vmax = max(vmax, 0.1)  # floor so the colorbar is not degenerate

    fig_height = max(6, 0.45 * n_rows + 2)
    fig, ax = plt.subplots(figsize=(10, fig_height))

    im = ax.imshow(
        rho_matrix,
        aspect="auto",
        cmap="RdBu_r",
        vmin=-vmax,
        vmax=vmax,
        interpolation="nearest",
    )

    # Asterisks for FDR-significant cells
    for i in range(n_rows):
        for j in range(n_cols):
            if sig_matrix[i, j]:
                ax.text(
                    j, i, "*",
                    ha="center", va="center",
                    fontsize=14, fontweight="bold", color="black",
                )

    # Labels
    col_labels = [PROJECT_LABELS.get(p, p) for p in projects]
    ax.set_xticks(range(n_cols))
    ax.set_xticklabels(col_labels, rotation=45, ha="right", fontsize=9)
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(comparisons, fontsize=8)

    ax.set_title(
        "Spearman rho across 6 TCGA cohorts (* = FDR q < 0.05)",
        fontsize=12, fontweight="bold", pad=12,
    )

    cbar = fig.colorbar(im, ax=ax, shrink=0.6, pad=0.02)
    cbar.set_label("Spearman rho", fontsize=10)

    plt.tight_layout()
    out_path = os.path.join(FIGURE_DIR, "multicancer_correlation_heatmap.png")
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved heatmap: {out_path}")
    return out_path


# ── Figure 2: Forest plot -- Circadian CV vs CD274 ─────────────────────────
def fisher_z_ci(rho, n, alpha=0.05):
    """
    95% CI for Spearman rho via Fisher z-transform.
    Returns (rho_lower, rho_upper).
    """
    if n < 4:
        return (np.nan, np.nan)
    z = np.arctanh(rho)
    se = 1.0 / np.sqrt(n - 3)
    z_crit = stats.norm.ppf(1 - alpha / 2)
    z_lo = z - z_crit * se
    z_hi = z + z_crit * se
    return (np.tanh(z_lo), np.tanh(z_hi))


def plot_forest(results_df):
    """
    Forest plot of Circadian CV vs CD274 rho across 6 cancer types.
    Horizontal error bars show 95% CI via Fisher z-transform.
    """
    cv_rows = results_df[results_df["comparison"] == "CD274 vs Circadian_CV"].copy()
    cv_rows = cv_rows.sort_values("rho")

    labels = []
    rhos = []
    ci_lo = []
    ci_hi = []

    for _, row in cv_rows.iterrows():
        proj = row["project"]
        lbl = PROJECT_LABELS.get(proj, proj)
        n = int(row["n"])
        labels.append(f"{lbl} (n={n})")
        rhos.append(row["rho"])
        lo, hi = fisher_z_ci(row["rho"], n)
        ci_lo.append(lo)
        ci_hi.append(hi)

    rhos = np.array(rhos)
    ci_lo = np.array(ci_lo)
    ci_hi = np.array(ci_hi)
    y_pos = np.arange(len(labels))

    fig, ax = plt.subplots(figsize=(8, 5))

    # Horizontal error bars: xerr expects (2, N) array with lower/upper distances
    xerr_lo = rhos - ci_lo
    xerr_hi = ci_hi - rhos
    xerr = np.vstack([xerr_lo, xerr_hi])

    ax.errorbar(
        rhos, y_pos,
        xerr=xerr,
        fmt="o",
        color="steelblue",
        ecolor="gray",
        elinewidth=1.5,
        capsize=4,
        markersize=7,
        markeredgecolor="navy",
        markeredgewidth=0.8,
    )

    ax.axvline(0, color="black", linewidth=0.8, linestyle="--", alpha=0.6)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=9)
    ax.set_xlabel("Spearman rho (CD274 vs Circadian CV)", fontsize=10)
    ax.set_title(
        "Circadian Coherence vs PD-L1: Forest Plot (95% CI)",
        fontsize=12, fontweight="bold",
    )

    # Annotate rho values next to points
    for i, (r, q) in enumerate(
        zip(rhos, cv_rows["q_value"].values)
    ):
        sig_marker = " *" if q < 0.05 else ""
        ax.text(
            ci_hi[i] + 0.01, i,
            f"rho={r:+.3f}{sig_marker}",
            va="center", fontsize=8,
        )

    ax.invert_yaxis()
    plt.tight_layout()
    out_path = os.path.join(FIGURE_DIR, "multicancer_circadian_cv_forest.png")
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved forest plot: {out_path}")
    return out_path


# ── Main ───────────────────────────────────────────────────────────────────
def main():
    print("=" * 60)
    print("  Multi-Cancer Correlation Analysis (6 TCGA Cohorts)")
    print("=" * 60)

    # 1. Run correlations across all 6 cancer types
    results = run_all_cancers()

    # 2. Apply FDR correction (Benjamini-Hochberg) across all p-values
    results = apply_fdr(results)

    # 3. Save CSV
    csv_path = os.path.join(RESULTS_CSV_DIR, "tcga_multicancer_correlations.csv")
    results.to_csv(csv_path, index=False)
    print(f"\nCorrelations saved: {csv_path}")

    # 4. Generate figures
    print("\nGenerating figures...")
    os.makedirs(FIGURE_DIR, exist_ok=True)
    heatmap_path = plot_heatmap(results)
    forest_path = plot_forest(results)

    # 5. Upload to GCS
    print("\nUploading to GCS...")
    try:
        upload_to_gcs(csv_path, "analysis/tcga_multicancer_correlations.csv")
        upload_to_gcs(heatmap_path, "figures/multicancer_correlation_heatmap.png")
        upload_to_gcs(forest_path, "figures/multicancer_circadian_cv_forest.png")
    except Exception as e:
        print(f"  GCS upload failed: {e}")
        print("  (Files saved locally; upload manually if needed)")

    # 6. Summary table
    print("\n" + "=" * 60)
    print("  Summary: significant correlations by category")
    print("=" * 60)
    summary = (
        results.groupby("category")
        .agg(
            total=("significant", "size"),
            n_significant=("significant", "sum"),
            mean_rho=("rho", "mean"),
        )
        .reset_index()
    )
    summary["pct_significant"] = (
        100.0 * summary["n_significant"] / summary["total"]
    ).round(1)
    print(summary.to_string(index=False))

    print("\nDone.")


if __name__ == "__main__":
    main()
