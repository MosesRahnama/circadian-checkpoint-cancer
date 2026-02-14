"""
Immune Subtype Classification and Circadian Coherence
======================================================
Classifies tumors into boundary-failure subtypes:
  - Active Masking: high PD-L1 + high BMAL1 + low PER1
  - Decoherence:    low PD-L1 + low B2M
  - Mixed:          everything else

Then compares circadian coherence (CV) across subtypes.

Statistical tests:
  - Kruskal-Wallis across 3 groups
  - Pairwise Mann-Whitney U (Active_Masking vs Decoherence)
  - Benjamini-Hochberg FDR correction
"""
import sys
import os

sys.path.insert(0, os.path.dirname(__file__))

import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import false_discovery_control

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from tcga_config import (RESULTS_CSV_DIR,
    CHECKPOINT, MHC_I, GAP_JUNCTION, CIRCADIAN,
    ALL_PROJECTS, PROJECT_LABELS,
    log_transform, compute_circadian_cv, classify_boundary_failure,
    DATA_DIR, FIGURE_DIR,
    setup_plotting, upload_to_gcs,
)

# ── Constants ─────────────────────────────────────────────────────────────
CSV_PATH = os.path.join(DATA_DIR, r"tcga_expanded_tpm.csv")

SUBTYPE_COLORS = {
    "Active_Masking": "#D62728",   # red
    "Decoherence": "#1F77B4",      # blue
    "Mixed": "#999999",            # gray
}
SUBTYPE_ORDER = ["Active_Masking", "Decoherence", "Mixed"]


def load_tumor_data():
    """Load expression data and filter to tumor samples only."""
    df = pd.read_csv(CSV_PATH)
    print(f"Loaded {len(df)} total samples")

    tumor_mask = df["sample_type_name"].str.contains(
        "Tumor|Metastatic", case=False, na=False
    )
    df_tumor = df[tumor_mask].copy()
    print(f"Tumor samples: {len(df_tumor)}")
    return df_tumor


def sig_label(p):
    """Return significance stars."""
    if p < 0.001:
        return "***"
    elif p < 0.01:
        return "**"
    elif p < 0.05:
        return "*"
    else:
        return "ns"


def run_subtype_tests(df_tumor):
    """
    For each cancer type:
      - Kruskal-Wallis across 3 subtypes for circadian CV
      - Pairwise Mann-Whitney U: Active_Masking vs Decoherence
      - Compute group means for circadian CV and BMAL1/PER1 ratio
    """
    results = []

    for project in ALL_PROJECTS:
        proj_df = df_tumor[df_tumor["project_short_name"] == project].copy()
        if len(proj_df) < 20:
            print(f"  Skipping {project}: only {len(proj_df)} tumor samples")
            continue

        label = PROJECT_LABELS.get(project, project)

        # Compute circadian CV
        proj_df["circadian_cv"] = compute_circadian_cv(proj_df)

        # Get subtype labels
        proj_df["subtype"] = classify_boundary_failure(proj_df).values

        # Group counts
        counts = proj_df["subtype"].value_counts()
        print(f"\n  {project} ({label}): {len(proj_df)} tumors")
        for st in SUBTYPE_ORDER:
            print(f"    {st}: {counts.get(st, 0)}")

        # Extract groups
        am_cv = proj_df.loc[
            proj_df["subtype"] == "Active_Masking", "circadian_cv"
        ].dropna()
        dc_cv = proj_df.loc[
            proj_df["subtype"] == "Decoherence", "circadian_cv"
        ].dropna()
        mx_cv = proj_df.loc[
            proj_df["subtype"] == "Mixed", "circadian_cv"
        ].dropna()

        # Kruskal-Wallis (need >= 2 groups with >= 5 samples)
        groups_for_kw = [g for g in [am_cv, dc_cv, mx_cv] if len(g) >= 5]
        if len(groups_for_kw) >= 2:
            kw_stat, kw_p = stats.kruskal(*groups_for_kw)
        else:
            kw_stat, kw_p = np.nan, np.nan

        # Pairwise Mann-Whitney: Active_Masking vs Decoherence
        if len(am_cv) >= 5 and len(dc_cv) >= 5:
            mw_stat, mw_p = stats.mannwhitneyu(
                am_cv, dc_cv, alternative="two-sided"
            )
        else:
            mw_stat, mw_p = np.nan, np.nan

        # Compute mean BMAL1/PER1 ratio per group (log scale)
        bmal1_per1_ratios = {}
        for st in SUBTYPE_ORDER:
            sub = proj_df[proj_df["subtype"] == st]
            if len(sub) > 0:
                arntl_log = log_transform(sub["ARNTL"])
                per1_log = log_transform(sub["PER1"])
                # log ratio = log(ARNTL) - log(PER1) in log2 space
                ratio = (arntl_log - per1_log).mean()
                bmal1_per1_ratios[st] = float(ratio)
            else:
                bmal1_per1_ratios[st] = np.nan

        results.append({
            "project": project,
            "project_label": label,
            "n_total": len(proj_df),
            "n_active_masking": int(counts.get("Active_Masking", 0)),
            "n_decoherence": int(counts.get("Decoherence", 0)),
            "n_mixed": int(counts.get("Mixed", 0)),
            "cv_mean_active_masking": float(am_cv.mean()) if len(am_cv) > 0 else np.nan,
            "cv_mean_decoherence": float(dc_cv.mean()) if len(dc_cv) > 0 else np.nan,
            "cv_mean_mixed": float(mx_cv.mean()) if len(mx_cv) > 0 else np.nan,
            "bmal1_per1_ratio_am": bmal1_per1_ratios.get("Active_Masking", np.nan),
            "bmal1_per1_ratio_dc": bmal1_per1_ratios.get("Decoherence", np.nan),
            "bmal1_per1_ratio_mx": bmal1_per1_ratios.get("Mixed", np.nan),
            "kruskal_stat": float(kw_stat) if np.isfinite(kw_stat) else np.nan,
            "kruskal_p": float(kw_p) if np.isfinite(kw_p) else np.nan,
            "mw_am_vs_dc_stat": float(mw_stat) if np.isfinite(mw_stat) else np.nan,
            "mw_am_vs_dc_p": float(mw_p) if np.isfinite(mw_p) else np.nan,
        })

    return pd.DataFrame(results)


def apply_fdr(results_df):
    """Apply BH-FDR to Kruskal-Wallis and Mann-Whitney p-values."""
    for col in ["kruskal_p", "mw_am_vs_dc_p"]:
        pvals = results_df[col].values
        fdr_col = col.replace("_p", "_fdr")
        finite_mask = np.isfinite(pvals)
        adjusted = np.full_like(pvals, np.nan)
        if finite_mask.sum() > 0:
            adjusted[finite_mask] = false_discovery_control(
                pvals[finite_mask], method="bh"
            )
        results_df[fdr_col] = adjusted
    return results_df


def plot_circadian_cv_by_subtype(df_tumor, plt_handle):
    """
    Figure 1: 2x3 grid of box + strip plots.
    Circadian CV by boundary-failure subtype, per cancer type.
    """
    fig, axes = plt_handle.subplots(2, 3, figsize=(16, 10))
    fig.suptitle(
        "Circadian CV by Boundary-Failure Subtype",
        fontsize=14, fontweight="bold",
    )
    axes_flat = axes.flatten()

    for idx, project in enumerate(ALL_PROJECTS):
        ax = axes_flat[idx]
        proj_df = df_tumor[df_tumor["project_short_name"] == project].copy()
        label = PROJECT_LABELS.get(project, project)

        if len(proj_df) < 20:
            ax.set_visible(False)
            continue

        proj_df["circadian_cv"] = compute_circadian_cv(proj_df)
        proj_df["subtype"] = classify_boundary_failure(proj_df).values

        # Prepare data for boxplot
        box_data = []
        box_labels = []
        box_colors = []
        for st in SUBTYPE_ORDER:
            vals = proj_df.loc[
                proj_df["subtype"] == st, "circadian_cv"
            ].dropna().values
            box_data.append(vals)
            n_st = len(vals)
            box_labels.append(f"{st.replace('_', ' ')}\n(n={n_st})")
            box_colors.append(SUBTYPE_COLORS[st])

        # Box plots
        bp = ax.boxplot(
            box_data,
            positions=range(1, len(SUBTYPE_ORDER) + 1),
            widths=0.6,
            patch_artist=True,
            showfliers=False,
        )
        for patch, color in zip(bp["boxes"], box_colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.5)

        # Strip (jitter) overlay
        for i, (vals, color) in enumerate(zip(box_data, box_colors)):
            if len(vals) == 0:
                continue
            jitter = np.random.default_rng(42).normal(0, 0.08, size=len(vals))
            ax.scatter(
                np.full(len(vals), i + 1) + jitter,
                vals,
                color=color, alpha=0.3, s=8, edgecolors="none",
            )

        ax.set_xticks(range(1, len(SUBTYPE_ORDER) + 1))
        ax.set_xticklabels(box_labels, fontsize=8)
        ax.set_ylabel("Circadian CV")

        # Kruskal-Wallis annotation
        groups_for_kw = [d for d in box_data if len(d) >= 5]
        if len(groups_for_kw) >= 2:
            _, kw_p = stats.kruskal(*groups_for_kw)
            ax.set_title(
                f"{label}\nKruskal-Wallis p={kw_p:.2e} {sig_label(kw_p)}",
                fontsize=10,
            )
        else:
            ax.set_title(f"{label}", fontsize=10)

    plt_handle.tight_layout()
    out_path = os.path.join(FIGURE_DIR, "immune_subtype_circadian_cv.png")
    fig.savefig(out_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt_handle.close(fig)
    print(f"  Saved: {out_path}")
    return out_path


def plot_boundary_scatter(df_tumor, plt_handle):
    """
    Figure 2: 2D scatter of log2(CD274+1) vs log2(B2M+1),
    colored by subtype, for the largest cohort (BRCA).
    Quadrant lines at within-project medians.
    """
    # Use BRCA as the largest cohort
    brca = df_tumor[df_tumor["project_short_name"] == "TCGA-BRCA"].copy()
    if len(brca) < 20:
        print("  BRCA cohort too small for scatter plot -- skipping")
        return None

    brca["subtype"] = classify_boundary_failure(brca).values
    brca["cd274_log"] = log_transform(brca["CD274"])
    brca["b2m_log"] = log_transform(brca["B2M"])

    cd274_med = brca["cd274_log"].median()
    b2m_med = brca["b2m_log"].median()

    fig, ax = plt_handle.subplots(figsize=(9, 8))

    for st in SUBTYPE_ORDER:
        sub = brca[brca["subtype"] == st]
        ax.scatter(
            sub["cd274_log"], sub["b2m_log"],
            c=SUBTYPE_COLORS[st],
            label=f"{st.replace('_', ' ')} (n={len(sub)})",
            alpha=0.45, s=18, edgecolors="none",
        )

    # Quadrant lines at medians
    ax.axvline(cd274_med, color="black", linestyle="--", linewidth=0.8, alpha=0.5)
    ax.axhline(b2m_med, color="black", linestyle="--", linewidth=0.8, alpha=0.5)

    ax.set_xlabel("log2(CD274 + 1)  [PD-L1]", fontsize=11)
    ax.set_ylabel("log2(B2M + 1)  [MHC-I]", fontsize=11)
    ax.set_title(
        "Boundary-Failure Subtypes: PD-L1 vs MHC-I (BRCA)",
        fontsize=13, fontweight="bold",
    )
    ax.legend(loc="upper left", frameon=True, fontsize=9)

    # Annotate quadrants
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    offset_x = (xlim[1] - xlim[0]) * 0.02
    offset_y = (ylim[1] - ylim[0]) * 0.02
    ax.text(
        cd274_med + offset_x, b2m_med - offset_y,
        "Active Masking\nregion",
        fontsize=8, alpha=0.6, color=SUBTYPE_COLORS["Active_Masking"],
        va="top",
    )
    ax.text(
        cd274_med - offset_x, b2m_med - offset_y,
        "Decoherence\nregion",
        fontsize=8, alpha=0.6, color=SUBTYPE_COLORS["Decoherence"],
        ha="right", va="top",
    )

    plt_handle.tight_layout()
    out_path = os.path.join(FIGURE_DIR, "immune_subtype_boundary_scatter.png")
    fig.savefig(out_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt_handle.close(fig)
    print(f"  Saved: {out_path}")
    return out_path


def main():
    plt_handle = setup_plotting()
    df_tumor = load_tumor_data()

    # Print per-project breakdown
    for proj in ALL_PROJECTS:
        n = (df_tumor["project_short_name"] == proj).sum()
        print(f"  {proj}: {n} tumor samples")

    # Run statistical tests
    print("\n" + "=" * 60)
    print("  Running subtype classification and statistical tests")
    print("=" * 60)
    results_df = run_subtype_tests(df_tumor)

    if results_df.empty:
        print("No results to process.")
        return

    results_df = apply_fdr(results_df)

    # Print summary table
    print("\n-- Subtype Circadian CV Summary (FDR-corrected) --")
    for _, row in results_df.iterrows():
        kw_fdr = row.get("kruskal_fdr", np.nan)
        mw_fdr = row.get("mw_am_vs_dc_fdr", np.nan)
        print(
            f"  {row['project_label']:15s}  "
            f"AM={row['cv_mean_active_masking']:.3f}  "
            f"DC={row['cv_mean_decoherence']:.3f}  "
            f"Mx={row['cv_mean_mixed']:.3f}  "
            f"KW-FDR={kw_fdr:.2e}  "
            f"MW(AM-DC)-FDR={mw_fdr:.2e}"
        )

    # Save CSV
    csv_path = os.path.join(RESULTS_CSV_DIR, "immune_subtype_comparison.csv")
    results_df.to_csv(csv_path, index=False)
    print(f"\nResults saved to: {csv_path}")

    # Generate figures
    print("\nGenerating figures...")
    fig1_path = plot_circadian_cv_by_subtype(df_tumor, plt_handle)
    fig2_path = plot_boundary_scatter(df_tumor, plt_handle)

    # Upload
    print("\nUploading to GCS...")
    try:
        upload_to_gcs(csv_path, "analysis/immune_subtype_comparison.csv")
        upload_to_gcs(fig1_path, "figures/immune_subtype_circadian_cv.png")
        if fig2_path:
            upload_to_gcs(
                fig2_path, "figures/immune_subtype_boundary_scatter.png"
            )
    except Exception as e:
        print(f"  GCS upload failed (non-fatal): {e}")

    print("\nDone.")


if __name__ == "__main__":
    main()
