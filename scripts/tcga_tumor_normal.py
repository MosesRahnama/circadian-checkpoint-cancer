"""
Tumor vs Matched Normal Comparison
====================================
Compares circadian CV, gap junction expression, PD-L1, and MHC-I
between matched tumor and normal tissue samples across cancer types
with sufficient normal-tissue representation (>= 20 normals).

Statistical tests:
  - Wilcoxon signed-rank (paired)
  - Mann-Whitney U (unpaired)
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
    log_transform, compute_circadian_cv,
    DATA_DIR, FIGURE_DIR,
    setup_plotting, upload_to_gcs,
)

# ── Constants ─────────────────────────────────────────────────────────────
CSV_PATH = os.path.join(DATA_DIR, r"tcga_expanded_tpm.csv")
MIN_NORMALS = 20

# Cancer types eligible for tumor-vs-normal (exclude SKCM, only 1 normal)
ELIGIBLE_PROJECTS = [
    "TCGA-LUAD",
    "TCGA-BRCA",
    "TCGA-COAD",
    "TCGA-HNSC",
    "TCGA-LUSC",
]

# Analytes to compare between tumor and normal
COMPARISON_GENES = ["CD274", "B2M"] + GAP_JUNCTION


def load_data():
    """Load expression data and print summary."""
    df = pd.read_csv(CSV_PATH)
    print(f"Loaded {len(df)} samples from {CSV_PATH}")
    return df


def get_matched_pairs(df, project):
    """
    For a given project, find patients (case_barcodes) that have both
    a Solid Tissue Normal sample and a Tumor/Metastatic sample.
    Returns (tumor_df, normal_df) aligned by case_barcode.
    """
    proj_df = df[df["project_short_name"] == project].copy()

    normal = proj_df[proj_df["sample_type_name"] == "Solid Tissue Normal"]
    tumor = proj_df[proj_df["sample_type_name"].str.contains(
        "Tumor|Metastatic", case=False, na=False
    )].copy()

    # Find shared case barcodes
    shared_cases = sorted(
        set(normal["case_barcode"]) & set(tumor["case_barcode"])
    )

    if len(shared_cases) < 5:
        print(f"  {project}: only {len(shared_cases)} matched pairs -- skipping")
        return None, None, None, None, []

    # For patients with multiple tumor samples, pick a deterministic
    # biologically-prioritized sample (Primary > Recurrent > Metastatic).
    tumor["_sample_priority"] = 4
    sample_type = tumor["sample_type_name"].fillna("").astype(str)
    tumor.loc[sample_type.str.contains("Primary Tumor", case=False), "_sample_priority"] = 1
    tumor.loc[sample_type.str.contains("Recurrent Tumor", case=False), "_sample_priority"] = 2
    tumor.loc[sample_type.str.contains("Metastatic", case=False), "_sample_priority"] = 3

    normal_matched = (
        normal[normal["case_barcode"].isin(shared_cases)]
        .sort_values(["case_barcode", "sample_barcode"])
        .drop_duplicates(subset="case_barcode", keep="first")
        .set_index("case_barcode")
        .loc[shared_cases]
    )
    tumor_matched = (
        tumor[tumor["case_barcode"].isin(shared_cases)]
        .sort_values(["case_barcode", "_sample_priority", "sample_barcode"])
        .drop_duplicates(subset="case_barcode", keep="first")
        .set_index("case_barcode")
        .loc[shared_cases]
        .drop(columns=["_sample_priority"])
    )

    return tumor_matched, normal_matched, tumor, normal, shared_cases


def run_paired_tests(tumor_df, normal_df, tumor_unpaired_df, normal_unpaired_df, shared_cases):
    """
    Run Wilcoxon signed-rank on matched pairs and Mann-Whitney U
    on unpaired project-level tumor vs normal samples.
    """
    results = []

    # -- Circadian CV --
    tumor_cv = compute_circadian_cv(tumor_df).values
    normal_cv = compute_circadian_cv(normal_df).values

    # Drop NaN pairs
    valid = np.isfinite(tumor_cv) & np.isfinite(normal_cv)
    tumor_cv_unpaired = compute_circadian_cv(tumor_unpaired_df).values
    normal_cv_unpaired = compute_circadian_cv(normal_unpaired_df).values
    if valid.sum() >= 5:
        stat_w, p_w = stats.wilcoxon(tumor_cv[valid], normal_cv[valid])
    else:
        stat_w, p_w = np.nan, np.nan

    if np.isfinite(tumor_cv_unpaired).sum() >= 5 and np.isfinite(normal_cv_unpaired).sum() >= 5:
        stat_u, p_u = stats.mannwhitneyu(
            tumor_cv_unpaired[np.isfinite(tumor_cv_unpaired)],
            normal_cv_unpaired[np.isfinite(normal_cv_unpaired)],
            alternative="two-sided",
        )
    else:
        stat_u, p_u = np.nan, np.nan

    if valid.sum() >= 5 or (
        np.isfinite(tumor_cv_unpaired).sum() >= 5 and np.isfinite(normal_cv_unpaired).sum() >= 5
    ):
        results.append({
            "analyte": "Circadian_CV",
            "n_pairs": int(valid.sum()),
            "n_tumor_unpaired": int(np.isfinite(tumor_cv_unpaired).sum()),
            "n_normal_unpaired": int(np.isfinite(normal_cv_unpaired).sum()),
            "tumor_mean": float(np.nanmean(tumor_cv[valid])),
            "normal_mean": float(np.nanmean(normal_cv[valid])),
            "tumor_median": float(np.nanmedian(tumor_cv[valid])),
            "normal_median": float(np.nanmedian(normal_cv[valid])),
            "wilcoxon_stat": float(stat_w),
            "wilcoxon_p": float(p_w),
            "mannwhitney_stat": float(stat_u),
            "mannwhitney_p": float(p_u),
        })

    # -- Individual genes --
    for gene in COMPARISON_GENES:
        if gene not in tumor_df.columns:
            continue
        t_vals = log_transform(tumor_df[gene]).values
        n_vals = log_transform(normal_df[gene]).values
        valid = np.isfinite(t_vals) & np.isfinite(n_vals)
        if valid.sum() >= 5:
            stat_w, p_w = stats.wilcoxon(t_vals[valid], n_vals[valid])
        else:
            stat_w, p_w = np.nan, np.nan

        t_vals_unpaired = log_transform(tumor_unpaired_df[gene]).values
        n_vals_unpaired = log_transform(normal_unpaired_df[gene]).values
        valid_t_unpaired = np.isfinite(t_vals_unpaired)
        valid_n_unpaired = np.isfinite(n_vals_unpaired)
        if valid_t_unpaired.sum() >= 5 and valid_n_unpaired.sum() >= 5:
            stat_u, p_u = stats.mannwhitneyu(
                t_vals_unpaired[valid_t_unpaired],
                n_vals_unpaired[valid_n_unpaired],
                alternative="two-sided",
            )
        else:
            stat_u, p_u = np.nan, np.nan

        if not np.isfinite(p_w) and not np.isfinite(p_u):
            continue

        results.append({
            "analyte": gene,
            "n_pairs": int(valid.sum()),
            "n_tumor_unpaired": int(valid_t_unpaired.sum()),
            "n_normal_unpaired": int(valid_n_unpaired.sum()),
            "tumor_mean": float(np.nanmean(t_vals[valid])),
            "normal_mean": float(np.nanmean(n_vals[valid])),
            "tumor_median": float(np.nanmedian(t_vals[valid])),
            "normal_median": float(np.nanmedian(n_vals[valid])),
            "wilcoxon_stat": float(stat_w),
            "wilcoxon_p": float(p_w),
            "mannwhitney_stat": float(stat_u),
            "mannwhitney_p": float(p_u),
        })

    return results


def apply_fdr(results_df):
    """Apply Benjamini-Hochberg FDR to Wilcoxon and Mann-Whitney p-values."""
    for col in ["wilcoxon_p", "mannwhitney_p"]:
        pvals = results_df[col].values
        # scipy false_discovery_control returns adjusted p-values
        fdr_col = col.replace("_p", "_fdr")
        if len(pvals) > 0 and np.all(np.isfinite(pvals)):
            results_df[fdr_col] = false_discovery_control(pvals, method="bh")
        else:
            results_df[fdr_col] = np.nan
    return results_df


def sig_label(p):
    """Return significance stars for a p-value."""
    if p < 0.001:
        return "***"
    elif p < 0.01:
        return "**"
    elif p < 0.05:
        return "*"
    else:
        return "ns"


def plot_paired_circadian_cv(all_paired_data, plt_handle):
    """
    Figure 1: 2x3 grid, one panel per cancer type.
    Each panel = boxplot (tumor vs normal) with paired lines.
    """
    fig, axes = plt_handle.subplots(2, 3, figsize=(16, 10))
    fig.suptitle(
        "Circadian CV: Tumor vs Matched Normal",
        fontsize=14, fontweight="bold",
    )
    axes_flat = axes.flatten()

    for idx, (project, data) in enumerate(all_paired_data.items()):
        ax = axes_flat[idx]
        tumor_cv = data["tumor_cv"]
        normal_cv = data["normal_cv"]
        n_pairs = data["n_pairs"]
        p_val = data["wilcoxon_p"]

        # Boxplots
        bp = ax.boxplot(
            [normal_cv, tumor_cv],
            positions=[1, 2],
            widths=0.5,
            patch_artist=True,
            showfliers=False,
        )
        bp["boxes"][0].set_facecolor("#4DBEEE")
        bp["boxes"][0].set_alpha(0.6)
        bp["boxes"][1].set_facecolor("#D95319")
        bp["boxes"][1].set_alpha(0.6)

        # Paired lines
        for i in range(len(normal_cv)):
            ax.plot(
                [1, 2],
                [normal_cv[i], tumor_cv[i]],
                color="gray", alpha=0.15, linewidth=0.5,
            )

        ax.set_xticks([1, 2])
        ax.set_xticklabels(["Normal", "Tumor"])
        ax.set_ylabel("Circadian CV")

        label = PROJECT_LABELS.get(project, project)
        ax.set_title(
            f"{label}\nn={n_pairs}, Wilcoxon p={p_val:.2e} {sig_label(p_val)}",
            fontsize=10,
        )

    # Hide unused panel(s)
    for idx in range(len(all_paired_data), len(axes_flat)):
        axes_flat[idx].set_visible(False)

    plt_handle.tight_layout()
    out_path = os.path.join(FIGURE_DIR, "tumor_normal_circadian_cv_paired.png")
    fig.savefig(out_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt_handle.close(fig)
    print(f"  Saved: {out_path}")
    return out_path


def plot_cross_cancer_summary(all_paired_data, plt_handle):
    """
    Figure 2: Grouped bar chart -- mean circadian CV normal vs tumor
    per cancer type.
    """
    projects = list(all_paired_data.keys())
    labels = [PROJECT_LABELS.get(p, p) for p in projects]
    normal_means = [all_paired_data[p]["normal_mean"] for p in projects]
    tumor_means = [all_paired_data[p]["tumor_mean"] for p in projects]

    x = np.arange(len(projects))
    width = 0.35

    fig, ax = plt_handle.subplots(figsize=(10, 6))
    bars_n = ax.bar(x - width / 2, normal_means, width, label="Normal",
                    color="#4DBEEE", alpha=0.8, edgecolor="black", linewidth=0.5)
    bars_t = ax.bar(x + width / 2, tumor_means, width, label="Tumor",
                    color="#D95319", alpha=0.8, edgecolor="black", linewidth=0.5)

    ax.set_ylabel("Mean Circadian CV")
    ax.set_title(
        "Circadian Coherence: Normal vs Tumor Across Cancer Types",
        fontsize=13, fontweight="bold",
    )
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=10)
    ax.legend(frameon=True)
    ax.set_ylim(0, max(normal_means + tumor_means) * 1.25)

    # Annotate with significance
    for i, proj in enumerate(projects):
        p_val = all_paired_data[proj]["wilcoxon_p"]
        ymax = max(normal_means[i], tumor_means[i])
        ax.text(
            x[i], ymax * 1.05,
            sig_label(p_val),
            ha="center", fontsize=11, fontweight="bold",
        )

    plt_handle.tight_layout()
    out_path = os.path.join(FIGURE_DIR, "tumor_normal_circadian_cv_summary.png")
    fig.savefig(out_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt_handle.close(fig)
    print(f"  Saved: {out_path}")
    return out_path


def main():
    plt_handle = setup_plotting()
    df = load_data()

    all_results = []
    all_paired_data = {}  # for plotting

    for project in ELIGIBLE_PROJECTS:
        proj_df = df[df["project_short_name"] == project]
        n_normals = (proj_df["sample_type_name"] == "Solid Tissue Normal").sum()
        if n_normals < MIN_NORMALS:
            print(f"  Skipping {project}: only {n_normals} normals")
            continue

        print(f"\n{'='*60}")
        label = PROJECT_LABELS.get(project, project)
        print(f"  {project} ({label})")
        print(f"{'='*60}")

        tumor_matched, normal_matched, tumor_unpaired, normal_unpaired, shared_cases = get_matched_pairs(
            df, project
        )
        if tumor_matched is None:
            continue

        print(f"  {len(shared_cases)} matched tumor-normal pairs")

        # Run statistical tests
        test_results = run_paired_tests(
            tumor_matched, normal_matched, tumor_unpaired, normal_unpaired, shared_cases
        )

        for r in test_results:
            r["project"] = project
            r["project_label"] = label

        all_results.extend(test_results)

        # Print summary
        for r in test_results:
            print(
                f"  {r['analyte']:20s}  "
                f"tumor={r['tumor_median']:.3f}  "
                f"normal={r['normal_median']:.3f}  "
                f"Wilcoxon p={r['wilcoxon_p']:.2e} {sig_label(r['wilcoxon_p'])}"
            )

        # Store data for plotting
        tumor_cv = compute_circadian_cv(tumor_matched).values
        normal_cv = compute_circadian_cv(normal_matched).values
        valid = np.isfinite(tumor_cv) & np.isfinite(normal_cv)

        cv_row = next(
            (r for r in test_results if r["analyte"] == "Circadian_CV"), None
        )

        all_paired_data[project] = {
            "tumor_cv": tumor_cv[valid],
            "normal_cv": normal_cv[valid],
            "n_pairs": int(valid.sum()),
            "wilcoxon_p": cv_row["wilcoxon_p"] if cv_row else 1.0,
            "tumor_mean": float(np.nanmean(tumor_cv[valid])),
            "normal_mean": float(np.nanmean(normal_cv[valid])),
        }

    # Combine results and apply FDR
    if not all_results:
        print("No results to process.")
        return

    results_df = pd.DataFrame(all_results)
    results_df = apply_fdr(results_df)

    # Reorder columns
    col_order = [
        "project", "project_label", "analyte", "n_pairs",
        "n_tumor_unpaired", "n_normal_unpaired",
        "tumor_mean", "normal_mean", "tumor_median", "normal_median",
        "wilcoxon_stat", "wilcoxon_p", "wilcoxon_fdr",
        "mannwhitney_stat", "mannwhitney_p", "mannwhitney_fdr",
    ]
    results_df = results_df[[c for c in col_order if c in results_df.columns]]

    # Save CSV
    csv_path = os.path.join(RESULTS_CSV_DIR, "tumor_normal_comparison.csv")
    results_df.to_csv(csv_path, index=False)
    print(f"\nResults saved to: {csv_path}")

    # Print FDR-corrected summary
    print("\n-- FDR-corrected Wilcoxon results (Circadian CV) --")
    cv_rows = results_df[results_df["analyte"] == "Circadian_CV"]
    for _, row in cv_rows.iterrows():
        print(
            f"  {row['project_label']:15s}  "
            f"tumor={row['tumor_median']:.3f}  normal={row['normal_median']:.3f}  "
            f"p={row['wilcoxon_p']:.2e}  FDR={row['wilcoxon_fdr']:.2e}"
        )

    # Generate figures
    print("\nGenerating figures...")
    fig1_path = plot_paired_circadian_cv(all_paired_data, plt_handle)
    fig2_path = plot_cross_cancer_summary(all_paired_data, plt_handle)

    # Upload
    print("\nUploading to GCS...")
    try:
        upload_to_gcs(csv_path, "analysis/tumor_normal_comparison.csv")
        upload_to_gcs(fig1_path, "figures/tumor_normal_circadian_cv_paired.png")
        upload_to_gcs(fig2_path, "figures/tumor_normal_circadian_cv_summary.png")
    except Exception as e:
        print(f"  GCS upload failed (non-fatal): {e}")

    print("\nDone.")


if __name__ == "__main__":
    main()
