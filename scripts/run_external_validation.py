"""
External Validation: GSE91061 (Riaz et al. 2017, Nivolumab Melanoma)
====================================================================
Uses pre-downloaded FPKM data and XML metadata to test the framework's
prospective prediction in an independent immunotherapy cohort.

Output:
  - experiments/tcga/external_validation_results.csv
  - results/external_validation_geo.png
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

warnings.filterwarnings("ignore")

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_PAPER_ROOT = os.path.dirname(SCRIPT_DIR)
RESULTS_DIR = os.path.join(_PAPER_ROOT, "results", "figures")
CACHE_DIR = os.path.join(SCRIPT_DIR, "geo_cache")
TARGET_EXPR_PATH = os.path.join(CACHE_DIR, "gse91061_target_expression.csv")
TARGET_META_PATH = os.path.join(CACHE_DIR, "gse91061_sample_meta.csv")


def ensure_cache_files():
    """Ensure required GEO cache files exist; rebuild if missing."""
    missing = [p for p in [TARGET_EXPR_PATH, TARGET_META_PATH] if not os.path.exists(p)]
    if not missing:
        return

    print("Missing GEO cache files detected; attempting to rebuild cache.")
    from external_validation_geo import prepare_gse91061_cache

    prepare_gse91061_cache(force_refresh=False, allow_download=True)

    still_missing = [p for p in [TARGET_EXPR_PATH, TARGET_META_PATH] if not os.path.exists(p)]
    if still_missing:
        raise FileNotFoundError(
            "Could not prepare required GEO cache files: "
            + ", ".join(still_missing)
        )


def main():
    print("=" * 65)
    print("  External Validation: GSE91061 (Nivolumab Melanoma)")
    print("=" * 65)

    ensure_cache_files()

    # Load expression and metadata (committed by default; rebuilt when missing)
    expr = pd.read_csv(TARGET_EXPR_PATH, index_col=0)
    meta = pd.read_csv(TARGET_META_PATH, index_col=0)

    required_meta_cols = {"response", "timepoint", "patient"}
    missing_meta_cols = [c for c in required_meta_cols if c not in meta.columns]
    if missing_meta_cols:
        raise ValueError(
            "Missing expected columns in gse91061_sample_meta.csv: "
            + ", ".join(missing_meta_cols)
        )

    # Response annotations are now embedded in the committed sample_meta CSV.
    # No XML parsing needed; the CSV was enriched from the GEO XML during
    # initial data preparation (see external_validation_geo.py).
    expr["response"] = meta["response"]
    expr["timepoint"] = meta["timepoint"]
    expr["patient"] = meta["patient"]

    # Pre-treatment samples with known response
    pre = expr[expr["timepoint"] == "Pre"].copy()
    pre = pre[pre["response"].isin(["PD", "SD", "PRCR"])].copy()
    print(f"Pre-treatment samples with known response: {len(pre)}")
    print(f"Response distribution: {pre['response'].value_counts().to_dict()}")

    # Compute circadian CV
    clock_genes = ["ARNTL", "CLOCK", "PER1", "PER2", "CRY1", "CRY2"]
    clock_log = np.log2(pre[clock_genes] + 1)
    pre["circ_cv"] = clock_log.std(axis=1) / clock_log.mean(axis=1)
    pre["cd274_log"] = np.log2(pre["CD274"] + 1)
    pre["arntl_log"] = np.log2(pre["ARNTL"] + 1)
    pre["per1_log"] = np.log2(pre["PER1"] + 1)
    pre["b2m_log"] = np.log2(pre["B2M"] + 1)

    # === Test 1: Circadian CV - PD-L1 correlation ===
    rho_cv, p_cv = stats.spearmanr(pre["circ_cv"], pre["cd274_log"])
    print(f"\n=== Test 1: Circadian CV vs PD-L1 ===")
    direction = "NEGATIVE (replicates TCGA)" if rho_cv < 0 else "POSITIVE (does NOT replicate)"
    sig = "***" if p_cv < 0.001 else "**" if p_cv < 0.01 else "*" if p_cv < 0.05 else "ns"
    print(f"  Spearman rho = {rho_cv:.3f}, p = {p_cv:.4f} {sig}")
    print(f"  Direction: {direction}")

    # === Test 2: Classify AM/DC ===
    cd274_med = pre["cd274_log"].median()
    arntl_med = pre["arntl_log"].median()
    per1_med = pre["per1_log"].median()
    b2m_med = pre["b2m_log"].median()

    pre["bf_class"] = "Mixed"
    am_mask = (
        (pre["cd274_log"] > cd274_med)
        & (pre["arntl_log"] > arntl_med)
        & (pre["per1_log"] < per1_med)
    )
    dc_mask = (pre["cd274_log"] < cd274_med) & (pre["b2m_log"] < b2m_med)
    pre.loc[am_mask, "bf_class"] = "Active_Masking"
    pre.loc[dc_mask, "bf_class"] = "Decoherence"

    print(f"\n=== Test 2: AM/DC Classification ===")
    print(f"  Active Masking: {am_mask.sum()}")
    print(f"  Decoherence: {dc_mask.sum()}")
    print(f"  Mixed: {(~am_mask & ~dc_mask).sum()}")

    # === Test 3: Response rate by subtype ===
    pre["responder"] = pre["response"] == "PRCR"
    pre["progressor"] = pre["response"] == "PD"

    print(f"\n=== Test 3: Response by Boundary-Failure Subtype ===")
    for bf in ["Active_Masking", "Decoherence", "Mixed"]:
        sub = pre[pre["bf_class"] == bf]
        if len(sub) == 0:
            continue
        resp_rate = sub["responder"].mean() * 100
        prog_rate = sub["progressor"].mean() * 100
        print(f"  {bf:18s}: n={len(sub):2d}  response={resp_rate:.1f}%  progression={prog_rate:.1f}%")

    am_resp = pre[pre["bf_class"] == "Active_Masking"]
    dc_resp = pre[pre["bf_class"] == "Decoherence"]

    fisher_or, fisher_p = np.nan, np.nan
    if len(am_resp) >= 3 and len(dc_resp) >= 3:
        table = np.array([
            [am_resp["responder"].sum(), (~am_resp["responder"]).sum()],
            [dc_resp["responder"].sum(), (~dc_resp["responder"]).sum()],
        ])
        fisher_or, fisher_p = stats.fisher_exact(table)
        print(f"\n  Fisher exact (AM vs DC, responder): OR={fisher_or:.2f}, p={fisher_p:.4f}")

    # === Test 4: Circadian CV by response ===
    print(f"\n=== Test 4: Circadian CV by Response ===")
    for resp in ["PRCR", "SD", "PD"]:
        sub = pre[pre["response"] == resp]
        print(f"  {resp}: n={len(sub)}, CV mean={sub['circ_cv'].mean():.3f}")

    # Save results
    results = {
        "dataset": "GSE91061 (Riaz et al. 2017)",
        "n_pre_treatment": len(pre),
        "cv_pdl1_rho": round(rho_cv, 4),
        "cv_pdl1_p": round(p_cv, 6),
        "replicates_tcga_direction": rho_cv < 0,
        "n_am": int(am_mask.sum()),
        "n_dc": int(dc_mask.sum()),
        "am_response_pct": round(am_resp["responder"].mean() * 100, 1) if len(am_resp) > 0 else None,
        "dc_response_pct": round(dc_resp["responder"].mean() * 100, 1) if len(dc_resp) > 0 else None,
        "fisher_or": round(fisher_or, 3) if np.isfinite(fisher_or) else None,
        "fisher_p": round(fisher_p, 4) if np.isfinite(fisher_p) else None,
    }
    res_df = pd.DataFrame([results])
    csv_path = os.path.join(_PAPER_ROOT, "results", "csv", "external_validation_results.csv")
    res_df.to_csv(csv_path, index=False)
    print(f"\nSaved: {csv_path}")

    # === FIGURE ===
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    # Panel 1: CV vs PD-L1
    ax = axes[0]
    ax.scatter(pre["cd274_log"], pre["circ_cv"], s=25, alpha=0.5, c="#2c7fb8")
    z = np.polyfit(pre["cd274_log"].values, pre["circ_cv"].values, 1)
    x_line = np.linspace(pre["cd274_log"].min(), pre["cd274_log"].max(), 100)
    ax.plot(x_line, np.polyval(z, x_line), "r--", linewidth=1.5)
    ax.set_xlabel("log2(PD-L1 + 1)")
    ax.set_ylabel("Circadian CV")
    ax.set_title(f"Circadian CV vs PD-L1\nrho={rho_cv:.3f}, p={p_cv:.3f} {sig}")
    ax.grid(alpha=0.2)

    # Panel 2: Response rate by subtype
    ax = axes[1]
    subtypes = ["Active_Masking", "Mixed", "Decoherence"]
    resp_rates = []
    ns = []
    for bf in subtypes:
        sub = pre[pre["bf_class"] == bf]
        resp_rates.append(sub["responder"].mean() * 100 if len(sub) > 0 else 0)
        ns.append(len(sub))

    colors_bf = ["#1b9e77", "#7570b3", "#d95f02"]
    ax.bar(range(len(subtypes)), resp_rates, color=colors_bf, alpha=0.85,
           edgecolor="black", linewidth=0.5)
    ax.set_xticks(range(len(subtypes)))
    ax.set_xticklabels([s.replace("_", " ") for s in subtypes], fontsize=9)
    ax.set_ylabel("Response Rate (PR/CR, %)")
    ax.set_title("Response by Boundary Mode\n(Pre-treatment)")
    for i, (r, n) in enumerate(zip(resp_rates, ns)):
        ax.text(i, r + 1.5, f"n={n}", ha="center", fontsize=9)
    ax.grid(axis="y", alpha=0.2)

    # Panel 3: CV by response
    ax = axes[2]
    resp_order = ["PRCR", "SD", "PD"]
    resp_labels = ["PR/CR", "SD", "PD"]
    resp_colors = ["#1b9e77", "#7570b3", "#d95f02"]
    data_by_resp = [pre[pre["response"] == r]["circ_cv"].values for r in resp_order]
    bp = ax.boxplot(data_by_resp, tick_labels=resp_labels, patch_artist=True, showfliers=False)
    for patch, color in zip(bp["boxes"], resp_colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    ax.set_ylabel("Circadian CV")
    ax.set_title("Circadian CV by\nImmunotherapy Response")
    ax.grid(axis="y", alpha=0.2)

    fig.suptitle(
        "External Validation: GSE91061 (Riaz et al. 2017, Nivolumab Melanoma, n="
        + str(len(pre)) + " pre-treatment)",
        fontsize=12, fontweight="bold",
    )
    fig.tight_layout()
    fig_path = os.path.join(RESULTS_DIR, "external_validation_geo.png")
    fig.savefig(fig_path, dpi=300, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {fig_path}")

    print("\nExternal validation complete.")


if __name__ == "__main__":
    main()
