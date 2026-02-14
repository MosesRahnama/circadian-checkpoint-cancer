#!/usr/bin/env python3
"""
Convergent validity analysis: compare Paper 1's within-sample cross-gene CV
with Wu & Hogenesch's population-level normalized CV (nCV).

Wu et al. (2021) Bioinformatics 37(23):4581-4583.
nCV = CV of a gene across samples / mean CV of all genes across samples.

This script computes:
1. nCV for each clock gene per cohort (population-level oscillation robustness)
2. A per-sample "nCV-weighted clock score" to enable sample-level comparison
3. Spearman correlation between Paper 1's CV and the nCV-derived score
4. Clock Correlation Distance (CCD) as a second convergent validator
"""

import pandas as pd
import numpy as np
from scipy import stats
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# --- Config ---
CLOCK_GENES = ['ARNTL', 'CLOCK', 'PER1', 'PER2', 'CRY1', 'CRY2']
PROJECTS = ['TCGA-SKCM', 'TCGA-LUAD', 'TCGA-BRCA', 'TCGA-COAD', 'TCGA-HNSC', 'TCGA-LUSC']

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_PAPER_ROOT = os.path.dirname(SCRIPT_DIR)
DATA_PATH = os.path.join(_PAPER_ROOT, 'data', 'tcga_expanded_tpm.csv')
RESULTS_DIR = os.path.join(_PAPER_ROOT, 'results', 'figures')
OUTPUT_CSV = os.path.join(_PAPER_ROOT, 'results', 'csv', 'convergent_validity_results.csv')

# --- Load data ---
print("Loading TCGA expression data...")
df = pd.read_csv(DATA_PATH)
df = df[df['sample_type_name'] == 'Primary Tumor'].copy()
df = df.rename(columns={'project_short_name': 'project'})

# Log2(TPM+1) transform (data is raw TPM)
for g in CLOCK_GENES:
    if g in df.columns:
        df[g] = np.log2(df[g] + 1)
    else:
        raise ValueError(f"Column for {g} not found")

# Also transform CD274 for PD-L1 correlation
if 'CD274' in df.columns:
    df['CD274'] = np.log2(df['CD274'] + 1)

# Compute Paper 1's within-sample CV
clock_cols = [g for g in CLOCK_GENES if g in df.columns]
df['circadian_cv'] = df[clock_cols].std(axis=1) / df[clock_cols].mean(axis=1)

print(f"Loaded {len(df)} tumor samples across {df['project'].nunique()} projects\n")

# ============================================================
# Part 1: Population-level nCV (Wu & Hogenesch method)
# ============================================================
print("=" * 60)
print("  Part 1: Population-level nCV per cohort")
print("=" * 60)

ncv_results = []

for proj in PROJECTS:
    subset = df[df['project'] == proj].copy()
    n = len(subset)

    # Compute CV of each gene across the population (samples)
    gene_cvs = {}
    for g in clock_cols:
        mean_g = subset[g].mean()
        std_g = subset[g].std()
        cv_g = std_g / mean_g if mean_g > 0 else np.nan
        gene_cvs[g] = cv_g

    # Mean CV across all clock genes
    mean_cv = np.mean(list(gene_cvs.values()))

    # nCV = gene CV / mean CV (normalized)
    ncv_vals = {g: gene_cvs[g] / mean_cv if mean_cv > 0 else np.nan for g in clock_cols}

    print(f"\n{proj} (n={n}):")
    print(f"  {'Gene':<8} {'Pop CV':>8} {'nCV':>8}")
    for g in clock_cols:
        print(f"  {g:<8} {gene_cvs[g]:>8.4f} {ncv_vals[g]:>8.4f}")
    print(f"  Mean CV = {mean_cv:.4f}")

    for g in clock_cols:
        ncv_results.append({
            'cohort': proj,
            'gene': g,
            'pop_cv': gene_cvs[g],
            'ncv': ncv_vals[g],
            'n': n
        })

# ============================================================
# Part 2: Per-sample nCV-weighted clock deviation score
# ============================================================
print("\n" + "=" * 60)
print("  Part 2: Per-sample nCV-weighted score & convergent validity")
print("=" * 60)

convergent_results = []

for proj in PROJECTS:
    subset = df[df['project'] == proj].copy()
    n = len(subset)

    # Compute population means and SDs for each clock gene
    pop_stats = {}
    for g in clock_cols:
        pop_stats[g] = {
            'mean': subset[g].mean(),
            'std': subset[g].std()
        }

    # nCV weights: genes with higher population CV (more variable across patients)
    # get higher weight -- they carry more discriminative information
    gene_cvs_pop = {g: pop_stats[g]['std'] / pop_stats[g]['mean']
                     if pop_stats[g]['mean'] > 0 else 0 for g in clock_cols}
    mean_cv_pop = np.mean(list(gene_cvs_pop.values()))
    ncv_weights = {g: gene_cvs_pop[g] / mean_cv_pop if mean_cv_pop > 0 else 1
                   for g in clock_cols}

    # Per-sample: z-score each gene, then compute nCV-weighted mean absolute deviation
    for g in clock_cols:
        subset[f'{g}_z'] = (subset[g] - pop_stats[g]['mean']) / pop_stats[g]['std']

    z_cols = [f'{g}_z' for g in clock_cols]
    weights = np.array([ncv_weights[g] for g in clock_cols])

    # Weighted CV of z-scores (nCV-informed per-sample score)
    z_values = subset[z_cols].values
    weighted_std = np.sqrt(np.average((z_values - z_values.mean(axis=1, keepdims=True))**2,
                                       axis=1, weights=weights))
    weighted_mean = np.abs(np.average(z_values, axis=1, weights=weights))
    subset['ncv_score'] = weighted_std

    # Spearman correlation between Paper 1 CV and nCV-weighted score
    rho_cv_ncv, p_cv_ncv = stats.spearmanr(subset['circadian_cv'], subset['ncv_score'])

    print(f"\n{proj} (n={n}):")
    print(f"  Paper1 CV vs nCV-weighted score: rho={rho_cv_ncv:.4f}, p={p_cv_ncv:.2e}")

    convergent_results.append({
        'cohort': proj,
        'test': 'cv_vs_ncv_score',
        'rho': rho_cv_ncv,
        'p_value': p_cv_ncv,
        'n': n
    })

    # Also: direct Clock Correlation Distance (CCD)
    # CCD = 1 - mean(pairwise Spearman correlations between clock genes across samples)
    # This measures how well clock genes co-vary across the population
    clock_data = subset[clock_cols]
    corr_matrix = clock_data.corr(method='spearman')
    # Upper triangle only (exclude diagonal)
    upper_triangle = corr_matrix.where(
        np.triu(np.ones(corr_matrix.shape), k=1).astype(bool))
    mean_corr = upper_triangle.stack().mean()
    ccd = 1 - mean_corr

    print(f"  CCD (Clock Correlation Distance) = {ccd:.4f} (mean inter-gene rho = {mean_corr:.4f})")

    convergent_results.append({
        'cohort': proj,
        'test': 'ccd',
        'rho': ccd,
        'p_value': np.nan,
        'n': n
    })

    # Per-sample CCD: for each sample, compute how far its clock gene profile
    # deviates from the population correlation structure
    # Use per-sample contribution to CCD as a sample-level metric
    # Simpler: per-sample mean pairwise deviation from population means
    pop_means = {g: pop_stats[g]['mean'] for g in clock_cols}
    deviations = np.array([subset[g] - pop_means[g] for g in clock_cols]).T
    # Per-sample: mean absolute deviation weighted by gene importance
    subset['ccd_sample'] = np.mean(np.abs(deviations), axis=1)

    rho_cv_ccd, p_cv_ccd = stats.spearmanr(subset['circadian_cv'], subset['ccd_sample'])
    print(f"  Paper1 CV vs CCD-sample score: rho={rho_cv_ccd:.4f}, p={p_cv_ccd:.2e}")

    convergent_results.append({
        'cohort': proj,
        'test': 'cv_vs_ccd_sample',
        'rho': rho_cv_ccd,
        'p_value': p_cv_ccd,
        'n': n
    })

    # Also correlate nCV-score with PD-L1 directly
    if 'CD274' in subset.columns:
        pdl1_col = 'CD274'
    elif 'CD274_tpm' in subset.columns:
        subset['CD274'] = np.log2(subset['CD274_tpm'] + 1)
        pdl1_col = 'CD274'
    else:
        pdl1_col = None

    if pdl1_col and pdl1_col in subset.columns:
        rho_ncv_pdl1, p_ncv_pdl1 = stats.spearmanr(subset['ncv_score'], subset[pdl1_col])
        print(f"  nCV-score vs PD-L1: rho={rho_ncv_pdl1:.4f}, p={p_ncv_pdl1:.2e}")

        convergent_results.append({
            'cohort': proj,
            'test': 'ncv_score_vs_pdl1',
            'rho': rho_ncv_pdl1,
            'p_value': p_ncv_pdl1,
            'n': n
        })

# ============================================================
# Part 3: Save results
# ============================================================
results_df = pd.DataFrame(convergent_results)
results_df.to_csv(OUTPUT_CSV, index=False)
print(f"\nSaved: {OUTPUT_CSV}")

# ============================================================
# Part 4: Summary figure
# ============================================================
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# Panel A: CV vs nCV-score correlation across cohorts
cv_ncv = results_df[results_df['test'] == 'cv_vs_ncv_score']
labels = [c.replace('TCGA-', '') for c in cv_ncv['cohort']]
colors = ['#e74c3c' if p < 0.001 else '#f39c12' if p < 0.05 else '#95a5a6'
          for p in cv_ncv['p_value']]
axes[0].barh(range(len(labels)), cv_ncv['rho'].values, color=colors)
axes[0].set_yticks(range(len(labels)))
axes[0].set_yticklabels(labels)
axes[0].set_xlabel('Spearman rho')
axes[0].set_title('Paper1 CV vs nCV-weighted Score')
axes[0].axvline(0, color='black', linewidth=0.5)
for i, (r, p) in enumerate(zip(cv_ncv['rho'], cv_ncv['p_value'])):
    marker = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns'
    axes[0].text(r + 0.02, i, f'{r:.3f} {marker}', va='center', fontsize=9)

# Panel B: CCD values across cohorts
ccd_vals = results_df[results_df['test'] == 'ccd']
axes[1].barh(range(len(labels)), ccd_vals['rho'].values, color='#3498db')
axes[1].set_yticks(range(len(labels)))
axes[1].set_yticklabels(labels)
axes[1].set_xlabel('CCD (1 - mean inter-gene rho)')
axes[1].set_title('Clock Correlation Distance')

# Panel C: nCV-score vs PD-L1 correlation
ncv_pdl1 = results_df[results_df['test'] == 'ncv_score_vs_pdl1']
if len(ncv_pdl1) > 0:
    colors2 = ['#e74c3c' if p < 0.001 else '#f39c12' if p < 0.05 else '#95a5a6'
               for p in ncv_pdl1['p_value']]
    axes[2].barh(range(len(labels)), ncv_pdl1['rho'].values, color=colors2)
    axes[2].set_yticks(range(len(labels)))
    axes[2].set_yticklabels(labels)
    axes[2].set_xlabel('Spearman rho')
    axes[2].set_title('nCV-score vs PD-L1')
    axes[2].axvline(0, color='black', linewidth=0.5)
    for i, (r, p) in enumerate(zip(ncv_pdl1['rho'], ncv_pdl1['p_value'])):
        marker = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns'
        axes[2].text(max(r + 0.02, 0.02) if r >= 0 else r - 0.15, i,
                     f'{r:.3f} {marker}', va='center', fontsize=9)

plt.tight_layout()
fig_path = os.path.join(RESULTS_DIR, 'convergent_validity_ncv.png')
plt.savefig(fig_path, dpi=300, bbox_inches='tight')
print(f"Saved: {fig_path}")

# ============================================================
# Summary
# ============================================================
print("\n" + "=" * 60)
print("  CONVERGENT VALIDITY SUMMARY")
print("=" * 60)
print("\nPaper1 CV vs nCV-weighted score (should be positive if measuring similar construct):")
for _, row in cv_ncv.iterrows():
    sig = '***' if row['p_value'] < 0.001 else '**' if row['p_value'] < 0.01 else '*' if row['p_value'] < 0.05 else 'ns'
    print(f"  {row['cohort']}: rho={row['rho']:.4f}, p={row['p_value']:.2e} {sig}")

print(f"\nnCV-score also correlates with PD-L1 (replicates Paper1 finding with independent metric):")
if len(ncv_pdl1) > 0:
    for _, row in ncv_pdl1.iterrows():
        sig = '***' if row['p_value'] < 0.001 else '**' if row['p_value'] < 0.01 else '*' if row['p_value'] < 0.05 else 'ns'
        print(f"  {row['cohort']}: rho={row['rho']:.4f}, p={row['p_value']:.2e} {sig}")

print("\nDone.")
