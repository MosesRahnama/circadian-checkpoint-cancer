#!/usr/bin/env python3
"""
Thorsson Immune Subtype Benchmarking + CCD Convergent Validity
==============================================================
Tests whether AM/DC subtypes add information beyond Thorsson C1-C6 immune
subtypes, and computes Clock Correlation Distance (CCD) as a second
convergent validity metric alongside nCV.

Analyses:
1. AM/DC distribution across Thorsson C1-C6 (is AM just C2?)
2. CV-PD-L1 correlation WITHIN each Thorsson immune subtype
3. AM/DC survival stratification WITHIN Thorsson C2 (IFN-gamma dominant)
4. Multivariable model: does CV predict PD-L1 beyond immune subtype?
5. CCD (Clock Correlation Distance) - Shilts et al. 2018 method
6. CCD-PD-L1 correlation as independent replication of CV-PD-L1 finding

Output:
  results/csv/thorsson_benchmarking_results.csv
  results/figures/thorsson_benchmarking.png
"""
import os
import sys
sys.path.insert(0, os.path.dirname(__file__))

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from tcga_config import (
    DATA_DIR, RESULTS_CSV_DIR, CIRCADIAN, ALL_PROJECTS, PROJECT_LABELS,
    log_transform, compute_circadian_cv, classify_boundary_failure
)

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_PAPER_ROOT = os.path.dirname(SCRIPT_DIR)
FIG_DIR = os.path.join(_PAPER_ROOT, 'results', 'figures')
os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(RESULTS_CSV_DIR, exist_ok=True)

# Load data
print("Loading data...")
expr = pd.read_csv(os.path.join(DATA_DIR, 'tcga_expanded_tpm.csv'))
cov = pd.read_csv(os.path.join(DATA_DIR, 'tcga_purity_immune_covariates.csv'))

# Filter to primary tumors
expr = expr[expr['sample_type_name'] == 'Primary Tumor'].copy()

# Compute circadian CV and boundary-failure subtypes
expr['Circadian_CV'] = compute_circadian_cv(expr)
expr['boundary_mode'] = classify_boundary_failure(expr)

# Merge with Thorsson covariates
df = expr.merge(cov, on=['case_barcode', 'project_short_name'], how='inner')
df['CD274_log'] = log_transform(df['CD274'])

# Log-transform clock genes
for g in CIRCADIAN:
    df[f'{g}_log'] = log_transform(df[g])

print(f"Merged dataset: {len(df)} samples with Thorsson immune subtypes")
print(f"Immune subtype distribution:\n{df['immune_subtype'].value_counts().sort_index()}\n")

results = []

# ============================================================
# Part 1: AM/DC distribution across Thorsson C1-C6
# ============================================================
print("=" * 60)
print("  Part 1: AM/DC vs Thorsson Immune Subtypes")
print("=" * 60)

cross = pd.crosstab(df['boundary_mode'], df['immune_subtype'], margins=True)
print("\nCross-tabulation (AM/DC x Thorsson C1-C6):")
print(cross.to_string())

# Chi-squared test
ct = pd.crosstab(df['boundary_mode'], df['immune_subtype'])
chi2, p_chi2, dof, _ = stats.chi2_contingency(ct)
print(f"\nChi-squared: {chi2:.1f}, dof={dof}, p={p_chi2:.2e}")
print("  (If significant: AM/DC is NOT uniformly distributed across immune subtypes)")

# What fraction of AM is in C2?
am_in_c2 = cross.loc['Active_Masking', 'C2'] if 'C2' in cross.columns else 0
am_total = cross.loc['Active_Masking', 'All']
print(f"\nAM in C2: {am_in_c2}/{am_total} ({100*am_in_c2/am_total:.1f}%)")
dc_in_c2 = cross.loc['Decoherence', 'C2'] if 'C2' in cross.columns else 0
dc_total = cross.loc['Decoherence', 'All']
print(f"DC in C2: {dc_in_c2}/{dc_total} ({100*dc_in_c2/dc_total:.1f}%)")

results.append({
    'analysis': 'am_dc_vs_thorsson',
    'test': 'chi_squared',
    'statistic': chi2, 'p_value': p_chi2, 'n': len(df),
    'note': f'AM_in_C2={am_in_c2}/{am_total} ({100*am_in_c2/am_total:.1f}%)'
})

# ============================================================
# Part 2: CV-PD-L1 within each Thorsson subtype
# ============================================================
print("\n" + "=" * 60)
print("  Part 2: CV-PD-L1 Correlation WITHIN Thorsson Subtypes")
print("=" * 60)

for subtype in sorted(df['immune_subtype'].dropna().unique()):
    sub = df[df['immune_subtype'] == subtype]
    if len(sub) < 20:
        continue
    rho, p = stats.spearmanr(sub['Circadian_CV'], sub['CD274_log'])
    sig = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns'
    print(f"  {subtype} (n={len(sub):4d}): rho={rho:+.4f}, p={p:.2e} {sig}")
    results.append({
        'analysis': 'cv_pdl1_within_thorsson',
        'test': f'spearman_{subtype}',
        'statistic': rho, 'p_value': p, 'n': len(sub),
        'note': f'CV-PD-L1 within {subtype}'
    })

# ============================================================
# Part 3: AM/DC survival within C2 (IFN-gamma dominant)
# ============================================================
print("\n" + "=" * 60)
print("  Part 3: AM vs DC Survival WITHIN Thorsson C2")
print("=" * 60)

clin = pd.read_csv(os.path.join(DATA_DIR, 'tcga_clinical.csv'))
# Construct OS and OS.time from clinical data
clin['OS.time'] = clin['days_to_death'].fillna(clin['days_to_last_follow_up'])
clin['OS'] = (clin['vital_status'] == 'Dead').astype(int)
df_surv = df.merge(clin[['case_barcode', 'OS', 'OS.time']].drop_duplicates('case_barcode'),
                    on='case_barcode', how='inner')
df_surv = df_surv.dropna(subset=['OS', 'OS.time'])
df_surv = df_surv[df_surv['OS.time'] > 0]

for subtype in ['C2', 'C1', 'C3']:
    sub = df_surv[(df_surv['immune_subtype'] == subtype) &
                  (df_surv['boundary_mode'].isin(['Active_Masking', 'Decoherence']))]
    am = sub[sub['boundary_mode'] == 'Active_Masking']
    dc = sub[sub['boundary_mode'] == 'Decoherence']

    if len(am) < 10 or len(dc) < 10:
        print(f"  {subtype}: n_AM={len(am)}, n_DC={len(dc)} -- too small, skipping")
        results.append({
            'analysis': f'survival_within_{subtype}',
            'test': 'logrank',
            'statistic': np.nan, 'p_value': np.nan,
            'n': len(am) + len(dc),
            'note': f'n_AM={len(am)}, n_DC={len(dc)}, underpowered'
        })
        continue

    # Log-rank test (manual implementation)
    from scipy.stats import norm
    times = sorted(set(sub['OS.time']))
    O_am, E_am = 0, 0
    for t in times:
        at_risk_am = (am['OS.time'] >= t).sum()
        at_risk_dc = (dc['OS.time'] >= t).sum()
        at_risk_total = at_risk_am + at_risk_dc
        if at_risk_total == 0:
            continue
        events_am = ((am['OS.time'] == t) & (am['OS'] == 1)).sum()
        events_dc = ((dc['OS.time'] == t) & (dc['OS'] == 1)).sum()
        events_total = events_am + events_dc
        if events_total == 0:
            continue
        expected = at_risk_am * events_total / at_risk_total
        O_am += events_am
        E_am += expected

    # Simplified log-rank z-score
    if E_am > 0:
        z = (O_am - E_am) / np.sqrt(E_am) if E_am > 1 else 0
        p_lr = 2 * norm.sf(abs(z))
    else:
        z, p_lr = 0, 1.0

    sig = '***' if p_lr < 0.001 else '**' if p_lr < 0.01 else '*' if p_lr < 0.05 else 'ns'
    direction = "AM better" if z < 0 else "DC better"
    print(f"  {subtype} (n_AM={len(am)}, n_DC={len(dc)}): z={z:.3f}, p={p_lr:.4f} {sig} ({direction})")

    results.append({
        'analysis': f'survival_within_{subtype}',
        'test': 'logrank',
        'statistic': z, 'p_value': p_lr,
        'n': len(am) + len(dc),
        'note': f'n_AM={len(am)}, n_DC={len(dc)}, {direction}'
    })

# ============================================================
# Part 4: Multivariable: does CV predict PD-L1 beyond immune subtype?
# ============================================================
print("\n" + "=" * 60)
print("  Part 4: Incremental Value of CV Beyond Immune Subtype")
print("=" * 60)

from scipy.stats import pearsonr

for proj in ['TCGA-SKCM', 'TCGA-LUAD', 'TCGA-BRCA', 'TCGA-HNSC']:
    sub = df[(df['project_short_name'] == proj) & df['immune_subtype'].notna()].copy()
    if len(sub) < 50:
        continue

    # Model 1: PD-L1 ~ immune_subtype (R2)
    dummies = pd.get_dummies(sub['immune_subtype'], prefix='C', drop_first=True)
    X1 = dummies.values
    y = sub['CD274_log'].values

    # Manual R2 via correlation
    from numpy.linalg import lstsq
    X1_const = np.column_stack([np.ones(len(y)), X1])
    beta1, _, _, _ = lstsq(X1_const, y, rcond=None)
    pred1 = X1_const @ beta1
    ss_res1 = np.sum((y - pred1)**2)
    ss_tot = np.sum((y - y.mean())**2)
    r2_immune = 1 - ss_res1 / ss_tot

    # Model 2: PD-L1 ~ immune_subtype + CV (R2)
    X2_const = np.column_stack([X1_const, sub['Circadian_CV'].values])
    beta2, _, _, _ = lstsq(X2_const, y, rcond=None)
    pred2 = X2_const @ beta2
    ss_res2 = np.sum((y - pred2)**2)
    r2_full = 1 - ss_res2 / ss_tot

    delta_r2 = r2_full - r2_immune

    # F-test for the added CV term (1 df)
    n = len(y)
    k_full = X2_const.shape[1]
    f_stat = (ss_res1 - ss_res2) / (ss_res2 / (n - k_full))
    from scipy.stats import f as f_dist
    p_f = 1 - f_dist.cdf(f_stat, 1, n - k_full)

    sig = '***' if p_f < 0.001 else '**' if p_f < 0.01 else '*' if p_f < 0.05 else 'ns'
    label = PROJECT_LABELS.get(proj, proj)
    print(f"  {label:12s} (n={n:4d}): R2_immune={r2_immune:.4f}, R2_full={r2_full:.4f}, "
          f"delta_R2={delta_r2:.4f}, F={f_stat:.2f}, p={p_f:.2e} {sig}")

    results.append({
        'analysis': 'incremental_r2',
        'test': f'f_test_{proj}',
        'statistic': f_stat, 'p_value': p_f, 'n': n,
        'note': f'R2_immune={r2_immune:.4f}, R2_full={r2_full:.4f}, delta_R2={delta_r2:.4f}'
    })

# ============================================================
# Part 5: CCD (Clock Correlation Distance) - Shilts et al. 2018
# ============================================================
print("\n" + "=" * 60)
print("  Part 5: CCD (Clock Correlation Distance)")
print("=" * 60)

# CCD = 1 - mean(pairwise Spearman rho between clock genes across samples)
# Per Shilts et al. 2018 (PeerJ)
clock_log_cols = [f'{g}_log' for g in CIRCADIAN]

for proj in ALL_PROJECTS:
    sub = df[df['project_short_name'] == proj][clock_log_cols].dropna()
    corr = sub.corr(method='spearman')
    upper = corr.where(np.triu(np.ones(corr.shape), k=1).astype(bool))
    mean_rho = upper.stack().mean()
    ccd = 1 - mean_rho

    label = PROJECT_LABELS.get(proj, proj)
    print(f"  {label:12s} (n={len(sub):4d}): CCD={ccd:.4f} (mean inter-gene rho={mean_rho:.4f})")

    results.append({
        'analysis': 'ccd_population',
        'test': f'ccd_{proj}',
        'statistic': ccd, 'p_value': np.nan, 'n': len(sub),
        'note': f'mean_inter_gene_rho={mean_rho:.4f}'
    })

# ============================================================
# Part 6: Per-sample CCD proxy vs PD-L1
# ============================================================
print("\n" + "=" * 60)
print("  Part 6: Per-Sample CCD Proxy vs PD-L1")
print("=" * 60)

# Per-sample CCD: for each sample, compute pairwise correlations of its
# clock gene values against the population mean profile
for proj in ALL_PROJECTS:
    sub = df[df['project_short_name'] == proj].copy()
    clock_vals = sub[clock_log_cols].values
    pop_means = clock_vals.mean(axis=0)

    # Per-sample: mean absolute deviation from population means
    deviations = clock_vals - pop_means
    sub['ccd_sample'] = np.mean(np.abs(deviations), axis=1)

    rho, p = stats.spearmanr(sub['ccd_sample'], sub['CD274_log'])
    sig = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns'
    label = PROJECT_LABELS.get(proj, proj)
    print(f"  {label:12s} (n={len(sub):4d}): CCD-PD-L1 rho={rho:+.4f}, p={p:.2e} {sig}")

    results.append({
        'analysis': 'ccd_vs_pdl1',
        'test': f'spearman_{proj}',
        'statistic': rho, 'p_value': p, 'n': len(sub),
        'note': f'CCD_sample vs PD-L1'
    })

    # Also correlate CCD-sample with Paper1 CV
    rho_cv, p_cv = stats.spearmanr(sub['ccd_sample'], sub['Circadian_CV'])
    results.append({
        'analysis': 'ccd_vs_cv',
        'test': f'spearman_{proj}',
        'statistic': rho_cv, 'p_value': p_cv, 'n': len(sub),
        'note': f'CCD_sample vs Paper1_CV (convergent validity)'
    })

# ============================================================
# Save results
# ============================================================
results_df = pd.DataFrame(results)
csv_path = os.path.join(RESULTS_CSV_DIR, 'thorsson_benchmarking_results.csv')
results_df.to_csv(csv_path, index=False)
print(f"\nSaved: {csv_path}")

# ============================================================
# Summary figure
# ============================================================
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Panel A: AM/DC distribution across Thorsson subtypes
ct_plot = pd.crosstab(df['immune_subtype'], df['boundary_mode'], normalize='index')
ct_plot = ct_plot[['Active_Masking', 'Mixed', 'Decoherence']]
ct_plot.plot(kind='bar', stacked=True, ax=axes[0, 0],
             color=['#2ecc71', '#95a5a6', '#e74c3c'], edgecolor='white')
axes[0, 0].set_title('A) AM/DC Distribution by Thorsson Immune Subtype')
axes[0, 0].set_xlabel('Thorsson Immune Subtype')
axes[0, 0].set_ylabel('Proportion')
axes[0, 0].legend(title='Boundary Mode', fontsize=8)
axes[0, 0].tick_params(axis='x', rotation=0)

# Panel B: CV-PD-L1 within Thorsson subtypes
within_data = results_df[results_df['analysis'] == 'cv_pdl1_within_thorsson']
subtypes = within_data['test'].str.replace('spearman_', '')
rhos = within_data['statistic'].values
ns = within_data['n'].values
colors = ['#e74c3c' if p < 0.05 else '#95a5a6' for p in within_data['p_value']]
axes[0, 1].barh(range(len(subtypes)), rhos, color=colors)
axes[0, 1].set_yticks(range(len(subtypes)))
axes[0, 1].set_yticklabels([f'{s} (n={n})' for s, n in zip(subtypes, ns)])
axes[0, 1].set_xlabel('Spearman rho (CV vs PD-L1)')
axes[0, 1].set_title('B) CV-PD-L1 Within Each Thorsson Subtype')
axes[0, 1].axvline(0, color='black', linewidth=0.5)
for i, (r, p) in enumerate(zip(rhos, within_data['p_value'])):
    m = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns'
    axes[0, 1].text(r + 0.01 if r >= 0 else r - 0.08, i, f'{r:.3f} {m}', va='center', fontsize=8)

# Panel C: Incremental R2
incr = results_df[results_df['analysis'] == 'incremental_r2']
labels_c = [t.replace('f_test_TCGA-', '') for t in incr['test']]
r2_immune = [float(n.split('R2_immune=')[1].split(',')[0]) for n in incr['note']]
delta_r2 = [float(n.split('delta_R2=')[1]) for n in incr['note']]
x = range(len(labels_c))
axes[1, 0].bar(x, r2_immune, label='Immune subtype R2', color='#3498db', alpha=0.7)
axes[1, 0].bar(x, delta_r2, bottom=r2_immune, label='+ CV (delta R2)', color='#e67e22', alpha=0.9)
axes[1, 0].set_xticks(x)
axes[1, 0].set_xticklabels(labels_c)
axes[1, 0].set_ylabel('R-squared (PD-L1 prediction)')
axes[1, 0].set_title('C) Incremental R2: CV Beyond Immune Subtype')
axes[1, 0].legend(fontsize=8)
for i, (dr, p) in enumerate(zip(delta_r2, incr['p_value'])):
    m = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns'
    axes[1, 0].text(i, r2_immune[i] + dr + 0.005, f'+{dr:.3f} {m}', ha='center', fontsize=8)

# Panel D: CCD vs PD-L1
ccd_pdl1 = results_df[results_df['analysis'] == 'ccd_vs_pdl1']
labels_d = [t.replace('spearman_TCGA-', '') for t in ccd_pdl1['test']]
rhos_d = ccd_pdl1['statistic'].values
colors_d = ['#e74c3c' if p < 0.05 else '#95a5a6' for p in ccd_pdl1['p_value']]
axes[1, 1].barh(range(len(labels_d)), rhos_d, color=colors_d)
axes[1, 1].set_yticks(range(len(labels_d)))
axes[1, 1].set_yticklabels(labels_d)
axes[1, 1].set_xlabel('Spearman rho')
axes[1, 1].set_title('D) CCD-Sample vs PD-L1 (Independent Clock Metric)')
axes[1, 1].axvline(0, color='black', linewidth=0.5)
for i, (r, p) in enumerate(zip(rhos_d, ccd_pdl1['p_value'])):
    m = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns'
    axes[1, 1].text(r + 0.01 if r >= 0 else r - 0.08, i, f'{r:.3f} {m}', va='center', fontsize=8)

plt.tight_layout()
fig_path = os.path.join(FIG_DIR, 'thorsson_benchmarking.png')
plt.savefig(fig_path, dpi=300, bbox_inches='tight')
print(f"Saved: {fig_path}")

print("\n" + "=" * 60)
print("  DONE: Thorsson benchmarking + CCD analysis complete")
print("=" * 60)
