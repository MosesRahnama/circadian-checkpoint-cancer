# Supplementary Upload Package

This directory contains the organized supplementary materials prepared for journal submission. Files are numbered for clear cross-referencing with the manuscript.

---

## Directory Structure

```
supplementary_upload/
|-- main_figures/     <-- 16 numbered main-text figures
|-- figures/          <-- 8 numbered supplementary figures
|-- data/             <-- 21 numbered supplementary data files
```

---

## Main Figures (`main_figures/`)

These are the 16 figures referenced in the main text of the manuscript, numbered to match figure references.

| File | Figure | Description |
|------|--------|-------------|
| `1_multicancer_circadian_cv_forest.png` | Fig 1 | Core finding: circadian CV-PD-L1 coupling across six cancer types |
| `2_multicancer_correlation_heatmap.png` | Fig 2 | Cross-cancer correlation structure for key gene pairs |
| `3_survival_boundary_failure.png` | Fig 3 | Kaplan-Meier survival by boundary-failure subtype |
| `4_tumor_normal_circadian_cv_paired.png` | Fig 4 | Tumor vs matched normal circadian CV |
| `5_tumor_normal_circadian_cv_summary.png` | Fig 5 | Cohort-level tumor-normal CV shift summary |
| `6_immune_subtype_prevalence_stacked.png` | Fig 6 | AM/DC/Mixed subtype composition across cohorts |
| `7_immune_subtype_circadian_cv.png` | Fig 7 | Circadian CV distributions by subtype |
| `8_immune_subtype_boundary_scatter.png` | Fig 8 | Checkpoint-clock boundary space scatter |
| `9_survival_test_qvalue_summary.png` | Fig 9 | FDR-adjusted survival evidence summary |
| `10_robustness_primary_hr_forest.png` | Fig 10 | Hazard-ratio forest plot for robustness models |
| `11_rmst_am_vs_dc.png` | Fig 11 | RMST difference at four time horizons |
| `12_stage_circadian_cv_boxplot.png` | Fig 12 | Circadian CV by pathologic tumor stage |
| `13_immune_residualization_scatter.png` | Fig 13 | CV vs PD-L1 after immune residualization |
| `14_purity_stratified_correlations.png` | Fig 14 | Purity-stratified CV-PD-L1 correlations |
| `15_external_validation_geo.png` | Fig 15 | GSE91061 external validation |
| `16_thorsson_benchmarking.png` | Fig 16 | Thorsson immune subtype benchmarking and CCD |

---

## Supplementary Figures (`figures/`)

| File | Figure | Description |
|------|--------|-------------|
| `S1_hypothesis1_TCGA_SKCM.png` | S1 | Discovery cohort detail: melanoma (SKCM) |
| `S2_hypothesis1_TCGA_LUAD.png` | S2 | Discovery cohort detail: lung adenocarcinoma (LUAD) |
| `S3_survival_circadian_quartile.png` | S3 | CV quartile survival analysis (null result) |
| `S4_threshold_sensitivity_sweep.png` | S4 | AM/DC threshold sensitivity sweep |
| `S5_convergent_validity_ncv.png` | S5 | nCV convergent validity across six cohorts |
| `S6_observability_index_survival.png` | S6 | Independent observability index survival |
| `S7_control_budget_combined.png` | S7 | Proliferation vs circadian coherence |
| `S8_external_validation_gse78220.png` | S8 | GSE78220 melanoma anti-PD-1 replication |

---

## Supplementary Data (`data/`)

All CSV files are numbered to match supplementary data references in the manuscript.

| File | Description | Rows |
|------|-------------|------|
| `01_tcga_multicancer_correlations.csv` | All Spearman correlations (6 cohorts x 17 comparisons) | 102 |
| `02_survival_logrank_results.csv` | All 12 pre-specified survival tests with FDR | 12 |
| `03_robustness_primary_tests.csv` | Cox model results and PH diagnostics | 8 |
| `04_robustness_cox_coefficients.csv` | Full coefficient tables for all Cox models | ~200 |
| `05_robustness_partial_correlations.csv` | Partial correlations after covariate adjustment | ~10 |
| `06_robustness_ph_diagnostics.csv` | Schoenfeld residual PH assumption tests | ~200 |
| `07_rmst_results.csv` | RMST at 4 horizons with bootstrap CIs | 8 |
| `08_immune_residualization_results.csv` | Residualized and purity-stratified results | ~16 |
| `09_immune_subtype_comparison.csv` | AM/DC prevalence and CV statistics | 6 |
| `10_observability_index_results.csv` | Independent observability index results | 6 |
| `11_external_validation_results.csv` | GSE91061 replication results | 1 |
| `12_external_validation_gse78220.csv` | GSE78220 replication results | 1 |
| `13_external_validation_gse115821.csv` | GSE115821 replication results | 1 |
| `14_tumor_normal_comparison.csv` | Tumor vs normal comparisons | ~35 |
| `15_threshold_sensitivity_results.csv` | Threshold sweep (22 tests) | 22 |
| `16_convergent_validity_results.csv` | nCV convergent validity (6 cohorts) | 6 |
| `17_thorsson_benchmarking_results.csv` | Thorsson subtype benchmarking and CCD | ~30 |
| `18_master_fdr_results.csv` | Global FDR across all 244+ tests | 244 |
| `19_evidence_ledger_prespecified.csv` | 32-endpoint evidence audit | 32 |
| `20_stage_analysis_results.csv` | Stage-stratified analysis results | ~48 |
| `21_hypothesis1_correlations.csv` | Discovery cohort detailed correlations | ~50 |
