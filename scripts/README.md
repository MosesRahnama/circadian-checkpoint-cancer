# Analysis Scripts

This directory contains all Python scripts that produce the results reported in the manuscript. Every script reads from `../data/` and writes outputs to `../results/csv/` and `../results/figures/`.

---

## How to Run

### Run everything at once

```bash
cd scripts
python run_all.py
```

Total runtime: ~15-20 minutes. Executes all 14 analysis scripts in the order listed below.

### Run individual scripts

Each script is self-contained and can be run independently:

```bash
cd scripts
python <script_name>.py
```

---

## Master Runner

### `run_all.py`
Orchestrates the full analysis pipeline. Runs all 14 scripts in sequence, reports pass/fail status and total runtime. This is the single entry point for full reproduction.

---

## Shared Configuration

### `tcga_config.py`
Shared configuration file used by all analysis scripts. Defines:
- File paths to source data (`../data/`) and output directories (`../results/csv/`, `../results/figures/`)
- Cohort definitions and sample counts
- Gene lists (clock genes, checkpoint genes, MHC-I genes, etc.)
- Common utility functions for loading data, computing circadian CV, and formatting output

---

## Primary Analysis Scripts

These produce the core results reported in the main text.

### `tcga_multicancer.py`
**Manuscript section:** Results 3.1
**What it does:** Computes Spearman correlations between all target gene pairs across six TCGA cohorts. Tests the primary hypothesis: circadian CV is negatively correlated with PD-L1 (CD274).
**Output CSV:** `tcga_multicancer_correlations.csv` (102 rows: 17 comparisons x 6 cohorts)
**Output figures:** `multicancer_correlation_heatmap.png`, `multicancer_circadian_cv_forest.png`

### `tcga_tumor_normal.py`
**Manuscript section:** Results 3.4
**What it does:** Compares circadian CV between tumor and matched normal tissue using Wilcoxon signed-rank (paired) and Mann-Whitney U (unpaired) tests. Tests the "locked, not broken" hypothesis.
**Output CSV:** `tumor_normal_comparison.csv`
**Output figures:** `tumor_normal_circadian_cv_paired.png`, `tumor_normal_circadian_cv_summary.png`

### `tcga_immune_subtype.py`
**Manuscript section:** Results 3.5
**What it does:** Classifies tumors into Active Masking (AM), Temporal Decoherence (DC), and Mixed subtypes. Computes Kruskal-Wallis and Mann-Whitney statistics for circadian CV differences between subtypes.
**Output CSV:** `immune_subtype_comparison.csv`
**Output figures:** `immune_subtype_circadian_cv.png`, `immune_subtype_boundary_scatter.png`, `immune_subtype_prevalence_stacked.png`

### `tcga_survival.py`
**Manuscript section:** Results 3.6
**What it does:** Kaplan-Meier survival analysis with log-rank tests for 12 pre-specified contrasts (6 AM-vs-DC + 6 circadian CV quartile Q1-vs-Q4). Applies Benjamini-Hochberg FDR correction across all 12 tests.
**Output CSV:** `survival_logrank_results.csv`
**Output figures:** `survival_boundary_failure.png`, `survival_circadian_quartile.png`

### `tcga_stage_analysis.py`
**Manuscript section:** Results 3.8
**What it does:** Tests whether circadian CV varies by pathologic tumor stage (Kruskal-Wallis, Spearman trend). Computes global FDR across all 244+ tests in the study.
**Output CSVs:** `stage_analysis_results.csv`, `master_fdr_results.csv`
**Output figure:** `stage_circadian_cv_boxplot.png`

---

## Sensitivity Analysis Scripts

These test robustness of the primary findings under alternative modeling choices and potential confounders.

### `robustness_check.py`
**Manuscript section:** Results 3.7
**What it does:** Fits multivariable Cox proportional hazards models with clinical and microenvironment covariates (age, sex, stage, purity, immune fractions). Tests 8 models across SKCM and LUAD, including PD-L1 x clock interaction terms. Checks PH assumptions.
**Output CSVs:** `robustness_primary_tests.csv`, `robustness_cox_coefficients.csv`, `robustness_partial_correlations.csv`, `robustness_ph_diagnostics.csv`
**Output figure:** `robustness_primary_hr_forest.png`

### `sensitivity_rmst.py`
**Manuscript section:** Results 3.7
**What it does:** Computes Restricted Mean Survival Time (RMST) as a PH-assumption-free alternative to Cox models. Tests AM vs DC at 4 time horizons (36, 60, 84, 120 months) with bootstrap 95% CIs (n=1,000 resamples, seed=42).
**Output CSV:** `rmst_results.csv`
**Output figure:** `rmst_am_vs_dc.png`

### `sensitivity_threshold_sweep.py`
**Manuscript section:** Results 3.6
**What it does:** Sweeps the AM/DC classification percentile threshold from 25th to 75th in 5% steps. Re-runs log-rank survival tests at each threshold to assess stability. Demonstrates that the SKCM signal is robust across 9/11 thresholds.
**Output CSV:** `threshold_sensitivity_results.csv`
**Output figure:** `threshold_sensitivity_sweep.png`

### `sensitivity_immune_residualization.py`
**Manuscript section:** Results 3.9
**What it does:** Regresses circadian CV on five immune covariates (leukocyte fraction, lymphocyte infiltration, IFN-gamma response, tumor purity, stromal fraction) via OLS. Correlates residuals with PD-L1 to test whether the coupling survives immune-fraction adjustment. Also tests purity-stratified correlations within tertiles.
**Output CSV:** `immune_residualization_results.csv`
**Output figures:** `immune_residualization_scatter.png`, `purity_stratified_correlations.png`

### `composite_observability_index.py`
**Manuscript section:** Results 3.10
**What it does:** Constructs an independent observability index from features deliberately excluded from the AM/DC classification (MHC-I composite without B2M, gap junction composite, orthogonal circadian CV without ARNTL/PER1). Tests whether this independent metric also separates subtypes and associates with survival.
**Output CSV:** `observability_index_results.csv`
**Output figure:** `observability_index_survival.png`

---

## External Validation Scripts

These test the circadian-checkpoint coupling in independent cohorts outside TCGA.

### `run_external_validation.py`
**Manuscript section:** Results 3.11
**What it does:** Tests the CV-PD-L1 coupling in GSE91061 (Riaz et al. 2017, nivolumab melanoma, n=49 pre-treatment samples). Also tests AM/DC response rate comparison (Fisher exact test).
**Output CSV:** `external_validation_results.csv`
**Output figure:** `external_validation_geo.png`

### `external_validation_geo.py`
Helper script for downloading and processing GEO datasets. Used by `run_external_validation.py`.

### `external_validation_gse78220.py`
**Manuscript section:** Results 3.11
**What it does:** Tests the CV-PD-L1 coupling in GSE78220 (Hugo et al. 2016, anti-PD-1 melanoma, n=27 baseline samples).
**Output CSV:** `external_validation_gse78220.csv`
**Output figure:** `external_validation_gse78220.png`

### `external_validation_gse136961.py`
Attempted external validation on GSE136961. Included for completeness.

---

## Convergent Validity Scripts

These test whether the circadian CV metric aligns with established clock-disruption measures from the literature.

### `convergent_validity_ncv.py`
**Manuscript section:** Methods (convergent validity paragraph)
**What it does:** Computes a per-sample nCV-weighted score based on the Wu & Hogenesch (2021) normalized coefficient of variation method. Correlates it with our within-sample CV to establish that both metrics partially capture the same clock-structure variation. Reports Spearman correlations in all 6 TCGA cohorts.
**Output CSV:** `convergent_validity_results.csv`
**Output figure:** `convergent_validity_ncv.png`

### `thorsson_benchmarking.py`
**Manuscript section:** Results 3.12
**What it does:** Three analyses in one:
1. Maps each tumor to Thorsson immune subtypes (C1-C6) and tests whether CV-PD-L1 coupling persists within each subtype.
2. Tests incremental variance (delta R-squared) of circadian CV beyond immune subtype for predicting PD-L1.
3. Computes a per-sample Clock Correlation Distance (CCD) proxy based on Shilts et al. (2018) and tests whether CCD also shows the negative PD-L1 coupling.
**Output CSV:** `thorsson_benchmarking_results.csv`
**Output figure:** `thorsson_benchmarking.png`

---

## Diagnostic Scripts

### `check_cv_mean_artifact.py`
Tests whether the CV-PD-L1 coupling is an artifact of mean clock-gene expression level. Computes partial correlations after residualizing both CV and PD-L1 on mean log2 clock-gene expression.

### `check_partial_correlation.py`
General partial correlation utility for testing associations after controlling for covariates.

---

## Cached External Data

The `geo_cache/` subdirectory contains pre-downloaded external validation data:

| File | Description |
|------|-------------|
| `gse91061_sample_meta.csv` | Sample metadata for GSE91061 (Riaz et al. 2017, nivolumab melanoma) |
| `gse91061_target_expression.csv` | Target gene expression for GSE91061 (49 pre-treatment samples) |
| `GSE78220_PatientFPKM.xlsx` | FPKM expression matrix for GSE78220 (Hugo et al. 2016, anti-PD-1 melanoma) |
| `GSE115821_MGH_counts.csv.gz` | Count matrix for GSE115821 (Auslander et al. 2018, melanoma) |
