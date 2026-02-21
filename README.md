# Circadian Clock-Gene Structure Couples with PD-L1 Across Six Cancer Types and Defines Survival-Associated Subtypes in Melanoma

**Author:** Moses Rahnama, Mina Analytics (moses@minaanalytics.com)
**Date:** February 2026
**Preprint:** [Research Square](https://doi.org/10.21203/rs.3.rs-6176498/v1)

---

## Overview

This repository contains the complete data, code, and outputs for the manuscript *"Circadian Clock-Gene Structure Couples with PD-L1 Across Six Cancer Types and Defines Survival-Associated Subtypes in Melanoma."*

Every statistic, figure, and table in the manuscript can be reproduced from the files in this repository. No external accounts, API keys, or internet access is required. All source data is committed.

**148 statistics in the manuscript have been independently verified against the output CSVs with zero discrepancies.** See [`CHECKER_REPORT.md`](CHECKER_REPORT.md) for the full verification audit.

---

## Summary of Findings

Across six TCGA cancer types (n=3,611 tumors), the within-sample coefficient of variation (CV) across six core clock genes is **negatively correlated** with PD-L1 (CD274) expression in all six cohorts (Spearman rho = -0.125 to -0.381, all FDR q < 0.05). Two boundary-failure subtypes are defined:

- **Active Masking (AM):** High PD-L1 + high BMAL1 + low PER1. The tumor retains a restructured clock and deploys PD-L1 to suppress immune clearance. Associated with better overall survival in melanoma.
- **Temporal Decoherence (DC):** Low PD-L1 + low B2M. The tumor has lost both checkpoint engagement and antigen-presentation capacity. Associated with worse overall survival in melanoma.

The survival association is melanoma-dominant (SKCM log-rank p=0.0011, FDR q=0.0132; RMST +24.2 months at 10 years). The signal survives immune-fraction residualization, purity stratification, and threshold sweeps. The negative direction is consistent across 9/9 datasets (6 TCGA + 3 external melanoma immunotherapy cohorts).

---

## Repository Structure

```
paper1-empirical/
|
|-- README.md                          <-- This file
|-- requirements.txt                   <-- Python dependencies
|-- CHECKER_REPORT.md                  <-- 148-statistic verification report
|
|-- data/                              <-- Source data (see data/README.md)
|-- scripts/                           <-- Analysis scripts (see scripts/README.md)
|-- results/                           <-- Generated outputs
|   |-- csv/                           <-- Statistical results (see results/csv/README.md)
|   |-- figures/                       <-- Figures (see results/figures/README.md)
|-- supplementary_upload/              <-- Journal submission package (see supplementary_upload/README.md)
|-- tables/                            <-- LaTeX table sources (see tables/README.md)
```

Each subdirectory contains its own `README.md` with file-by-file descriptions.

---

## Reproducing All Results

### Step 1: Install Dependencies

```bash
pip install -r requirements.txt
```

**Requirements:** Python 3.10+, numpy, pandas, scipy, matplotlib, seaborn, lifelines, statsmodels, GEOparse.

### Step 2: Run All Analyses

```bash
cd scripts
python run_all.py
```

This executes all 14 analysis scripts in sequence (~15-20 minutes total). All outputs are written to `results/csv/` (statistical tables) and `results/figures/` (figures).

### Step 3: Run Individual Scripts (Optional)

Each script is self-contained and can be run independently from the `scripts/` directory:

```bash
cd scripts

# --- Primary Analyses ---
python tcga_multicancer.py                   # Cross-cancer Spearman correlations (102 tests, 6 cohorts)
python tcga_tumor_normal.py                  # Tumor vs matched normal comparison
python tcga_immune_subtype.py                # AM/DC subtype classification and statistics
python tcga_survival.py                      # Kaplan-Meier survival (12 pre-specified tests)
python tcga_stage_analysis.py                # Stage-stratified analysis + global FDR

# --- Sensitivity Analyses ---
python robustness_check.py                   # Multivariable Cox PH models (8 tests)
python sensitivity_rmst.py                   # RMST at 4 time horizons (bootstrap n=1000)
python sensitivity_threshold_sweep.py        # AM/DC threshold sweep (25th-75th percentile)
python sensitivity_immune_residualization.py  # Immune-fraction residualization + purity tertiles
python composite_observability_index.py       # Independent observability index

# --- External Validation ---
python run_external_validation.py            # GSE91061 nivolumab melanoma (n=49)
python external_validation_gse78220.py       # GSE78220 anti-PD-1 melanoma (n=27)

# --- Convergent Validity ---
python convergent_validity_ncv.py            # nCV convergent validity (Wu & Hogenesch method)
python thorsson_benchmarking.py              # Thorsson immune subtype benchmarking + CCD
```

### Step 4: Verify

Compare generated CSVs in `results/csv/` against the statistics in the manuscript. The [`CHECKER_REPORT.md`](CHECKER_REPORT.md) documents the full stat-by-stat comparison.

---

## Data Provenance

All source data is committed to this repository. No external access is required for reproduction.

| File | Description | Original Source |
|------|-------------|-----------------|
| `data/tcga_expanded_tpm.csv` | Gene expression: 3,920 samples x 25 columns (21 target genes + metadata) | ISB-CGC BigQuery, TCGA RNA-seq hg38 GDC release r35 |
| `data/tcga_clinical.csv` | Clinical and survival data: 3,646 patients x 12 columns | ISB-CGC BigQuery |
| `data/tcga_purity_immune_covariates.csv` | Immune landscape covariates from Thorsson et al. 2018 | Public supplement |
| `scripts/geo_cache/gse91061_*` | GSE91061 nivolumab melanoma (Riaz et al. 2017) | GEO public |
| `scripts/geo_cache/GSE78220_*` | GSE78220 anti-PD-1 melanoma (Hugo et al. 2016) | GEO public |
| `scripts/geo_cache/GSE115821_*` | GSE115821 melanoma (Auslander et al. 2018) | GEO public |

---

## Manuscript

The manuscript is available as a preprint on [Research Square](https://doi.org/10.21203/rs.3.rs-6176498/v1).

---

## Figure-to-Manuscript Mapping

### Main Text Figures

| Figure | File | Section |
|--------|------|---------|
| 1 | `multicancer_circadian_cv_forest.png` | Introduction: core finding across 6 cancer types |
| 2 | `multicancer_correlation_heatmap.png` | Results 3.1: cross-cancer correlation structure |
| 3 | `survival_boundary_failure.png` | Introduction: Kaplan-Meier by subtype |
| 4 | `tumor_normal_circadian_cv_paired.png` | Results 3.4: tumor vs matched normal CV |
| 5 | `tumor_normal_circadian_cv_summary.png` | Results 3.4: cohort-level summary |
| 6 | `immune_subtype_prevalence_stacked.png` | Results 3.5: subtype composition |
| 7 | `immune_subtype_circadian_cv.png` | Results 3.5: CV by boundary-failure subtype |
| 8 | `immune_subtype_boundary_scatter.png` | Results 3.5: checkpoint-clock boundary space |
| 9 | `survival_test_qvalue_summary.png` | Results 3.6: FDR-adjusted survival evidence |
| 10 | `robustness_primary_hr_forest.png` | Results 3.7: hazard-ratio forest plot |
| 11 | `rmst_am_vs_dc.png` | Results 3.7: RMST difference at 4 time horizons |
| 12 | `stage_circadian_cv_boxplot.png` | Results 3.8: CV by pathologic stage |
| 13 | `immune_residualization_scatter.png` | Results 3.9: CV vs PD-L1 after immune residualization |
| 14 | `purity_stratified_correlations.png` | Results 3.9: purity-stratified correlations |
| 15 | `external_validation_geo.png` | Results 3.11: GSE91061 external replication |
| 16 | `thorsson_benchmarking.png` | Results 3.12: Thorsson immune subtype + CCD |

### Supplementary Figures

| Figure | File | Content |
|--------|------|---------|
| S1 | `hypothesis1_TCGA_SKCM.png` | Discovery cohort detail (melanoma) |
| S2 | `hypothesis1_TCGA_LUAD.png` | Discovery cohort detail (lung adenocarcinoma) |
| S3 | `survival_circadian_quartile.png` | CV quartile survival analysis (null result) |
| S4 | `threshold_sensitivity_sweep.png` | AM/DC threshold sensitivity (25th-75th percentile) |
| S5 | `convergent_validity_ncv.png` | nCV convergent validity across 6 cohorts |
| S6 | `observability_index_survival.png` | Observability index survival analysis |
| S7 | `control_budget_combined.png` | Proliferation vs circadian coherence |
| S8 | `external_validation_gse78220.png` | GSE78220 melanoma replication |

---

## Key Results at a Glance

| Result | Value |
|--------|-------|
| Primary coupling (CD274 vs Circadian CV) | Negative in 6/6 TCGA cohorts, all FDR q < 0.05 |
| Strongest effect | SKCM: rho = -0.381, p = 9.6e-18 |
| Melanoma survival (AM vs DC) | Log-rank p = 0.0011, FDR q = 0.0132 |
| RMST advantage (AM over DC, 10yr) | +24.2 months (95% CI: 10.8-36.8) |
| PD-L1 x clock interaction (SKCM) | Cox HR = 0.217, FDR q = 0.045 |
| Immune residualization | Signal persists in 4/4 cohorts |
| Purity stratification | 9/9 tertile-level tests significant (SKCM, LUAD, BRCA) |
| nCV convergent validity | rho = 0.30-0.60, all p < 0.003 (6/6 cohorts) |
| CCD convergent validity | PD-L1 coupling replicates in 4/6 cohorts |
| Thorsson within-subtype | CV-PD-L1 persists in 4/5 immune subtypes |
| External replication direction | Negative in 9/9 datasets |
| Pre-specified evidence audit | 4 FDR-significant, 3 nominal, 25 non-significant (32 total) |

---

## Verification

The [`CHECKER_REPORT.md`](CHECKER_REPORT.md) documents a comprehensive automated verification:

- **148 statistics** cross-checked against output CSVs
- **0 discrepancies**
- **12/12 scripts** re-executed without error
- **11/11 CSVs** audited
- **Split compliance** verified (no theoretical content in empirical paper)

---

## Citation

> Rahnama, M. (2026). Circadian Clock-Gene Structure Couples with PD-L1 Across Six Cancer Types and Defines Survival-Associated Subtypes in Melanoma. *Research Square* (Preprint). https://doi.org/10.21203/rs.3.rs-6176498/v1

---

## Contact

Moses Rahnama
Mina Analytics
moses@minaanalytics.com
