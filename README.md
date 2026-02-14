# Paper 1: Circadian-Checkpoint Coupling Defines Boundary-Failure Subtypes with Survival Associations Across Six Cancer Types

**Author:** Moses Rahnama, Mina Analytics (moses@minaanalytics.com)
**Repository:** https://github.com/MosesRahnama/cancer-paper-repository
**Date:** February 2026

---

## Manuscript

| File | Location |
|------|----------|
| LaTeX source | `paper/paper1-empirical/paper1_circadian_checkpoint.tex` |
| Compiled PDF | `paper/paper1-empirical/paper1_circadian_checkpoint.pdf` |
| Original mega-manuscript | `paper/Cancer_As_Boundary_Logic_Failure.tex` |
| Split instructions | `paper/manuscript_split_instructions.md` |

---

## Source Data (inputs to all analyses)

All committed to the repository. No GCP account or internet access required to reproduce.

| File | Description | Location |
|------|-------------|----------|
| `tcga_expanded_tpm.csv` | 3,920 samples x 25 columns (21 target genes + metadata) | `experiments/tcga/` |
| `tcga_clinical.csv` | 3,646 patients x 12 columns (survival, staging, demographics) | `experiments/tcga/` |
| `tcga_purity_immune_covariates.csv` | Thorsson et al. immune landscape covariates | `experiments/tcga/` |
| GSE91061 cache | GEO nivolumab melanoma (downloaded by script on first run) | `experiments/tcga/geo_cache/` |

**Data provenance:** ISB-CGC BigQuery `isb-cgc-bq.TCGA.RNAseq_hg38_gdc_r35` (release-locked to r35). Extraction script: `tcga_extract_expanded.py`.

---

## Analysis Scripts

All scripts are in `experiments/tcga/`. Run from that directory.

### Primary Analyses (Paper 1 results)

| Script | What it produces | Output CSV | Output figures |
|--------|-----------------|------------|----------------|
| `tcga_multicancer.py` | Cross-cancer Spearman correlations (102 tests, 6 cohorts) | `tcga_multicancer_correlations.csv` | `multicancer_correlation_heatmap.png`, `multicancer_circadian_cv_forest.png` |
| `tcga_tumor_normal.py` | Tumor vs matched normal (Wilcoxon, Mann-Whitney) | `tumor_normal_comparison.csv` | `tumor_normal_circadian_cv_paired.png`, `tumor_normal_circadian_cv_summary.png` |
| `tcga_immune_subtype.py` | AM/DC classification, Kruskal-Wallis, Mann-Whitney | `immune_subtype_comparison.csv` | `immune_subtype_circadian_cv.png`, `immune_subtype_boundary_scatter.png`, `immune_subtype_prevalence_stacked.png` |
| `tcga_survival.py` | Kaplan-Meier, log-rank (12 pre-specified tests) | `survival_logrank_results.csv` | `survival_boundary_failure.png`, `survival_circadian_quartile.png` |
| `tcga_stage_analysis.py` | Stage-stratified analysis + global FDR (244 tests) | `stage_analysis_results.csv`, `master_fdr_results.csv` | `stage_circadian_cv_boxplot.png` |

### Sensitivity Analyses

| Script | What it produces | Output CSV | Output figures |
|--------|-----------------|------------|----------------|
| `robustness_check.py` | Multivariable Cox PH (8 models, SKCM+LUAD) | `robustness_primary_tests.csv`, `robustness_cox_coefficients.csv`, `robustness_partial_correlations.csv`, `robustness_ph_diagnostics.csv` | `robustness_primary_hr_forest.png` |
| `sensitivity_rmst.py` | RMST at 4 time horizons (bootstrap n=1000, seed=42) | `rmst_results.csv` | `rmst_am_vs_dc.png` |
| `sensitivity_threshold_sweep.py` | AM/DC threshold 25th-75th percentile (22 tests) | `threshold_sensitivity_results.csv` | `threshold_sensitivity_sweep.png` |
| `sensitivity_immune_residualization.py` | OLS residualization + purity tertiles | `immune_residualization_results.csv` | `immune_residualization_scatter.png`, `purity_stratified_correlations.png` |
| `composite_observability_index.py` | Independent observability index (orthogonal features) | `observability_index_results.csv` | `observability_index_survival.png` |

### External Validation

| Script | What it produces | Output CSV | Output figures |
|--------|-----------------|------------|----------------|
| `run_external_validation.py` | GSE91061 nivolumab melanoma (n=49) | `external_validation_results.csv` | `external_validation_geo.png` |
| `external_validation_geo.py` | GEO data download + processing helper | (used by above) | |

### Convergent Validity (added post-review)

| Script | What it produces | Output CSV | Output figures |
|--------|-----------------|------------|----------------|
| `convergent_validity_ncv.py` | Wu & Hogenesch nCV comparison (6 cohorts) | `convergent_validity_results.csv` | `convergent_validity_ncv.png` |

### Utility Scripts

| Script | Purpose |
|--------|---------|
| `tcga_config.py` | Shared configuration and utility functions |
| `tcga_extract_expanded.py` | BigQuery data extraction (optional; data already committed) |
| `generate_manuscript_figures.py` | Manuscript figure generation |
| `sync_manuscript_artifacts.py` | Deterministic artifact sync + provenance |
| `ci_validate_committed_outputs.py` | CI/CD validation of committed outputs |

---

## Figures

### Main paper figures (19 files)

Located in `paper/paper1-empirical/figures/` (copies) and `results/` (originals).

| Figure | File | Section |
|--------|------|---------|
| 1 (Intro) | `multicancer_circadian_cv_forest.png` | Introduction |
| 2 | `multicancer_correlation_heatmap.png` | Results 3.1 |
| 3 | `tumor_normal_circadian_cv_paired.png` | Results 3.4 |
| 4 | `tumor_normal_circadian_cv_summary.png` | Results 3.4 |
| 5 | `immune_subtype_prevalence_stacked.png` | Results 3.5 |
| 6 | `immune_subtype_circadian_cv.png` | Results 3.5 |
| 7 | `immune_subtype_boundary_scatter.png` | Results 3.5 |
| 8 | `survival_boundary_failure.png` | Results 3.6 |
| 9 | `survival_test_qvalue_summary.png` | Results 3.6 |
| 10 | `robustness_primary_hr_forest.png` | Results 3.7 |
| 11 | `rmst_am_vs_dc.png` | Results 3.7 |
| 12 | `stage_circadian_cv_boxplot.png` | Results 3.8 |
| 13 | `immune_residualization_scatter.png` | Results 3.9 |
| 14 | `purity_stratified_correlations.png` | Results 3.9 |
| 15 | `external_validation_geo.png` | Results 3.11 |

### Supplementary figures (8 files)

Located in `results/csv/` and `results/figures/`.

| Figure | File | Content |
|--------|------|---------|
| S1 | `hypothesis1_TCGA_SKCM.png` | Discovery cohort detail (SKCM) |
| S2 | `hypothesis1_TCGA_LUAD.png` | Discovery cohort detail (LUAD) |
| S3 | `survival_circadian_quartile.png` | CV quartile survival (null) |
| S4 | `threshold_sensitivity_sweep.png` | Threshold sensitivity |
| S5 | `stage_circadian_spearman_summary.png` | Stage trend summary |
| S6 | `observability_index_survival.png` | Observability index survival |
| S7 | `control_budget_combined.png` | Proliferation vs coherence |
| S8 | `convergent_validity_ncv.png` | nCV convergent validity |

---

## Supplementary Data Files (16 files)

Located in `results/csv/` and `results/figures/`.

| File | Description | Rows |
|------|-------------|------|
| `tcga_multicancer_correlations.csv` | All Spearman correlations (6 cohorts x 17 comparisons) | 102 |
| `survival_logrank_results.csv` | All 12 survival tests with FDR | 12 |
| `robustness_primary_tests.csv` | Cox model results + PH diagnostics | 8 |
| `robustness_cox_coefficients.csv` | Full coefficient tables | ~200 |
| `robustness_partial_correlations.csv` | Partial correlations after adjustment | ~10 |
| `rmst_results.csv` | RMST at 4 horizons with bootstrap CIs | 8 |
| `immune_residualization_results.csv` | Residualized + purity-stratified results | ~16 |
| `immune_subtype_comparison.csv` | AM/DC prevalence and CV stats | 6 |
| `observability_index_results.csv` | Independent index results | 6 |
| `external_validation_results.csv` | GSE91061 replication | 1 |
| `tumor_normal_comparison.csv` | Tumor vs normal comparisons | ~35 |
| `threshold_sensitivity_results.csv` | Threshold sweep (22 tests) | 22 |
| `stage_analysis_results.csv` | Stage-stratified analysis | ~48 |
| `master_fdr_results.csv` | Global FDR across 244 tests | 244 |
| `evidence_ledger_prespecified.csv` | 32-endpoint evidence audit | 32 |
| `convergent_validity_results.csv` | nCV convergent validity | ~24 |

---

## Verification Reports

| File | Description | Location |
|------|-------------|----------|
| `CHECKER_REPORT.md` | 80+ statistics cross-checked against CSVs | `paper/paper1-empirical/` |
| `OBJECTIVE_ANALYSIS.md` | Independent novelty/literature analysis | `paper/paper1-empirical/` |

---

## Quick Reproduction

```bash
cd experiments/tcga

# Core analyses (~10 min total)
python tcga_multicancer.py
python tcga_tumor_normal.py
python tcga_immune_subtype.py
python tcga_survival.py
python tcga_stage_analysis.py

# Sensitivity analyses (~10 min total)
python robustness_check.py
python sensitivity_rmst.py
python sensitivity_threshold_sweep.py
python sensitivity_immune_residualization.py
python composite_observability_index.py

# External validation
python run_external_validation.py

# Convergent validity
python convergent_validity_ncv.py
```

**Requirements:** Python 3.10+, pandas, numpy, scipy, matplotlib, lifelines, statsmodels, seaborn, GEOparse.

```bash
pip install -r requirements-tcga.txt
```

---

## References (29 total in paper)

| Category | Count | Key citations |
|----------|-------|---------------|
| Circadian disruption | 5 | IARC 2019, Schernhammer 2006, Finnish 2023, Mello/Masri/Lamia 2026, Masri 2016 |
| Immune checkpoint | 2 | Lin et al. 2024, Nowicki et al. 2018 |
| Immunoediting | 2 | Dunn et al. 2004, Schreiber et al. 2011 |
| Immune landscape | 1 | Thorsson et al. 2018 |
| Hallmarks | 2 | Hanahan 2011, 2022 |
| Circadian-checkpoint mechanistic | 3 | Liu/RORA 2024, Quan/CLOCK 2025, Fortin et al. 2024 |
| Prior art | 2 | de Assis & Kinker 2018, Wu et al. 2019 |
| Clock metric | 1 | Wu & Hogenesch nCV 2021 |
| TCGA/ISB-CGC | 1 | Reynolds et al. 2017 |
| Gap junctions | 2 | Connexins 2018, Neijssen 2005 |
| Bioelectric | 2 | Levin 2021, 2013 |
| Stemness/entropy | 2 | Teschendorff 2017, Malta 2018 |
| External validation | 1 | Riaz et al. 2017 |
| Self-citations | 2 | Rahnama 2026 framework, Rahnama 2026 therapeutic |

---

## Split Plan

This is Paper 1 of a planned 4-paper split. See `paper/manuscript_split_instructions.md` for full details.

| Paper | Focus | Status |
|-------|-------|--------|
| **1 (this)** | Empirical: TCGA circadian-checkpoint coupling | Written |
| 2 | Theoretical: boundary logic framework | Planned |
| 3 | Simulation: tissue-graph therapeutic operators | Planned |
| 4 | Perspective: forced distinction design principle | Planned |
