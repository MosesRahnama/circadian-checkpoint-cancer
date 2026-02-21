# Source Data

This directory contains all source data required to reproduce the analyses in the manuscript. All files are committed to the repository. No external accounts, API keys, or internet access is required.

---

## Files

### `tcga_expanded_tpm.csv`

**Description:** Bulk RNA-seq gene expression data from six TCGA cancer cohorts.

- **Rows:** 3,920 samples (3,611 tumor + 309 matched normal)
- **Columns:** 25 (21 target genes + sample metadata)
- **Expression values:** Transcripts per million (TPM), log2(TPM+1) transformed in analysis scripts
- **Target genes (21):**
  - Immune checkpoint: CD274 (PD-L1), PDCD1LG2, PDCD1
  - MHC class I: HLA-A, HLA-B, HLA-C, B2M
  - Gap junction: GJA1 (Cx43), GJB2 (Cx26), GJA5 (Cx40), GJB6 (Cx30)
  - Circadian clock: ARNTL (BMAL1), CLOCK, PER1, PER2, CRY1, CRY2
  - Differentiation markers: CDH1, VIM, MYC, TP53
- **Cohorts included:**
  - TCGA-SKCM (melanoma): 472 tumor, 1 normal
  - TCGA-LUAD (lung adenocarcinoma): 530 tumor, 59 normal
  - TCGA-BRCA (breast): 1,113 tumor, 113 normal
  - TCGA-COAD (colon): 473 tumor, 41 normal
  - TCGA-HNSC (head/neck squamous cell): 522 tumor, 44 normal
  - TCGA-LUSC (lung squamous cell): 501 tumor, 51 normal
- **Original source:** ISB-CGC BigQuery table `isb-cgc-bq.TCGA.RNAseq_hg38_gdc_r35` (GDC release r35, version-locked)

---

### `tcga_clinical.csv`

**Description:** Clinical and survival data for TCGA patients across six cohorts.

- **Rows:** 3,646 patients
- **Columns:** 12
- **Key fields:** patient barcode, cancer type, overall survival time, vital status, age at diagnosis, sex, pathologic tumor stage
- **Used for:** Kaplan-Meier survival analysis, Cox proportional hazards models, stage-stratified analysis
- **Original source:** ISB-CGC BigQuery (TCGA clinical tables)

---

### `tcga_purity_immune_covariates.csv`

**Description:** Immune landscape and tumor purity covariates from Thorsson et al. (2018), "The Immune Landscape of Human Tumors," *Immunity* 48(4):812-830.

- **Key fields:** tumor purity, leukocyte fraction, lymphocyte infiltration, IFN-gamma response, TGF-beta response, stromal fraction, myeloid infiltration, immune subtype (C1-C6)
- **Used for:** Immune-fraction residualization, purity stratification, Thorsson immune subtype benchmarking, multivariable Cox models
- **Original source:** Thorsson et al. 2018 public supplement

---

## Data Provenance Note

All three files were extracted once from their original public sources and committed to this repository. The extraction script (`scripts/tcga_extract_expanded.py`, referenced in the manuscript) documents the BigQuery SQL used. Reviewers can verify provenance by comparing against the original ISB-CGC BigQuery tables or the Thorsson et al. supplement.
