# LaTeX Tables

This directory contains LaTeX table source files used in the manuscript.

---

## Files

### `evidence_table_prespecified.tex`

**Description:** Pre-specified endpoint evidence table in LaTeX format. Lists all 32 pre-specified endpoint families used for inferential claims in the manuscript:
- 12 survival tests (6 AM-vs-DC + 6 Q1-vs-Q4, BH-corrected across all survival tests)
- 8 primary robustness Cox tests (BH-corrected across primary robustness tests)
- 12 circadian-stage tests (q-values from global BH correction)

Each endpoint is categorized as:
- **FDR-significant** (q < 0.05)
- **Nominal** (p < 0.05 with q >= 0.05)
- **Non-significant**

**Source data:** Generated from `results/csv/evidence_ledger_prespecified.csv`
**Manuscript reference:** Evidence Summary Table (Table 2)
