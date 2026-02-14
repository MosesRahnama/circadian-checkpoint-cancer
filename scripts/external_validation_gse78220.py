#!/usr/bin/env python3
"""
External Validation #2: GSE78220 (Hugo et al. 2016, Cell)
=========================================================
Melanoma anti-PD-1 immunotherapy cohort, n=28 pre-treatment.
Tests CV-PD-L1 coupling in a second independent melanoma cohort.

Hugo et al. (2016). Genomic and Transcriptomic Features of Response
to Anti-PD-1 Therapy in Metastatic Melanoma. Cell, 165(1), 35-44.
"""
import os
import sys
import warnings
import urllib.request
import gzip
import io

sys.path.insert(0, os.path.dirname(__file__))
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_PAPER_ROOT = os.path.dirname(SCRIPT_DIR)
CACHE_DIR = os.path.join(SCRIPT_DIR, 'geo_cache')
FIG_DIR = os.path.join(_PAPER_ROOT, 'results', 'figures')
CSV_DIR = os.path.join(_PAPER_ROOT, 'results', 'csv')
os.makedirs(CACHE_DIR, exist_ok=True)
os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(CSV_DIR, exist_ok=True)

CLOCK_GENES = ['ARNTL', 'CLOCK', 'PER1', 'PER2', 'CRY1', 'CRY2']
GENE_ALIASES = {
    'BMAL1': 'ARNTL', 'ARNTL1': 'ARNTL',
}

# GSE78220 supplementary file URL
GSE78220_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE78nnn/GSE78220/suppl/GSE78220_PatientFPKM.xlsx"
GSE78220_CACHE = os.path.join(CACHE_DIR, "GSE78220_PatientFPKM.xlsx")


def download_gse78220():
    """Download the GSE78220 expression matrix if not cached."""
    if os.path.exists(GSE78220_CACHE):
        print(f"Using cached: {GSE78220_CACHE}")
        return GSE78220_CACHE

    print(f"Downloading GSE78220 from NCBI GEO...")
    try:
        urllib.request.urlretrieve(GSE78220_URL, GSE78220_CACHE)
        print(f"Saved: {GSE78220_CACHE}")
    except Exception as e:
        # Try alternative: CSV format
        alt_url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE78nnn/GSE78220/suppl/GSE78220_PatientFPKM.txt.gz"
        alt_cache = os.path.join(CACHE_DIR, "GSE78220_PatientFPKM.txt.gz")
        try:
            urllib.request.urlretrieve(alt_url, alt_cache)
            print(f"Saved (txt.gz): {alt_cache}")
            return alt_cache
        except Exception as e2:
            print(f"Download failed: {e2}")
            return None

    return GSE78220_CACHE


def load_gse78220(filepath):
    """Load and process GSE78220 expression data."""
    if filepath is None:
        return None

    print(f"Loading: {filepath}")
    if filepath.endswith('.xlsx'):
        df = pd.read_excel(filepath)
    elif filepath.endswith('.txt.gz'):
        df = pd.read_csv(filepath, sep='\t', compression='gzip')
    elif filepath.endswith('.txt'):
        df = pd.read_csv(filepath, sep='\t')
    else:
        df = pd.read_csv(filepath)

    print(f"Raw shape: {df.shape}")
    print(f"Columns: {list(df.columns[:5])}...")

    # The first column should be gene names/IDs
    gene_col = df.columns[0]
    df = df.set_index(gene_col)

    # Check which target genes are present
    all_genes = set(df.index)
    found = {}
    for g in CLOCK_GENES + ['CD274']:
        if g in all_genes:
            found[g] = g
        else:
            # Check aliases
            for alias, canonical in GENE_ALIASES.items():
                if alias in all_genes and canonical == g:
                    found[g] = alias
                    break

    print(f"\nGenes found: {list(found.keys())}")
    missing = [g for g in CLOCK_GENES + ['CD274'] if g not in found]
    if missing:
        print(f"Missing: {missing}")

    if 'CD274' not in found or len([g for g in CLOCK_GENES if g in found]) < 4:
        print("Insufficient genes for analysis")
        return None

    # Extract expression for target genes (baseline only, exclude on-treatment)
    sample_cols = [c for c in df.columns if c != gene_col and 'OnTx' not in c and 'on.tx' not in c.lower()]
    print(f"Baseline samples: {len(sample_cols)} (excluded on-treatment)")

    result = pd.DataFrame(index=sample_cols)
    for gene, idx_name in found.items():
        vals = df.loc[idx_name, sample_cols]
        if isinstance(vals, pd.DataFrame):
            vals = vals.iloc[0]
        result[gene] = vals.values.astype(float)

    return result


def main():
    print("=" * 60)
    print("  External Validation #2: GSE78220 (Hugo et al. 2016)")
    print("  Melanoma anti-PD-1, pre-treatment RNA-seq")
    print("=" * 60)

    filepath = download_gse78220()
    expr = load_gse78220(filepath)

    if expr is None:
        print("\nCould not load GSE78220. Skipping.")
        return

    # Log2(FPKM+1) transform
    for col in expr.columns:
        expr[col] = np.log2(expr[col] + 1)

    # Compute circadian CV
    clock_present = [g for g in CLOCK_GENES if g in expr.columns]
    print(f"\nClock genes available: {clock_present} ({len(clock_present)}/6)")
    expr['Circadian_CV'] = expr[clock_present].std(axis=1) / expr[clock_present].mean(axis=1)

    n = len(expr)
    print(f"Samples with valid CV: {n}")

    # CV-PD-L1 correlation
    if 'CD274' in expr.columns:
        rho, p = stats.spearmanr(expr['Circadian_CV'], expr['CD274'])
        sig = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns'
        direction = "negative (replicates TCGA)" if rho < 0 else "positive (does NOT replicate)"

        print(f"\n{'=' * 60}")
        print(f"  CV-PD-L1 Spearman correlation:")
        print(f"  rho = {rho:.4f}, p = {p:.4f} {sig}")
        print(f"  Direction: {direction}")
        print(f"  n = {n}")
        print(f"{'=' * 60}")

        # Save results
        results = pd.DataFrame([{
            'dataset': 'GSE78220',
            'reference': 'Hugo et al. 2016 Cell',
            'cancer_type': 'Melanoma (anti-PD-1)',
            'n': n,
            'cv_pdl1_rho': round(rho, 4),
            'cv_pdl1_p': round(p, 4),
            'direction': 'negative' if rho < 0 else 'positive',
            'replicates_tcga': rho < 0,
            'clock_genes_used': len(clock_present),
        }])

        csv_path = os.path.join(CSV_DIR, 'external_validation_gse78220.csv')
        results.to_csv(csv_path, index=False)
        print(f"\nSaved: {csv_path}")

        # Simple scatter plot
        fig, ax = plt.subplots(1, 1, figsize=(6, 5))
        ax.scatter(expr['CD274'], expr['Circadian_CV'], alpha=0.6, s=40, color='#e74c3c')
        ax.set_xlabel('CD274 (PD-L1) log2(FPKM+1)')
        ax.set_ylabel('Circadian CV')
        ax.set_title(f'GSE78220: CV vs PD-L1 (n={n}, rho={rho:.3f}, p={p:.3f})')

        # Add trend line
        z = np.polyfit(expr['CD274'], expr['Circadian_CV'], 1)
        x_line = np.linspace(expr['CD274'].min(), expr['CD274'].max(), 100)
        ax.plot(x_line, np.polyval(z, x_line), '--', color='gray', alpha=0.7)

        fig_path = os.path.join(FIG_DIR, 'external_validation_gse78220.png')
        plt.savefig(fig_path, dpi=300, bbox_inches='tight')
        print(f"Saved: {fig_path}")
    else:
        print("CD274 not found in dataset")


if __name__ == '__main__':
    main()
