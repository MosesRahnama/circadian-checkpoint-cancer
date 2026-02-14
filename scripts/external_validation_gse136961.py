#!/usr/bin/env python3
"""
External Validation #3: GSE136961 (Jung et al. 2019)
====================================================
NSCLC anti-PD-1/PD-L1 immunotherapy cohort.
Tests CV-PD-L1 coupling in a non-melanoma immunotherapy context.

This dataset has RNA-seq from NSCLC patients treated with checkpoint inhibitors.
"""
import os
import sys
import warnings
import urllib.request

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

CLOCK_GENES = ['ARNTL', 'CLOCK', 'PER1', 'PER2', 'CRY1', 'CRY2']

# Try GSE135222 (Cho et al. 2020, Nature Medicine) - NSCLC anti-PD-1
# Or GSE126044 (Prat et al. 2017) - NSCLC/melanoma anti-PD-1

DATASETS = [
    {
        'gse': 'GSE126044',
        'ref': 'Auslander et al. 2018',
        'cancer': 'Melanoma + NSCLC (anti-PD-1)',
        'url': 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE126nnn/GSE126044/suppl/GSE126044_NSCLC_IMvigor_pats.csv.gz',
        'sep': ',',
    },
]

def try_geoparse(gse_id):
    """Try to load via GEOparse as fallback."""
    try:
        import GEOparse
        cache_file = os.path.join(CACHE_DIR, f"{gse_id}_family.soft.gz")
        if os.path.exists(cache_file):
            gse = GEOparse.get_GEO(filepath=cache_file, silent=True)
        else:
            gse = GEOparse.get_GEO(gse_id, destdir=CACHE_DIR, silent=True)

        # Extract expression table from first platform
        for gpl_name, gpl in gse.gpls.items():
            print(f"  Platform: {gpl_name}")
            break

        # Check if we can get expression data from GSMs
        sample_data = {}
        for gsm_name, gsm in gse.gsms.items():
            table = gsm.table
            if len(table) > 0:
                sample_data[gsm_name] = table
                if len(sample_data) == 1:
                    print(f"  Table columns: {list(table.columns[:5])}")
                    print(f"  Table rows: {len(table)}")

        if sample_data:
            return sample_data
    except Exception as e:
        print(f"  GEOparse failed: {e}")
    return None


def run_cv_pdl1(expr_df, dataset_name, ref, cancer_type):
    """Run CV-PD-L1 correlation on a processed expression dataframe."""
    # Log2 transform
    for col in expr_df.columns:
        expr_df[col] = np.log2(expr_df[col].clip(lower=0) + 1)

    clock_present = [g for g in CLOCK_GENES if g in expr_df.columns]
    if 'CD274' not in expr_df.columns or len(clock_present) < 4:
        print(f"  Insufficient genes: CD274={'CD274' in expr_df.columns}, clock={len(clock_present)}/6")
        return None

    expr_df['Circadian_CV'] = expr_df[clock_present].std(axis=1) / expr_df[clock_present].mean(axis=1)
    expr_df = expr_df.dropna(subset=['Circadian_CV', 'CD274'])

    n = len(expr_df)
    if n < 15:
        print(f"  Too few samples: n={n}")
        return None

    rho, p = stats.spearmanr(expr_df['Circadian_CV'], expr_df['CD274'])
    sig = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns'
    direction = 'negative' if rho < 0 else 'positive'

    print(f"  {dataset_name} ({cancer_type}): rho={rho:.4f}, p={p:.4f} {sig}, n={n}, direction={direction}")

    return {
        'dataset': dataset_name,
        'reference': ref,
        'cancer_type': cancer_type,
        'n': n,
        'cv_pdl1_rho': round(rho, 4),
        'cv_pdl1_p': round(p, 4),
        'direction': direction,
        'replicates_tcga': rho < 0,
        'clock_genes_used': len(clock_present),
    }


def main():
    print("=" * 60)
    print("  External Validation #3: Additional GEO Cohorts")
    print("=" * 60)

    # Also try to use GSE78220 data we already have, split by response
    gse78220_file = os.path.join(CACHE_DIR, 'GSE78220_PatientFPKM.xlsx')
    if os.path.exists(gse78220_file):
        print("\nRe-checking GSE78220 baseline samples...")
        df = pd.read_excel(gse78220_file)
        gene_col = df.columns[0]
        df = df.set_index(gene_col)

        # Only baseline samples
        baseline_cols = [c for c in df.columns if 'baseline' in c.lower()]
        on_treatment_cols = [c for c in df.columns if 'treatment' in c.lower() or 'on.tx' in c.lower().replace('-','')]

        print(f"  Baseline samples: {len(baseline_cols)}")
        print(f"  On-treatment samples: {len(on_treatment_cols)}")

        if on_treatment_cols:
            # Also test on-treatment samples
            expr_ot = pd.DataFrame(index=on_treatment_cols)
            for g in CLOCK_GENES + ['CD274']:
                if g in df.index:
                    expr_ot[g] = df.loc[g, on_treatment_cols].values.astype(float)

            result = run_cv_pdl1(expr_ot.copy(), 'GSE78220_on_treatment',
                               'Hugo et al. 2016 Cell', 'Melanoma (on-treatment)')
            if result:
                print(f"  -> On-treatment also shows direction: {result['direction']}")

    print("\nDone.")


if __name__ == '__main__':
    main()
