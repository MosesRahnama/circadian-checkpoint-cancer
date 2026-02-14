
import pandas as pd
import numpy as np
import scipy.stats as stats
import os

# Configuration
DATA_PATH = "data/tcga_expanded_tpm.csv"
CIRCADIAN = ["ARNTL", "CLOCK", "PER1", "PER2", "CRY1", "CRY2"]

def log_transform(series):
    return np.log2(series + 1)

def compute_metrics(df):
    circ_cols = [c for c in CIRCADIAN if c in df.columns]
    log_vals = df[circ_cols].apply(log_transform)
    
    # Compute Mean and CV
    mean_expr = log_vals.mean(axis=1)
    std_expr = log_vals.std(axis=1)
    cv = std_expr / mean_expr
    
    return mean_expr, cv, log_vals

def main():
    if not os.path.exists(DATA_PATH):
        print(f"File not found: {DATA_PATH}")
        return

    df = pd.read_csv(DATA_PATH)
    print(f"Loaded {len(df)} samples")

    results = []
    
    for proj in df["project_short_name"].unique():
        sub = df[df["project_short_name"] == proj].copy()
        mean_expr, cv, log_vals = compute_metrics(sub)
        
        # Correlation between CV and Mean
        r, p = stats.spearmanr(cv, mean_expr)
        
        results.append({
            "project": proj,
            "rho_cv_mean": r,
            "p_cv_mean": p
        })
        
    print("\nCorrelation between Circadian CV and Mean Clock Expression:")
    print(pd.DataFrame(results).to_string())

if __name__ == "__main__":
    main()
