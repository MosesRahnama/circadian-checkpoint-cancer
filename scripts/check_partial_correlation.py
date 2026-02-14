
import pandas as pd
import numpy as np
import scipy.stats as stats
import os

# Configuration
DATA_PATH = "data/tcga_expanded_tpm.csv"
CIRCADIAN = ["ARNTL", "CLOCK", "PER1", "PER2", "CRY1", "CRY2"]

def log_transform(series):
    return np.log2(series + 1)

def partial_corr(x, y, covar):
    # Residualize x and y on covar
    slope_x, intercept_x, _, _, _ = stats.linregress(covar, x)
    resid_x = x - (slope_x * covar + intercept_x)
    
    slope_y, intercept_y, _, _, _ = stats.linregress(covar, y)
    resid_y = y - (slope_y * covar + intercept_y)
    
    return stats.spearmanr(resid_x, resid_y)

def main():
    if not os.path.exists(DATA_PATH):
        print(f"File not found: {DATA_PATH}")
        return

    df = pd.read_csv(DATA_PATH)
    print(f"Loaded {len(df)} samples")

    results = []
    
    for proj in ["TCGA-SKCM", "TCGA-LUAD", "TCGA-BRCA", "TCGA-COAD", "TCGA-HNSC", "TCGA-LUSC"]:
        sub = df[df["project_short_name"] == proj].copy()
        
        # Calculate Metrics
        circ_cols = [c for c in CIRCADIAN if c in sub.columns]
        log_vals = sub[circ_cols].apply(log_transform)
        mean_expr = log_vals.mean(axis=1)
        std_expr = log_vals.std(axis=1)
        cv = std_expr / mean_expr
        pdl1 = log_transform(sub["CD274"])
        
        # Filter NaNs
        mask = pdl1.notna() & cv.notna() & mean_expr.notna() & np.isfinite(cv)
        sub = sub[mask]
        pdl1 = pdl1[mask]
        cv = cv[mask]
        mean_expr = mean_expr[mask]
        
        # Original Correlation
        r_orig, p_orig = stats.spearmanr(pdl1, cv)
        
        # Partial Correlation controlling for Mean Expression
        r_partial, p_partial = partial_corr(pdl1, cv, mean_expr)
        
        results.append({
            "project": proj,
            "orig_rho": r_orig,
            "orig_p": p_orig,
            "partial_rho": r_partial,
            "partial_p": p_partial,
            "change": r_partial - r_orig
        })
        
    print("\nOriginal vs Partial Correlation (controlling for Mean Clock Expression):")
    print(pd.DataFrame(results).to_string())

if __name__ == "__main__":
    main()
