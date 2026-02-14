#!/usr/bin/env python3
"""
Master runner for Paper 1: reproduces all analyses from source data.

Run from the scripts/ directory:
    cd paper1-empirical/scripts
    python run_all.py

All source data is in ../data/ (committed, no internet needed).
All outputs go to ../results/csv/ and ../results/figures/.
"""
import subprocess
import sys
import os
import time

SCRIPTS = [
    ("tcga_multicancer.py",                "Cross-cancer correlations (102 tests, 6 cohorts)"),
    ("tcga_tumor_normal.py",               "Tumor vs matched normal comparison"),
    ("tcga_immune_subtype.py",             "Active Masking / Decoherence classification"),
    ("tcga_survival.py",                   "Kaplan-Meier survival (12 pre-specified tests)"),
    ("tcga_stage_analysis.py",             "Stage-stratified analysis + global FDR (245 tests)"),
    ("robustness_check.py",                "Multivariable Cox models (8 tests, SKCM+LUAD)"),
    ("sensitivity_rmst.py",                "RMST at 4 time horizons (bootstrap n=1000)"),
    ("sensitivity_threshold_sweep.py",     "AM/DC threshold sweep (25th-75th percentile)"),
    ("sensitivity_immune_residualization.py", "Immune-fraction residualization + purity tertiles"),
    ("composite_observability_index.py",   "Independent observability index"),
    ("run_external_validation.py",         "GSE91061 nivolumab melanoma replication"),
    ("convergent_validity_ncv.py",         "nCV convergent validity (Wu & Hogenesch method)"),
    ("thorsson_benchmarking.py",           "Thorsson immune subtype benchmarking + CCD"),
]

def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    paper_root = os.path.dirname(script_dir)

    # Ensure output directories exist
    os.makedirs(os.path.join(paper_root, "results", "csv"), exist_ok=True)
    os.makedirs(os.path.join(paper_root, "results", "figures"), exist_ok=True)

    print("=" * 70)
    print("  Paper 1: Full Reproduction Pipeline")
    print("=" * 70)
    print(f"  Source data: {os.path.join(paper_root, 'data')}")
    print(f"  Output CSVs: {os.path.join(paper_root, 'results', 'csv')}")
    print(f"  Output figs: {os.path.join(paper_root, 'results', 'figures')}")
    print("=" * 70)

    passed = 0
    failed = 0
    results = []

    for script, description in SCRIPTS:
        script_path = os.path.join(script_dir, script)
        print(f"\n{'-' * 70}")
        print(f"  [{passed + failed + 1}/{len(SCRIPTS)}] {description}")
        print(f"  Script: {script}")
        print(f"{'-' * 70}")

        t0 = time.time()
        try:
            result = subprocess.run(
                [sys.executable, script_path],
                cwd=script_dir,
                capture_output=True,
                text=True,
                timeout=600,
            )
            elapsed = time.time() - t0

            if result.returncode == 0:
                passed += 1
                status = "PASS"
                print(f"  Status: PASS ({elapsed:.1f}s)")
            else:
                failed += 1
                status = "FAIL"
                print(f"  Status: FAIL (exit code {result.returncode})")
                print(f"  stderr: {result.stderr[-500:]}")
        except subprocess.TimeoutExpired:
            failed += 1
            status = "TIMEOUT"
            elapsed = 600
            print(f"  Status: TIMEOUT (>600s)")
        except Exception as e:
            failed += 1
            status = "ERROR"
            elapsed = time.time() - t0
            print(f"  Status: ERROR ({e})")

        results.append((script, status, f"{elapsed:.1f}s"))

    print(f"\n{'=' * 70}")
    print(f"  SUMMARY: {passed}/{len(SCRIPTS)} passed, {failed} failed")
    print(f"{'=' * 70}")
    for script, status, elapsed in results:
        marker = "PASS" if status == "PASS" else "FAIL"
        print(f"  [{marker}] {script} ({elapsed})")

    if failed > 0:
        print(f"\n  WARNING: {failed} script(s) failed.")
        sys.exit(1)
    else:
        print(f"\n  All {passed} scripts completed successfully.")
        print(f"  Results are in: {os.path.join(paper_root, 'results')}")

if __name__ == "__main__":
    main()
