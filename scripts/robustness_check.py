"""
Robustness and Covariate-Control Analysis for TCGA Boundary-Failure Survival
=========================================================================

Implements a reviewer-resistant robustness layer on top of tcga_survival.py:
1) Per-cohort multivariable Cox model for Active Masking vs Decoherence
2) Stage as categorical covariates (I reference; II/III/IV + missing indicator)
3) Proportional hazards diagnostics via covariate x log(time) interaction screens
4) Continuous PD-L1 x clock interaction model
5) Optional microenvironment covariate merge hooks (purity/immune/stromal/IFN/TGFB)
6) BH correction across primary robustness tests

Dependencies: numpy, pandas, scipy (no lifelines/statsmodels required)
"""

import argparse
import os
import sys
import warnings

import numpy as np
import pandas as pd
from scipy import optimize, stats
from scipy.stats import false_discovery_control

sys.path.insert(0, os.path.dirname(__file__))

from tcga_config import (RESULTS_CSV_DIR,  # noqa: E402
    ALL_PROJECTS,
    DATA_DIR,
    classify_boundary_failure,
    compute_circadian_cv,
    log_transform,
    normalize_stage,
    upload_to_gcs,
)

warnings.filterwarnings("ignore", category=RuntimeWarning)

MIN_ROWS = 40
MIN_EVENTS = 10

PURITY_COL_CANDIDATES = [
    "purity",
    "tumor_purity",
    "cpe",
    "estimate_purity",
    "absolute_purity",
    "consensus_purity",
]

IMMUNE_LYMPH_CANDIDATES = [
    "lymphocyte_infiltration",
    "lymphocyte_score",
    "lymphocyte",
    "immune_lymphoid",
    "lymphoid_score",
]

IMMUNE_MYELOID_CANDIDATES = [
    "myeloid_infiltration",
    "myeloid_score",
    "myeloid",
    "immune_myeloid",
    "macrophage_regulation",
]

IMMUNE_LEUKOCYTE_CANDIDATES = [
    "leukocyte_fraction",
    "leukocyte_score",
    "immune_leukocyte",
]

STROMAL_COL_CANDIDATES = [
    "stromal_fraction",
    "stromal_score",
]

IFNG_COL_CANDIDATES = [
    "ifn_gamma_response",
    "ifng_response",
    "ifn-gamma_response",
]

TGFB_COL_CANDIDATES = [
    "tgfb_response",
    "tgf_beta_response",
    "tgf-beta_response",
]

CASE_ID_CANDIDATES = [
    "case_barcode",
    "bcr_patient_barcode",
    "submitter_id",
    "patient_id",
]

PROJECT_COL_CANDIDATES = [
    "project_short_name",
    "project_id",
    "cohort",
]

OPTIONAL_CONTROL_ORDER = [
    "purity",
    "immune_lymphoid",
    "immune_myeloid",
    "immune_leukocyte",
    "stromal_fraction",
    "ifn_gamma_response",
    "tgfb_response",
]


def deduplicate_one_tumor_sample_per_case(tumors: pd.DataFrame) -> pd.DataFrame:
    """Keep one tumor sample per (project_short_name, case_barcode)."""
    pr = np.full(len(tumors), 4, dtype=int)
    s = tumors["sample_type_name"].fillna("").astype(str)
    pr[s.str.contains("Primary Tumor", case=False)] = 1
    pr[s.str.contains("Recurrent Tumor", case=False)] = 2
    pr[s.str.contains("Metastatic", case=False)] = 3

    out = tumors.copy()
    out["_sample_priority"] = pr
    out = out.sort_values(
        ["project_short_name", "case_barcode", "_sample_priority", "sample_barcode"]
    )
    out = out.drop_duplicates(subset=["project_short_name", "case_barcode"], keep="first")
    return out.drop(columns=["_sample_priority"])


def compute_circadian_cv_orthogonal(df: pd.DataFrame) -> pd.Series:
    """
    Circadian CV sensitivity metric excluding ARNTL and PER1
    to reduce definitional coupling with subtype labels.
    """
    ortho_cols = [c for c in ["CLOCK", "PER2", "CRY1", "CRY2"] if c in df.columns]
    if len(ortho_cols) < 2:
        return pd.Series(np.nan, index=df.index)
    vals = df[ortho_cols].apply(log_transform)
    return vals.std(axis=1) / vals.mean(axis=1)


def standardize_project_name(x: str) -> str:
    if pd.isna(x):
        return x
    s = str(x).strip()
    if s.startswith("TCGA-"):
        return s
    if len(s) == 4:
        return f"TCGA-{s}"
    return s


def find_first_existing_column(df: pd.DataFrame, candidates) -> str | None:
    cols_lower = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in cols_lower:
            return cols_lower[cand.lower()]
    return None


def load_optional_covariates(path: str | None):
    """
    Load optional covariates file and map columns to standardized names:
      - case_barcode (required)
      - project_short_name (optional)
      - purity (optional)
      - immune_lymphoid (optional)
      - immune_myeloid (optional)
      - immune_leukocyte (optional)
      - stromal_fraction (optional)
      - ifn_gamma_response (optional)
      - tgfb_response (optional)
    """
    candidate_paths = []
    if path:
        candidate_paths.append(path)
    candidate_paths.extend([
        os.path.join(DATA_DIR, "tcga_purity_immune_covariates.csv"),
        os.path.join(DATA_DIR, "tcga_purity_covariates.csv"),
        os.path.join(DATA_DIR, "tcga_purity.csv"),
    ])

    chosen = None
    for p in candidate_paths:
        if p and os.path.exists(p):
            chosen = p
            break

    if chosen is None:
        return None, None

    extra = pd.read_csv(chosen)

    case_col = find_first_existing_column(extra, CASE_ID_CANDIDATES)
    if case_col is None:
        print(f"[covariates] Skipping {chosen}: no case identifier column found.")
        return None, chosen

    proj_col = find_first_existing_column(extra, PROJECT_COL_CANDIDATES)
    purity_col = find_first_existing_column(extra, PURITY_COL_CANDIDATES)
    lymph_col = find_first_existing_column(extra, IMMUNE_LYMPH_CANDIDATES)
    myeloid_col = find_first_existing_column(extra, IMMUNE_MYELOID_CANDIDATES)
    leuk_col = find_first_existing_column(extra, IMMUNE_LEUKOCYTE_CANDIDATES)
    stromal_col = find_first_existing_column(extra, STROMAL_COL_CANDIDATES)
    ifng_col = find_first_existing_column(extra, IFNG_COL_CANDIDATES)
    tgfb_col = find_first_existing_column(extra, TGFB_COL_CANDIDATES)

    mapped = pd.DataFrame()
    mapped["case_barcode"] = extra[case_col].astype(str).str.strip()

    if proj_col is not None:
        mapped["project_short_name"] = extra[proj_col].apply(standardize_project_name)

    if purity_col is not None:
        mapped["purity"] = pd.to_numeric(extra[purity_col], errors="coerce")
    if lymph_col is not None:
        mapped["immune_lymphoid"] = pd.to_numeric(extra[lymph_col], errors="coerce")
    if myeloid_col is not None:
        mapped["immune_myeloid"] = pd.to_numeric(extra[myeloid_col], errors="coerce")
    if leuk_col is not None:
        mapped["immune_leukocyte"] = pd.to_numeric(extra[leuk_col], errors="coerce")
    if stromal_col is not None:
        mapped["stromal_fraction"] = pd.to_numeric(extra[stromal_col], errors="coerce")
    if ifng_col is not None:
        mapped["ifn_gamma_response"] = pd.to_numeric(extra[ifng_col], errors="coerce")
    if tgfb_col is not None:
        mapped["tgfb_response"] = pd.to_numeric(extra[tgfb_col], errors="coerce")

    if "project_short_name" in mapped.columns:
        group_keys = ["case_barcode", "project_short_name"]
    else:
        group_keys = ["case_barcode"]

    val_cols = [c for c in mapped.columns if c not in group_keys]
    if val_cols:
        mapped = mapped.groupby(group_keys, as_index=False)[val_cols].mean(numeric_only=True)
    else:
        mapped = mapped[group_keys].drop_duplicates().copy()

    return mapped, chosen


def prepare_analysis_df(extra_covariates_path: str | None):
    expr = pd.read_csv(os.path.join(DATA_DIR, "tcga_expanded_tpm.csv"))
    clin = pd.read_csv(os.path.join(DATA_DIR, "tcga_clinical.csv"))

    tumors = expr[
        expr["sample_type_name"].str.contains("Tumor|Metastatic", case=False, na=False)
    ].copy()
    tumors = deduplicate_one_tumor_sample_per_case(tumors)

    merged = tumors.merge(clin, on="case_barcode", how="inner", suffixes=("", "_clin"))
    if "project_short_name_clin" in merged.columns:
        merged.drop(columns=["project_short_name_clin"], inplace=True)

    merged["days_to_death"] = pd.to_numeric(merged["days_to_death"], errors="coerce")
    merged["days_to_last_follow_up"] = pd.to_numeric(
        merged["days_to_last_follow_up"], errors="coerce"
    )

    merged["event"] = (merged["vital_status"] == "Dead").astype(int)
    merged["surv_time"] = np.where(
        merged["event"] == 1,
        merged["days_to_death"],
        merged["days_to_last_follow_up"],
    )
    merged["surv_time"] = pd.to_numeric(merged["surv_time"], errors="coerce")
    merged = merged.dropna(subset=["surv_time"])
    merged = merged[merged["surv_time"] > 0].copy()

    merged["circ_cv"] = compute_circadian_cv(merged)
    merged["circ_cv_ortho"] = compute_circadian_cv_orthogonal(merged)
    merged["bf_class"] = classify_boundary_failure(merged)
    merged["pdl1_log"] = log_transform(merged["CD274"])
    merged["b2m_log"] = log_transform(merged["B2M"])

    merged["age_at_index"] = pd.to_numeric(merged["age_at_index"], errors="coerce")
    merged["sex_male"] = merged["gender"].astype(str).str.lower().eq("male").astype(float)

    merged["stage_norm"] = merged["ajcc_stage"].apply(normalize_stage)
    merged["stage_missing"] = merged["stage_norm"].isna().astype(float)
    merged["stage_II"] = (merged["stage_norm"] == "II").astype(float)
    merged["stage_III"] = (merged["stage_norm"] == "III").astype(float)
    merged["stage_IV"] = (merged["stage_norm"] == "IV").astype(float)

    extra_df, extra_path = load_optional_covariates(extra_covariates_path)
    extra_cols = []
    if extra_df is not None:
        merge_keys = ["case_barcode", "project_short_name"] if "project_short_name" in extra_df.columns else ["case_barcode"]
        merged = merged.merge(extra_df, on=merge_keys, how="left")
        extra_cols = [c for c in OPTIONAL_CONTROL_ORDER if c in merged.columns]

    return merged, extra_cols, extra_path


def _cox_ll_grad_hess(beta, X, time, event):
    """
    Cox partial log-likelihood with Breslow ties.

    Returns: (loglik, grad, hess)
      grad = d(loglik)/d(beta)
      hess = d2(loglik)/d(beta)^2 (negative definite at optimum)
    """
    order = np.argsort(-time, kind="mergesort")
    t = time[order]
    e = event[order]
    Xs = X[order]

    eta = np.clip(Xs @ beta, -50, 50)
    exp_eta = np.exp(eta)

    cum_s0 = np.cumsum(exp_eta)
    cum_s1 = np.cumsum(exp_eta[:, None] * Xs, axis=0)
    outer = np.einsum("ni,nj,n->nij", Xs, Xs, exp_eta)
    cum_s2 = np.cumsum(outer, axis=0)

    loglik = 0.0
    grad = np.zeros(X.shape[1], dtype=float)
    hess = np.zeros((X.shape[1], X.shape[1]), dtype=float)

    event_times = np.unique(t[e == 1])

    for tt in event_times:
        evt_mask = (t == tt) & (e == 1)
        d = int(evt_mask.sum())
        if d == 0:
            continue

        risk_end = np.searchsorted(-t, -tt, side="right") - 1
        if risk_end < 0:
            continue

        s0 = cum_s0[risk_end]
        s1 = cum_s1[risk_end]
        s2 = cum_s2[risk_end]

        x_event_sum = Xs[evt_mask].sum(axis=0)
        mean_x = s1 / s0

        loglik += float(x_event_sum @ beta) - d * np.log(s0)
        grad += x_event_sum - d * mean_x
        hess -= d * (s2 / s0 - np.outer(mean_x, mean_x))

    return loglik, grad, hess


def fit_cox_breslow(
    df: pd.DataFrame,
    duration_col: str,
    event_col: str,
    covariates,
    l2_penalty: float = 1e-3,
):
    """Fit Cox PH model via penalized partial likelihood optimization."""
    cols = [duration_col, event_col] + list(covariates)
    work = df[cols].copy()
    work = work.replace([np.inf, -np.inf], np.nan).dropna()

    if len(work) < MIN_ROWS:
        return None, work, "too_few_rows"

    X_raw = work[covariates].astype(float).values
    time = work[duration_col].astype(float).values
    event = work[event_col].astype(int).values

    if int(event.sum()) < MIN_EVENTS:
        return None, work, "too_few_events"

    # Remove near-constant covariates to avoid singular fits.
    keep_idx = []
    kept_covars = []
    for i, c in enumerate(covariates):
        v = np.nanstd(X_raw[:, i])
        if np.isfinite(v) and v > 1e-10:
            keep_idx.append(i)
            kept_covars.append(c)

    if not keep_idx:
        return None, work, "all_covariates_constant"

    X_raw = X_raw[:, keep_idx]
    p = X_raw.shape[1]
    # Scale non-binary covariates for numerical stability.
    X = X_raw.copy()
    scale = np.ones(p, dtype=float)
    center = np.zeros(p, dtype=float)
    is_binary = np.array(
        [
            set(np.unique(X[:, j][np.isfinite(X[:, j])])) <= {0.0, 1.0}
            for j in range(p)
        ]
    )
    for j in range(p):
        if is_binary[j]:
            continue
        mu = float(np.nanmean(X[:, j]))
        sd = float(np.nanstd(X[:, j]))
        if sd > 1e-12:
            X[:, j] = (X[:, j] - mu) / sd
            center[j] = mu
            scale[j] = sd

    p = X.shape[1]

    def penalized_parts(beta):
        ll, grad, hess = _cox_ll_grad_hess(beta, X, time, event)
        ll_pen = ll - 0.5 * l2_penalty * float(beta @ beta)
        grad_pen = grad - l2_penalty * beta
        hess_pen = hess - l2_penalty * np.eye(p)
        return ll_pen, grad_pen, hess_pen

    # Newton-Raphson with backtracking line search (maximize penalized log-lik).
    beta = np.zeros(p, dtype=float)
    ll_old, grad_old, hess_old = penalized_parts(beta)
    converged = False
    message = "max_iter_reached"

    for _ in range(100):
        info = -hess_old
        try:
            step = np.linalg.solve(info, grad_old)
        except np.linalg.LinAlgError:
            step = np.linalg.pinv(info) @ grad_old

        step_norm = float(np.max(np.abs(step)))
        if step_norm < 1e-6:
            converged = True
            message = "converged_step_tol"
            break

        alpha = 1.0
        improved = False
        while alpha >= 1e-6:
            beta_new = beta + alpha * step
            ll_new, grad_new, hess_new = penalized_parts(beta_new)
            if np.isfinite(ll_new) and ll_new >= ll_old:
                beta = beta_new
                ll_old, grad_old, hess_old = ll_new, grad_new, hess_new
                improved = True
                break
            alpha *= 0.5

        if not improved:
            message = "line_search_failed"
            break

        if float(np.max(np.abs(grad_old))) < 1e-5:
            converged = True
            message = "converged_grad_tol"
            break

    beta_hat = beta
    hess_f = -hess_old
    if not converged:
        # Fallback to BFGS from current estimate.
        opt = optimize.minimize(
            fun=lambda b: -(penalized_parts(b)[0]),
            x0=beta_hat,
            jac=lambda b: -(penalized_parts(b)[1]),
            method="BFGS",
            options={"maxiter": 500, "gtol": 1e-6},
        )
        if opt.success:
            beta_hat = opt.x
            ll_old, grad_old, hess_old = penalized_parts(beta_hat)
            hess_f = -hess_old
            converged = True
            message = f"converged_bfgs:{opt.message}"
        else:
            return None, work, f"optimizer_failed:{message};{opt.message}"

    try:
        cov = np.linalg.pinv(hess_f)
    except np.linalg.LinAlgError:
        return None, work, "hessian_inversion_failed"

    se = np.sqrt(np.clip(np.diag(cov), 1e-15, None))
    # Back-transform coefficients to original covariate scales.
    beta_unscaled = beta_hat / scale
    se_unscaled = se / scale
    z = beta_unscaled / se_unscaled
    pvals = 2 * stats.norm.sf(np.abs(z))

    coef_df = pd.DataFrame(
        {
            "covariate": kept_covars,
            "beta": beta_unscaled,
            "se": se_unscaled,
            "z": z,
            "p_value": pvals,
            "hazard_ratio": np.exp(beta_unscaled),
            "ci95_low": np.exp(beta_unscaled - 1.96 * se_unscaled),
            "ci95_high": np.exp(beta_unscaled + 1.96 * se_unscaled),
            "n_rows": len(work),
            "n_events": int(event.sum()),
            "converged": converged,
        }
    )

    meta = {
        "n_rows": len(work),
        "n_events": int(event.sum()),
        "converged": converged,
        "message": message,
        "covariates_used": kept_covars,
    }
    return {"coef": coef_df, "meta": meta}, work, "ok"


def ph_time_interaction_screen(
    df: pd.DataFrame,
    duration_col: str,
    event_col: str,
    covariates,
    model_label: str,
    project: str,
):
    """
    Approximate PH check: add covariate x log(time) term one-at-a-time.
    Significant interaction suggests non-proportional hazards.
    """
    records = []

    base_cols = [duration_col, event_col] + list(covariates)
    base = df[base_cols].copy().replace([np.inf, -np.inf], np.nan).dropna()
    if len(base) < MIN_ROWS or base[event_col].sum() < MIN_EVENTS:
        return records

    logt = np.log(base[duration_col].clip(lower=1.0))

    for cov in covariates:
        inter_col = f"{cov}_x_logt"
        work = base.copy()
        work[inter_col] = work[cov].astype(float) * logt.values

        fit, _, status = fit_cox_breslow(
            work,
            duration_col=duration_col,
            event_col=event_col,
            covariates=list(covariates) + [inter_col],
        )

        rec = {
            "project": project,
            "model": model_label,
            "covariate": cov,
            "interaction_term": inter_col,
            "status": status,
            "n_rows": len(work),
            "n_events": int(work[event_col].sum()),
            "p_value_interaction": np.nan,
            "ph_violation_p_lt_0_05": False,
        }

        if fit is not None:
            row = fit["coef"].loc[fit["coef"]["covariate"] == inter_col]
            if len(row) == 1:
                p = float(row["p_value"].iloc[0])
                rec["p_value_interaction"] = p
                rec["ph_violation_p_lt_0_05"] = bool(p < 0.05)

        records.append(rec)

    return records


def partial_rank_correlation(df: pd.DataFrame, x_col: str, y_col: str, covars):
    """Partial Spearman-like correlation via rank residualization."""
    cols = [x_col, y_col] + list(covars)
    work = df[cols].copy().replace([np.inf, -np.inf], np.nan).dropna()
    n = len(work)
    m = len(covars)

    if n < max(30, m + 8):
        return None

    xr = work[x_col].rank(method="average").values.astype(float)
    yr = work[y_col].rank(method="average").values.astype(float)

    Z = np.column_stack(
        [
            np.ones(n),
            *[work[c].rank(method="average").values.astype(float) for c in covars],
        ]
    )

    bx, *_ = np.linalg.lstsq(Z, xr, rcond=None)
    by, *_ = np.linalg.lstsq(Z, yr, rcond=None)

    rx = xr - Z @ bx
    ry = yr - Z @ by

    r = np.corrcoef(rx, ry)[0, 1]
    r = float(np.clip(r, -0.999999, 0.999999))

    dof = n - m - 2
    if dof <= 0:
        return None

    t_stat = r * np.sqrt(dof / max(1e-12, 1 - r**2))
    p = 2 * stats.t.sf(abs(t_stat), df=dof)

    return {
        "n_rows": n,
        "r_partial": r,
        "t_stat": t_stat,
        "p_value": p,
    }


def add_primary_test(primary_records, project, test_name, fit, covariate):
    rec = {
        "project": project,
        "test_name": test_name,
        "covariate": covariate,
        "p_value": np.nan,
        "hazard_ratio": np.nan,
        "ci95_low": np.nan,
        "ci95_high": np.nan,
        "n_rows": np.nan,
        "n_events": np.nan,
        "status": "missing",
    }

    if fit is not None:
        coef = fit["coef"]
        row = coef.loc[coef["covariate"] == covariate]
        if len(row) == 1:
            rec.update(
                {
                    "p_value": float(row["p_value"].iloc[0]),
                    "hazard_ratio": float(row["hazard_ratio"].iloc[0]),
                    "ci95_low": float(row["ci95_low"].iloc[0]),
                    "ci95_high": float(row["ci95_high"].iloc[0]),
                    "n_rows": int(row["n_rows"].iloc[0]),
                    "n_events": int(row["n_events"].iloc[0]),
                    "status": "ok",
                }
            )

    primary_records.append(rec)


def run_for_project(df_all: pd.DataFrame, project: str, extra_cols):
    sub = df_all[df_all["project_short_name"] == project].copy()

    coef_records = []
    primary_records = []
    ph_records = []
    partial_records = []

    if len(sub) < MIN_ROWS:
        print(f"[{project}] Skipping: too few rows ({len(sub)})")
        return coef_records, primary_records, ph_records, partial_records

    stage_covars = ["stage_II", "stage_III", "stage_IV", "stage_missing"]
    base_covars = ["age_at_index", "sex_male"] + stage_covars

    # ------------------------------------------------------------------
    # Model 1: AM vs DC adjusted Cox
    # ------------------------------------------------------------------
    amdc = sub[sub["bf_class"].isin(["Active_Masking", "Decoherence"])].copy()
    amdc["boundary_mode"] = (amdc["bf_class"] == "Active_Masking").astype(float)
    extra_covars = [c for c in OPTIONAL_CONTROL_ORDER if c in extra_cols and c in amdc.columns]

    cov_m1 = ["boundary_mode"] + base_covars + extra_covars

    fit_m1, m1_data, status_m1 = fit_cox_breslow(
        amdc,
        duration_col="surv_time",
        event_col="event",
        covariates=cov_m1,
    )

    if fit_m1 is not None:
        tmp = fit_m1["coef"].copy()
        tmp["project"] = project
        tmp["model"] = "cox_am_vs_dc_adjusted"
        coef_records.append(tmp)

        ph_records.extend(
            ph_time_interaction_screen(
                m1_data,
                duration_col="surv_time",
                event_col="event",
                covariates=fit_m1["meta"]["covariates_used"],
                model_label="cox_am_vs_dc_adjusted",
                project=project,
            )
        )

    else:
        print(f"[{project}] cox_am_vs_dc_adjusted failed: {status_m1}")

    add_primary_test(primary_records, project, "am_vs_dc_adjusted", fit_m1, "boundary_mode")

    # ------------------------------------------------------------------
    # Model 1b: AM vs DC adjusted + continuous PD-L1/B2M sensitivity
    # ------------------------------------------------------------------
    cov_m1b = ["boundary_mode", "pdl1_log", "b2m_log"] + base_covars + extra_covars

    fit_m1b, m1b_data, status_m1b = fit_cox_breslow(
        amdc,
        duration_col="surv_time",
        event_col="event",
        covariates=cov_m1b,
    )

    if fit_m1b is not None:
        tmp = fit_m1b["coef"].copy()
        tmp["project"] = project
        tmp["model"] = "cox_am_vs_dc_plus_pdl1_b2m"
        coef_records.append(tmp)

        ph_records.extend(
            ph_time_interaction_screen(
                m1b_data,
                duration_col="surv_time",
                event_col="event",
                covariates=fit_m1b["meta"]["covariates_used"],
                model_label="cox_am_vs_dc_plus_pdl1_b2m",
                project=project,
            )
        )
    else:
        print(f"[{project}] cox_am_vs_dc_plus_pdl1_b2m failed: {status_m1b}")

    add_primary_test(
        primary_records,
        project,
        "am_vs_dc_plus_pdl1_b2m",
        fit_m1b,
        "boundary_mode",
    )

    # ------------------------------------------------------------------
    # Model 2: Continuous PD-L1 x clock interaction
    # ------------------------------------------------------------------
    sub2 = sub.copy()
    sub2["pdl1_x_clock"] = sub2["pdl1_log"] * sub2["circ_cv_ortho"]

    cov_m2 = [
        "pdl1_log",
        "circ_cv_ortho",
        "pdl1_x_clock",
    ] + base_covars + [c for c in extra_covars if c in sub2.columns]

    fit_m2, m2_data, status_m2 = fit_cox_breslow(
        sub2,
        duration_col="surv_time",
        event_col="event",
        covariates=cov_m2,
    )

    if fit_m2 is not None:
        tmp = fit_m2["coef"].copy()
        tmp["project"] = project
        tmp["model"] = "cox_pdl1_clock_interaction"
        coef_records.append(tmp)

        ph_records.extend(
            ph_time_interaction_screen(
                m2_data,
                duration_col="surv_time",
                event_col="event",
                covariates=fit_m2["meta"]["covariates_used"],
                model_label="cox_pdl1_clock_interaction",
                project=project,
            )
        )

    else:
        print(f"[{project}] cox_pdl1_clock_interaction failed: {status_m2}")

    add_primary_test(primary_records, project, "pdl1_x_clock_interaction", fit_m2, "pdl1_x_clock")

    # ------------------------------------------------------------------
    # Model 3: Within high PD-L1 subgroup, does clock still matter?
    # ------------------------------------------------------------------
    pdl1_cut = sub["pdl1_log"].median()
    high = sub[sub["pdl1_log"] >= pdl1_cut].copy()

    cov_m3 = ["circ_cv_ortho"] + base_covars + [c for c in extra_covars if c in high.columns]
    fit_m3, m3_data, status_m3 = fit_cox_breslow(
        high,
        duration_col="surv_time",
        event_col="event",
        covariates=cov_m3,
    )

    if fit_m3 is not None:
        tmp = fit_m3["coef"].copy()
        tmp["project"] = project
        tmp["model"] = "cox_high_pdl1_clock_effect"
        coef_records.append(tmp)

        ph_records.extend(
            ph_time_interaction_screen(
                m3_data,
                duration_col="surv_time",
                event_col="event",
                covariates=fit_m3["meta"]["covariates_used"],
                model_label="cox_high_pdl1_clock_effect",
                project=project,
            )
        )
    else:
        print(f"[{project}] cox_high_pdl1_clock_effect failed: {status_m3}")

    add_primary_test(primary_records, project, "high_pdl1_clock_effect", fit_m3, "circ_cv_ortho")

    # ------------------------------------------------------------------
    # Partial-correlation controls (using whatever controls are available)
    # ------------------------------------------------------------------
    available_controls = [c for c in OPTIONAL_CONTROL_ORDER if c in sub.columns]
    if available_controls:
        out = partial_rank_correlation(sub, "circ_cv", "pdl1_log", available_controls)
        if out is not None:
            partial_records.append(
                {
                    "project": project,
                    "test_name": "partial_corr_circ_cv_vs_pdl1_given_controls",
                    "covariates": ",".join(available_controls),
                    **out,
                }
            )

        out = partial_rank_correlation(sub, "circ_cv_ortho", "pdl1_log", available_controls)
        if out is not None:
            partial_records.append(
                {
                    "project": project,
                    "test_name": "partial_corr_ortho_cv_vs_pdl1_given_controls",
                    "covariates": ",".join(available_controls),
                    **out,
                }
            )

        if "ifn_gamma_response" in available_controls:
            out = partial_rank_correlation(
                sub,
                "circ_cv_ortho",
                "pdl1_log",
                ["ifn_gamma_response"],
            )
            if out is not None:
                partial_records.append(
                    {
                        "project": project,
                        "test_name": "partial_corr_ortho_cv_vs_pdl1_given_ifng",
                        "covariates": "ifn_gamma_response",
                        **out,
                    }
                )

    return coef_records, primary_records, ph_records, partial_records


def apply_bh(df: pd.DataFrame, p_col: str = "p_value") -> pd.DataFrame:
    out = df.copy()
    out["q_value"] = np.nan
    valid = out[p_col].notna().values
    if valid.sum() > 0:
        out.loc[valid, "q_value"] = false_discovery_control(
            out.loc[valid, p_col].astype(float).values,
            method="bh",
        )
    out["significant_fdr05"] = out["q_value"] < 0.05
    return out


def main():
    parser = argparse.ArgumentParser(description="Run TCGA robustness checks for survival claims.")
    parser.add_argument(
        "--projects",
        nargs="+",
        default=["TCGA-SKCM", "TCGA-LUAD"],
        help="Project IDs to analyze (default: TCGA-SKCM TCGA-LUAD)",
    )
    parser.add_argument(
        "--covariates",
        default=None,
        help="Optional path to purity/immune covariates CSV (merged by case_barcode).",
    )
    parser.add_argument(
        "--no-gcs",
        action="store_true",
        help="Skip GCS upload of robustness outputs.",
    )
    args = parser.parse_args()

    projects = [p for p in args.projects if p in ALL_PROJECTS]
    if not projects:
        raise SystemExit("No valid projects selected.")

    print("=" * 74)
    print("  TCGA Robustness Checks: Cox + PH + Microenvironment Controls + Interaction")
    print("=" * 74)

    df_all, extra_cols, extra_path = prepare_analysis_df(args.covariates)
    print(f"Rows with valid survival: {len(df_all)}")
    print(f"Selected projects: {projects}")

    if extra_path:
        print(f"Optional covariates loaded from: {extra_path}")
        if extra_cols:
            print(f"Covariate columns detected: {extra_cols}")
        else:
            print("No numeric control covariates detected in optional file.")
    else:
        print("No optional covariate file found. Running clinical-only controls.")

    coef_frames = []
    primary_rows = []
    ph_rows = []
    partial_rows = []

    for proj in projects:
        n_proj = int((df_all["project_short_name"] == proj).sum())
        print(f"\n[{proj}] rows in analysis frame: {n_proj}")
        c, p, ph, pc = run_for_project(df_all, proj, extra_cols)
        coef_frames.extend(c)
        primary_rows.extend(p)
        ph_rows.extend(ph)
        partial_rows.extend(pc)

    if coef_frames:
        coef_df = pd.concat(coef_frames, ignore_index=True)
    else:
        coef_df = pd.DataFrame(
            columns=[
                "covariate",
                "beta",
                "se",
                "z",
                "p_value",
                "hazard_ratio",
                "ci95_low",
                "ci95_high",
                "n_rows",
                "n_events",
                "converged",
                "project",
                "model",
            ]
        )

    primary_df = pd.DataFrame(primary_rows)

    if ph_rows:
        ph_df = pd.DataFrame(ph_rows)
    else:
        ph_df = pd.DataFrame(
            columns=[
                "project",
                "model",
                "covariate",
                "interaction_term",
                "status",
                "n_rows",
                "n_events",
                "p_value_interaction",
                "ph_violation_p_lt_0_05",
            ]
        )

    if partial_rows:
        partial_df = pd.DataFrame(partial_rows)
    else:
        partial_df = pd.DataFrame(
            columns=[
                "project",
                "test_name",
                "covariates",
                "n_rows",
                "r_partial",
                "t_stat",
                "p_value",
            ]
        )

    if not primary_df.empty:
        primary_df = apply_bh(primary_df, p_col="p_value")

    model_by_test = {
        "am_vs_dc_adjusted": "cox_am_vs_dc_adjusted",
        "am_vs_dc_plus_pdl1_b2m": "cox_am_vs_dc_plus_pdl1_b2m",
        "pdl1_x_clock_interaction": "cox_pdl1_clock_interaction",
        "high_pdl1_clock_effect": "cox_high_pdl1_clock_effect",
    }
    if not primary_df.empty:
        primary_df["model"] = primary_df["test_name"].map(model_by_test)
        if not ph_df.empty:
            ph_summary = (
                ph_df.groupby(["project", "model"], dropna=False)["ph_violation_p_lt_0_05"]
                .agg(["count", "sum"])
                .reset_index()
                .rename(columns={"count": "ph_screen_covariates", "sum": "ph_violation_count"})
            )
            primary_df = primary_df.merge(ph_summary, on=["project", "model"], how="left")
            primary_df["ph_screen_covariates"] = (
                primary_df["ph_screen_covariates"].fillna(0).astype(int)
            )
            primary_df["ph_violation_count"] = (
                primary_df["ph_violation_count"].fillna(0).astype(int)
            )
            with np.errstate(divide="ignore", invalid="ignore"):
                primary_df["ph_violation_fraction"] = np.where(
                    primary_df["ph_screen_covariates"] > 0,
                    primary_df["ph_violation_count"] / primary_df["ph_screen_covariates"],
                    np.nan,
                )
            primary_df["ph_caution"] = primary_df["ph_violation_fraction"] >= 0.25
        else:
            primary_df["ph_screen_covariates"] = 0
            primary_df["ph_violation_count"] = 0
            primary_df["ph_violation_fraction"] = np.nan
            primary_df["ph_caution"] = False

    out_coef = os.path.join(RESULTS_CSV_DIR, "robustness_cox_coefficients.csv")
    out_primary = os.path.join(RESULTS_CSV_DIR, "robustness_primary_tests.csv")
    out_ph = os.path.join(RESULTS_CSV_DIR, "robustness_ph_diagnostics.csv")
    out_partial = os.path.join(RESULTS_CSV_DIR, "robustness_partial_correlations.csv")

    coef_df.to_csv(out_coef, index=False)
    primary_df.to_csv(out_primary, index=False)
    ph_df.to_csv(out_ph, index=False)
    partial_df.to_csv(out_partial, index=False)

    print(f"\nSaved: {out_coef}")
    print(f"Saved: {out_primary}")
    print(f"Saved: {out_ph}")
    print(f"Saved: {out_partial}")

    if not primary_df.empty:
        cols = [
            "project",
            "test_name",
            "model",
            "covariate",
            "p_value",
            "q_value",
            "significant_fdr05",
            "ph_screen_covariates",
            "ph_violation_count",
            "ph_violation_fraction",
            "ph_caution",
            "n_rows",
            "n_events",
            "status",
        ]
        print("\nPrimary robustness tests (BH-adjusted):")
        print(primary_df[cols].to_string(index=False))

    if not args.no_gcs:
        try:
            upload_to_gcs(out_coef, "analysis/robustness_cox_coefficients.csv")
            upload_to_gcs(out_primary, "analysis/robustness_primary_tests.csv")
            upload_to_gcs(out_ph, "analysis/robustness_ph_diagnostics.csv")
            upload_to_gcs(out_partial, "analysis/robustness_partial_correlations.csv")
        except Exception as e:
            print(f"GCS upload skipped: {e}")

    print("\nRobustness analysis complete.")


if __name__ == "__main__":
    main()
