"""
Statistical helpers shared across roughness-scale notebooks.

Usage (in a notebook)
---------------------
    # src/ is already on sys.path after the labels import cell
    from stats import rmse, r2, mae, fit_ols, eval_yhat, best_residual_correlate
"""

import numpy as np
import pandas as pd


# ── basic scalar metrics ──────────────────────────────────────────────────────

def rmse(pred, obs):
    """Root-mean-square error, ignoring non-finite pairs."""
    pred, obs = np.asarray(pred, float), np.asarray(obs, float)
    m = np.isfinite(pred) & np.isfinite(obs)
    return float(np.sqrt(np.mean((pred[m] - obs[m]) ** 2))) if m.sum() > 0 else np.nan


def r2(pred, obs):
    """Pearson R² between pred and obs, ignoring non-finite pairs.
    Returns NaN if fewer than 3 valid pairs."""
    pred, obs = np.asarray(pred, float), np.asarray(obs, float)
    m = np.isfinite(pred) & np.isfinite(obs)
    if m.sum() < 3:
        return np.nan
    return float(np.corrcoef(pred[m], obs[m])[0, 1] ** 2)


def mae(pred, obs):
    """Mean absolute error, ignoring NaNs."""
    pred, obs = np.asarray(pred, float), np.asarray(obs, float)
    return float(np.nanmean(np.abs(pred - obs)))


# ── OLS helpers ───────────────────────────────────────────────────────────────

def fit_ols(X, y):
    """Fit OLS via least squares.

    Parameters
    ----------
    X : 2-D array  (must include intercept column)
    y : 1-D array

    Returns
    -------
    beta : coefficient vector
    r2   : R²
    adjr2 : adjusted R²
    yhat : fitted values
    """
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    yhat  = X @ beta
    resid = y - yhat
    n, p  = X.shape
    sse   = float(np.sum(resid ** 2))
    sst   = float(np.sum((y - y.mean()) ** 2))
    r2v   = 1.0 - (sse / sst) if sst > 0 else np.nan
    adjr2 = 1.0 - (1.0 - r2v) * (n - 1) / max(n - p, 1)
    return beta, r2v, adjr2, yhat


def eval_yhat(yhat, y):
    """Compute a summary dict of fit quality metrics.

    Returns
    -------
    dict with keys: bias, mae, rmse, r2
    """
    resid = np.asarray(yhat, float) - np.asarray(y, float)
    sse   = float(np.sum(resid ** 2))
    sst   = float(np.sum((y - np.mean(y)) ** 2))
    r2v   = 1.0 - (sse / sst) if sst > 0 else np.nan
    return {
        "bias": float(np.mean(resid)),
        "mae":  float(np.mean(np.abs(resid))),
        "rmse": float(np.sqrt(np.mean(resid ** 2))),
        "r2":   r2v,
    }


# ── residual–correlate scanners ───────────────────────────────────────────────

def best_residual_correlate(summary, resid, mask=None, sim_cols=None, min_n=20):
    """Scan simulation columns and return those most correlated with *resid*.

    Parameters
    ----------
    summary : DataFrame
    resid : 1-D array
        Pre-filtered residual array (length == mask.sum() if mask is given,
        else length == len(summary)).
    mask : bool array of len(summary), optional
        Applied to each sim column before correlating with *resid*.
    sim_cols : list of str, optional
        Columns to scan.  Default: standard set of design parameters.
    min_n : int
        Minimum valid pairs required to include a column.

    Returns
    -------
    list of (col, r) sorted by |r| descending.
    """
    if sim_cols is None:
        sim_cols = ['fV', 'sigma', 'p', 'tr', 'l', 'aniso',
                    'So', 'Ks_v', 'alpha_v', 'alpha_b']
    resid = np.asarray(resid, float)
    candidates = []
    for c in sim_cols:
        if c not in summary.columns:
            continue
        try:
            v = pd.to_numeric(summary[c], errors='coerce').to_numpy(float)
            if mask is not None:
                v = v[mask]
            m2 = np.isfinite(v) & np.isfinite(resid)
            if m2.sum() < min_n or v[m2].std() < 1e-10:
                continue
            r = float(np.corrcoef(v[m2], resid[m2])[0, 1])
            if np.isfinite(r):
                candidates.append((c, r))
        except Exception:
            pass
    candidates.sort(key=lambda x: abs(x[1]), reverse=True)
    return candidates


def best_mean_correlate(summary, resid_df, sim_cols=None, min_n=20):
    """Scan sim columns; return mean |r| across multiple residual columns.

    Parameters
    ----------
    summary : DataFrame
    resid_df : DataFrame
        Each column is a residual series (same index as summary).
    sim_cols : list of str, optional
    min_n : int

    Returns
    -------
    dict of {col: mean_abs_r} sorted by value descending.
    """
    if sim_cols is None:
        sim_cols = ['fV', 'sigma', 'p', 'tr', 'l', 'aniso',
                    'So', 'Ks_v', 'alpha_v', 'alpha_b']
    scores = {}
    for c in sim_cols:
        if c not in summary.columns:
            continue
        try:
            v = pd.to_numeric(summary[c], errors='coerce')
            rs = []
            for ec in resid_df.columns:
                resid = resid_df[ec]
                m = v.notna() & resid.notna() & np.isfinite(resid)
                if m.sum() < min_n or v[m].std() < 1e-10:
                    continue
                r = float(np.corrcoef(v[m], resid[m])[0, 1])
                if np.isfinite(r):
                    rs.append(abs(r))
            if rs:
                scores[c] = float(np.mean(rs))
        except Exception:
            pass
    return dict(sorted(scores.items(), key=lambda x: -x[1]))
