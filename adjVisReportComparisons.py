# -*- coding: utf-8 -*-
"""
Power-law goodness-of-fit (GoF) analyzer for adjacency/correlation matrices.

What it does
------------
1) For each file (normalization) and for every tissue block (A×A), every tissue pair (A×B),
   and the TS / CT aggregates:
   - Clauset-style continuous power-law tail fit (alpha, xmin)
   - KS-D on tail and a bootstrap p-value for power-law plausibility
   - Alternative models on the tail: truncated exponential, truncated lognormal
     -> log-likelihoods, AIC, ΔAIC vs power law, AIC weights
   - "Barabási-style" diagnostics:
       * Hist, CCDF, log–log CCDF
       * Degree vectors k (row+col sums) -> linear binning (12 bins) and log binning
         with linear regressions on:
           - (x, density)       [linear scale]
           - (log10 x, log10 y) [log–log]
           - (x, -log10 y)      [semi-log]
         -> R^2 reported as a quick (but naive) straight-line check

2) Across files:
   - Order preservation of GoF: rank tissue-pair GoF (by KS-D ascending or AIC_weight)
     and compute Spearman ρ between files (+ p-value)
   - Similarity preservation: build pairwise distance matrices between tissue-pair
     distributions within each file (Wasserstein and KS-D), then compare matrices
     across files via a Mantel test (permutation p-values)

3) Exports
   - PDF: summary tables + static figures
   - HTML: interactive plots and all tables

Assumptions
-----------
- CSV schema: first column is row id (e.g., "AC_ENSG..."), other columns are genes with
  prefix before the first "_ENSG" indicating the region/tissue (e.g., "AC_...").
- Values are continuous edge weights/correlations in [0, 1] (NaNs allowed).
- Files correspond to different normalizations.

References (methods inspiration)
--------------------------------
- Clauset, Shalizi, Newman (2009) "Power-law distributions in empirical data"
- Klaus, Yu, Plenz (2011) (notes on pitfalls / discrete cases)
- Barabási, *Network Science* (http://networksciencebook.com/chapter/4) for visual/naive checks

Author: Eden + ChatGPT
"""

from __future__ import annotations

import os
import re
import io
import math
import random
import tempfile
import html
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple, Iterable

import numpy as np
import pandas as pd
import polars as pl
from scipy import stats as sstats
from scipy.spatial.distance import squareform
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.io as pio

# ---- PDF (ReportLab)
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, Image as RLImage, Table, TableStyle, PageBreak
)
from reportlab.lib.utils import ImageReader

# ensure plotly static export works
pio.kaleido.scope.default_format = "png"

# -------------------------
# Small HTML helpers
# -------------------------
def _html_escape(s: str) -> str:
    return html.escape(s, quote=True)

def _table_html(rows: List[List[str]], caption: Optional[str] = None) -> str:
    if not rows:
        return ""
    thead = "".join(f"<th>{_html_escape(c)}</th>" for c in rows[0])
    body = []
    for r in rows[1:]:
        tds = "".join(f"<td>{_html_escape(str(c))}</td>" for c in r)
        body.append(f"<tr>{tds}</tr>")
    cap = f"<caption>{_html_escape(caption)}</caption>" if caption else ""
    return f"""
    <figure class="tbl">
      {cap}
      <table>
        <thead><tr>{thead}</tr></thead>
        <tbody>{''.join(body)}</tbody>
      </table>
    </figure>
    """

def _wrap_fig_html(fig, title: Optional[str] = None) -> str:
    t = f"<h3>{_html_escape(title)}</h3>" if title else ""
    return t + fig.to_html(full_html=False, include_plotlyjs=False, config={"responsive": True})

def _write_html_report(out_html_path: str, html_sections: List[str], title_text: str):
    os.makedirs(os.path.dirname(out_html_path) or ".", exist_ok=True)
    head = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8"/>
<meta name="viewport" content="width=device-width, initial-scale=1"/>
<title>{_html_escape(title_text)}</title>
<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
<style>
  :root {{
    --fg:#111; --muted:#666; --bg:#fff; --acc:#0b7285; --border:#ddd;
  }}
  body {{ font-family: system-ui,-apple-system,Segoe UI,Roboto,sans-serif; margin: 24px; color: var(--fg); background: var(--bg); }}
  h1 {{ margin: 0 0 8px 0; }}
  h2 {{ margin: 28px 0 6px 0; border-bottom: 1px solid var(--border); padding-bottom: 4px; }}
  h3 {{ margin: 18px 0 6px 0; color: var(--acc); }}
  p  {{ color: var(--fg); }}
  .muted {{ color: var(--muted); }}
  .section {{ margin-bottom: 28px; }}
  .tbl table {{ border-collapse: collapse; width: 100%; font-size: 14px; }}
  .tbl th, .tbl td {{ border: 1px solid var(--border); padding: 6px 8px; text-align: right; }}
  .tbl th:first-child, .tbl td:first-child {{ text-align: left; }}
  .tbl caption {{ text-align: left; margin-bottom: 6px; color: var(--muted); }}
  .divider {{ height: 1px; background: var(--border); margin: 24px 0; }}
  .chip {{ display:inline-block; border:1px solid var(--border); padding:2px 6px; border-radius:8px; margin:2px; font-size:12px; color: var(--muted); }}
  .desc {{ font-size: 14px; line-height: 1.4; }}
</style>
</head>
<body>
<header class="section">
  <h1>{_html_escape(title_text)}</h1>
  <p class="muted">Interactive companion to the PDF. Hover, zoom, and toggle legend entries on charts.</p>
</header>
<div class="divider"></div>
"""
    tail = "</body></html>"
    with open(out_html_path, "w", encoding="utf-8") as f:
        f.write(head)
        for sec in html_sections:
            f.write(sec)
        f.write(tail)

# -------------------------
# CSV structure helpers
# -------------------------
def _normalize(s: str) -> str:
    return s.replace("_", "").upper()

def _extract_region_prefix(col: str) -> Optional[str]:
    if "_ENSG" not in col:
        return None
    return col.split("_ENSG", 1)[0]

def _scan_prefix_to_cols(csv_path: str) -> Tuple[str, Dict[str, List[str]]]:
    lf = pl.scan_csv(csv_path)
    cols = lf.columns
    if not cols:
        raise ValueError(f"{csv_path} has no columns.")
    row_id_col = cols[0]
    m: Dict[str, List[str]] = {}
    for c in cols[1:]:
        pre = _extract_region_prefix(c)
        if pre:
            m.setdefault(pre, []).append(c)
    return row_id_col, m

def _scan_rows_by_prefix(csv_path: str, prefixes: List[str], row_id_col: str) -> Dict[str, List[str]]:
    out: Dict[str, List[str]] = {}
    lf = pl.scan_csv(csv_path)
    for pre in prefixes:
        rows = (
            lf.select(row_id_col)
              .filter(pl.col(row_id_col).str.starts_with(f"{pre}_"))
              .collect(streaming=True)[row_id_col]
              .to_list()
        )
        out[pre] = rows
    return out

def _build_canonical_samples(
    csv_paths: Sequence[str],
    requested_pairs: Sequence[Tuple[str, str]],
    max_rows: int,
    max_cols: int,
    seed: int = 42,
) -> Dict[Tuple[str, str], Tuple[List[str], List[str]]]:
    """
    For each (ra,rb), choose SAME row_ids (ra genes) and SAME col_names (rb genes)
    across ALL files (intersection, then sample).
    Keys are normalized tuple (A_norm, B_norm).
    """
    rng = random.Random(seed)
    per_file = []
    for p in csv_paths:
        row_id_col, pref2cols = _scan_prefix_to_cols(p)
        norm2raw = { _normalize(k): k for k in pref2cols.keys() }
        rows_by_pref = _scan_rows_by_prefix(p, list(pref2cols.keys()), row_id_col)
        per_file.append(dict(row_id_col=row_id_col, pref2cols=pref2cols,
                             rows_by_pref=rows_by_pref, norm2raw=norm2raw))
    canon: Dict[Tuple[str, str], Tuple[List[str], List[str]]] = {}
    for (ra, rb) in requested_pairs:
        nra, nrb = _normalize(ra), _normalize(rb)
        inter_rows, inter_cols, ok = None, None, True
        for meta in per_file:
            raw_a = meta["norm2raw"].get(nra)
            raw_b = meta["norm2raw"].get(nrb)
            if raw_a is None or raw_b is None:
                ok = False; break
            rows = set(meta["rows_by_pref"].get(raw_a, []))
            cols = set(meta["pref2cols"].get(raw_b, []))
            inter_rows = rows if inter_rows is None else (inter_rows & rows)
            inter_cols = cols if inter_cols is None else (inter_cols & cols)
        if not ok or not inter_rows or not inter_cols:
            continue
        prng = random.Random((hash((nra, nrb)) ^ seed) & 0x7FFFFFFF)
        rows_list, cols_list = list(inter_rows), list(inter_cols)
        if max_rows and len(rows_list) > max_rows:
            rows_list = prng.sample(rows_list, max_rows)
        if max_cols and len(cols_list) > max_cols:
            cols_list = prng.sample(cols_list, max_cols)
        canon[(nra, nrb)] = (sorted(rows_list), sorted(cols_list))
    return canon

def _fetch_subblock(csv_path: str, row_ids: List[str], col_names: List[str]) -> pd.DataFrame:
    lf = pl.scan_csv(csv_path)
    cols = lf.columns
    if not cols:
        return pd.DataFrame(index=[], columns=[])
    row_id_col = cols[0]
    keep = [c for c in col_names if c in cols]
    if not keep:
        return pd.DataFrame(index=[], columns=[])
    df = (
        lf.select([row_id_col] + keep)
          .filter(pl.col(row_id_col).is_in(row_ids))
          .collect(streaming=True)
          .to_pandas()
          .set_index(row_id_col)
    )
    idx = [r for r in row_ids if r in df.index]
    keep2 = [c for c in keep if c in df.columns]
    return df.loc[idx, keep2]

# -------------------------
# ECDF / CCDF helpers
# -------------------------
def _ecdf(x: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return np.array([]), np.array([])
    xs = np.sort(x)
    F = np.arange(1, xs.size + 1, dtype=float) / xs.size
    return xs, F

def _ccdf(x: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    xs, F = _ecdf(x)
    return xs, 1.0 - F

# -------------------------
# Power-law fits & tests
# -------------------------
def _alpha_mle_continuous(x: np.ndarray, xmin: float) -> float:
    # x > 0 assumed
    n = x.size
    return 1.0 + n / np.sum(np.log(x / xmin))

def _ksD_tail(x: np.ndarray, xmin: float, alpha: float) -> float:
    xs = np.sort(x)
    F_emp = np.arange(1, xs.size + 1, dtype=float) / xs.size
    # model CDF for continuous PL on [xmin, inf)
    F_model = 1.0 - (xs / float(xmin)) ** (1.0 - alpha)
    return float(np.max(np.abs(F_emp - F_model)))

def _fit_powerlaw_tail(x_all: np.ndarray,
                       xmin_grid: Optional[Sequence[float]] = None,
                       min_tail: int = 100) -> Dict[str, float]:
    """
    Clauset-like tail fit for continuous data: choose xmin minimizing KS-D.
    Returns: dict(alpha, xmin, ks_D, n_tail).
    """
    x = np.asarray(x_all, float)
    x = x[np.isfinite(x) & (x > 0)]
    if x.size < min_tail:
        return dict(alpha=np.nan, xmin=np.nan, ks_D=np.nan, n_tail=int(x.size))
    if xmin_grid is None:
        qs = np.linspace(0.70, 0.98, 15)
        xmin_grid = np.unique(np.quantile(x, qs))
    best = None
    for xmin in xmin_grid:
        tail = x[x >= xmin]
        n = tail.size
        if n < min_tail:
            continue
        alpha = _alpha_mle_continuous(tail, xmin)
        D = _ksD_tail(tail, xmin, alpha)
        cand = dict(alpha=float(alpha), xmin=float(xmin), ks_D=float(D), n_tail=int(n))
        if best is None or D < best["ks_D"]:
            best = cand
    return best or dict(alpha=np.nan, xmin=np.nan, ks_D=np.nan, n_tail=0)

def _sample_powerlaw_continuous(n: int, alpha: float, xmin: float, rng: np.random.Generator) -> np.ndarray:
    # inverse-CDF: x = xmin * (1 - U)^(-1/(alpha-1)), U~Uniform(0,1)
    U = rng.uniform(0, 1, size=n)
    return xmin * (1.0 - U) ** (-1.0 / (alpha - 1.0))

def _bootstrap_plausibility_pvalue(x_all: np.ndarray,
                                   B: int,
                                   min_tail: int = 100,
                                   seed: int = 0) -> Tuple[float, Dict[str, float]]:
    """
    Bootstrap p-value for power-law plausibility (Clauset et al.).
    Steps: fit on real data to get (alpha*, xmin*, D_obs).
           For b=1..B: sample tail size n* from PL(alpha*, xmin*), then refit on the sample
           (including searching xmin) -> D_b. p = mean[D_b >= D_obs].
    Returns (p_value, fit_dict).
    """
    fit = _fit_powerlaw_tail(x_all, min_tail=min_tail)
    alpha, xmin, D_obs, n = fit["alpha"], fit["xmin"], fit["ks_D"], fit["n_tail"]
    if not np.isfinite(alpha) or not np.isfinite(xmin) or n < min_tail:
        return np.nan, fit
    rng = np.random.default_rng(seed)
    xs = np.asarray(x_all, float)
    xs = xs[np.isfinite(xs) & (xs >= xmin)]  # keep only tail to set n for bootstrap
    n_tail = xs.size
    if n_tail < min_tail:
        return np.nan, fit
    ge = 0
    for _ in range(B):
        sim = _sample_powerlaw_continuous(n_tail, alpha, xmin, rng)
        # refit with the same procedure
        bfit = _fit_powerlaw_tail(sim, min_tail=min_tail)
        if np.isfinite(bfit["ks_D"]) and bfit["ks_D"] >= D_obs:
            ge += 1
    p = ge / float(B)
    return p, fit

# ---- alternative tail models on [xmin, inf) ----
def _ll_powerlaw_continuous(x: np.ndarray, xmin: float, alpha: float) -> float:
    # pdf: f(x) = (alpha-1)/xmin * (x/xmin)^(-alpha), x>=xmin
    return (x.size * (np.log(alpha - 1.0) - np.log(xmin))
            - alpha * np.sum(np.log(x / xmin)))

def _ll_trunc_exponential(x: np.ndarray, xmin: float) -> Tuple[float, float]:
    # f(x) = lambda * exp(-lambda (x - xmin)), x>=xmin
    lam = 1.0 / np.mean(x - xmin)
    ll = x.size * np.log(lam) - lam * np.sum(x - xmin)
    return ll, lam

def _ll_trunc_lognormal(x: np.ndarray, xmin: float) -> Tuple[float, float, float]:
    # Truncated at xmin: pdf(x) = (1/(x σ √(2π)) exp(-(ln x - μ)^2 / (2σ^2))) / (1 - Φ((ln xmin - μ)/σ))
    y = np.log(x)
    mu = np.mean(y)
    sigma = np.std(y, ddof=1) if y.size > 1 else np.nan
    if not np.isfinite(sigma) or sigma <= 0:
        return -np.inf, mu, sigma
    zmin = (np.log(xmin) - mu) / sigma
    tail_norm = 1.0 - sstats.norm.cdf(zmin)
    if tail_norm <= 0:
        return -np.inf, mu, sigma
    ll = np.sum(-np.log(x) - np.log(sigma) - 0.5*np.log(2*np.pi) - (y - mu)**2/(2*sigma**2) - np.log(tail_norm))
    return ll, mu, sigma

def _aic(ll: float, k: int) -> float:
    return 2*k - 2*ll

# -------------------------
# Binning & naive regressions (Barabási-style diagnostics)
# -------------------------
def _density_from_bins(x: np.ndarray, bins: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    # returns (x_mid, density) normalized to integrate ~1 over range
    hist, edges = np.histogram(x, bins=bins)
    widths = np.diff(edges)
    mids = 0.5*(edges[:-1] + edges[1:])
    dens = hist / (np.sum(hist) * widths + 1e-12)
    return mids, dens

def _linear_regression_R2(x: np.ndarray, y: np.ndarray) -> float:
    # simple OLS y ~ a + b x
    if x.size < 3 or y.size < 3 or np.all(y==y[0]):
        return np.nan
    X = np.vstack([np.ones_like(x), x]).T
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    yhat = X @ beta
    ss_res = np.sum((y - yhat)**2)
    ss_tot = np.sum((y - np.mean(y))**2)
    return 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan

def _binned_R2s_for_k(k: np.ndarray, bins_linear: int = 12, bins_log: int = 12) -> Dict[str, float]:
    k = np.asarray(k, float)
    k = k[np.isfinite(k) & (k > 0)]
    if k.size < max(10, bins_linear):
        return dict(R2_linear=np.nan, R2_loglog=np.nan, R2_semilog=np.nan)
    # linear binning
    lin_bins = np.linspace(np.min(k), np.max(k), bins_linear + 1)
    xm_lin, d_lin = _density_from_bins(k, lin_bins)
    # log binning
    kmin = np.min(k)
    kmax = np.max(k)
    if kmin <= 0:
        kmin = np.min(k[k>0])
    log_bins = np.logspace(np.log10(kmin), np.log10(kmax), bins_log + 1)
    xm_log, d_log = _density_from_bins(k, log_bins)
    # R^2: linear (x, d)
    R2_linear = _linear_regression_R2(xm_lin, d_lin)
    # R^2: log–log (log10 x, log10 d)
    m = (xm_log > 0) & (d_log > 0)
    R2_loglog = _linear_regression_R2(np.log10(xm_log[m]), np.log10(d_log[m]))
    # R^2: semi-log (x, -log10 d)
    m = d_lin > 0
    R2_semilog = _linear_regression_R2(xm_lin[m], -np.log10(d_lin[m]))
    return dict(R2_linear=R2_linear, R2_loglog=R2_loglog, R2_semilog=R2_semilog)

# -------------------------
# Figures
# -------------------------
def _subplots_block(v: np.ndarray, kvec: np.ndarray, title: str) -> go.Figure:
    """
    Make a 2x3 grid:
      (1) Histogram of v
      (2) CCDF(v)
      (3) CCDF(v) log–log
      (4) k linear binning scatter + fit title with R^2
      (5) k log binning scatter (log–log) + R^2
      (6) k linear binning (x, -log10 dens) + R^2
    """
    fig = make_subplots(rows=2, cols=3, subplot_titles=(
        "Histogram(v)", "CCDF(v)", "CCDF(v) log–log",
        "k: linear bin (R² lin)", "k: log bin (R² log–log)", "k: semilog (R²)"
    ))
    # 1: hist
    fig.add_trace(go.Histogram(x=v, nbinsx=60, name="v"), row=1, col=1)
    fig.update_xaxes(title_text="β", row=1, col=1)
    fig.update_yaxes(title_text="Count", row=1, col=1)

    # 2: CCDF
    xs, cc = _ccdf(v)
    fig.add_trace(go.Scatter(x=xs, y=cc, mode="lines", name="CCDF"), row=1, col=2)
    fig.update_xaxes(title_text="β", row=1, col=2)
    fig.update_yaxes(title_text="P(X ≥ x)", row=1, col=2)

    # 3: CCDF log–log
    mpos = xs > 0
    fig.add_trace(go.Scatter(x=xs[mpos], y=cc[mpos], mode="lines", name="CCDF"), row=1, col=3)
    fig.update_xaxes(type="log", title_text="β (log)", row=1, col=3)
    fig.update_yaxes(type="log", title_text="P(X ≥ x) (log)", row=1, col=3)

    # 4–6: k bins
    R2s = _binned_R2s_for_k(kvec)
    # 4: linear bin density
    kpos = kvec[kvec > 0]
    if kpos.size >= 12:
        lin_bins = np.linspace(np.min(kpos), np.max(kpos), 12 + 1)
        xm, d = _density_from_bins(kpos, lin_bins)
        fig.add_trace(go.Scatter(x=xm, y=d, mode="markers", name="dens"), row=2, col=1)
        # naive line fit
        if np.isfinite(R2s["R2_linear"]):
            X = np.vstack([np.ones_like(xm), xm]).T
            beta, *_ = np.linalg.lstsq(X, d, rcond=None); yhat = X @ beta
            fig.add_trace(go.Scatter(x=xm, y=yhat, mode="lines", name="lin fit"), row=2, col=1)
    fig.update_xaxes(title_text="k", row=2, col=1)
    fig.update_yaxes(title_text="density", row=2, col=1)
    fig.layout.annotations[3].text = f'k: linear bin (R²={R2s["R2_linear"]:.3f})'

    # 5: log bin log–log
    if kpos.size >= 12:
        log_bins = np.logspace(np.log10(np.min(kpos)), np.log10(np.max(kpos)), 12 + 1)
        xm2, d2 = _density_from_bins(kpos, log_bins)
        m = (xm2 > 0) & (d2 > 0)
        fig.add_trace(go.Scatter(x=np.log10(xm2[m]), y=np.log10(d2[m]), mode="markers", name="log bins"), row=2, col=2)
        if np.isfinite(R2s["R2_loglog"]):
            X = np.vstack([np.ones(m.sum()), np.log10(xm2[m])]).T
            beta, *_ = np.linalg.lstsq(X, np.log10(d2[m]), rcond=None); yhat = X @ beta
            fig.add_trace(go.Scatter(x=np.log10(xm2[m]), y=yhat, mode="lines", name="log fit"), row=2, col=2)
    fig.update_xaxes(title_text="log10 k", row=2, col=2)
    fig.update_yaxes(title_text="log10 density", row=2, col=2)
    fig.layout.annotations[4].text = f'k: log bin (R²={R2s["R2_loglog"]:.3f})'

    # 6: semilog
    if kpos.size >= 12:
        xm3, d3 = _density_from_bins(kpos, lin_bins if "lin_bins" in locals() else np.linspace(np.min(kpos), np.max(kpos), 12 + 1))
        m3 = d3 > 0
        fig.add_trace(go.Scatter(x=xm3[m3], y=-np.log10(d3[m3]), mode="markers", name="-log10 dens"), row=2, col=3)
        if np.isfinite(R2s["R2_semilog"]):
            X = np.vstack([np.ones(m3.sum()), xm3[m3]]).T
            beta, *_ = np.linalg.lstsq(X, -np.log10(d3[m3]), rcond=None); yhat = X @ beta
            fig.add_trace(go.Scatter(x=xm3[m3], y=yhat, mode="lines", name="semilog fit"), row=2, col=3)
    fig.update_xaxes(title_text="k", row=2, col=3)
    fig.update_yaxes(title_text="-log10 density", row=2, col=3)
    fig.layout.annotations[5].text = f'k: semilog (R²={R2s["R2_semilog"]:.3f})'

    fig.update_layout(title=title, width=1150, height=800,
                      margin=dict(l=60, r=20, t=70, b=60), legend=dict(orientation="h"))
    return fig

# -------------------------
# Metrics & matrices
# -------------------------
def _concat_finite(arrs: Iterable[np.ndarray]) -> np.ndarray:
    xs = []
    for a in arrs:
        a = np.asarray(a, float)
        if a.size:
            xs.append(a[np.isfinite(a)])
    return np.concatenate(xs) if xs else np.array([], dtype=float)

def _degrees_from_block(A: np.ndarray) -> np.ndarray:
    if A.size == 0:
        return np.array([], float)
    Z = A.copy()
    Z[~np.isfinite(Z)] = 0.0
    row_k = np.sum(Z, axis=1)
    col_k = np.sum(Z, axis=0)
    return _concat_finite([row_k, col_k])

def _distance_between_vectors(x: np.ndarray, y: np.ndarray) -> Dict[str, float]:
    m = {}
    if x.size and y.size:
        m["emd"] = float(sstats.wasserstein_distance(x, y))
        ks = sstats.ks_2samp(x, y, alternative="two-sided", method="auto")
        m["ks_D"] = float(ks.statistic)
    else:
        m["emd"] = np.nan; m["ks_D"] = np.nan
    return m

def _mantel_test(D1: np.ndarray, D2: np.ndarray, perms: int = 999, seed: int = 0) -> Tuple[float, float]:
    """
    Mantel test (Spearman) between two distance matrices, upper triangle only.
    Returns (rho, p).
    """
    rng = np.random.default_rng(seed)
    if D1.shape != D2.shape or D1.shape[0] != D1.shape[1]:
        return np.nan, np.nan
    iu = np.triu_indices(D1.shape[0], k=1)
    v1, v2 = D1[iu], D2[iu]
    ok = np.isfinite(v1) & np.isfinite(v2)
    if ok.sum() < 3:
        return np.nan, np.nan
    rho_obs, _ = sstats.spearmanr(v1[ok], v2[ok])
    # permutations of labels
    greater = 0
    for _ in range(perms):
        perm = rng.permutation(D2.shape[0])
        D2p = D2[perm][:, perm]
        v2p = D2p[iu]
        ok2 = np.isfinite(v1) & np.isfinite(v2p)
        if ok2.sum() < 3:
            continue
        rho_p, _ = sstats.spearmanr(v1[ok2], v2p[ok2])
        if rho_p >= rho_obs:
            greater += 1
    p = (greater + 1) / (perms + 1)
    return float(rho_obs), float(p)

# -------------------------
# Main per-block analysis
# -------------------------
@dataclass
class BlockGoF:
    file_label: str
    file_path: str
    region_a: str
    region_b: str
    n_finite: int
    tail_n: int
    alpha: float
    xmin: float
    ks_D_tail: float
    p_plaus: float
    ll_pl: float
    aic_pl: float
    ll_exp: float
    aic_exp: float
    ll_logn: float
    aic_logn: float
    dAIC_exp: float
    dAIC_logn: float
    aicw_pl: float
    aicw_exp: float
    aicw_logn: float
    R2_k_linear: float
    R2_k_loglog: float
    R2_k_semilog: float

def _aic_weights(aics: List[float]) -> List[float]:
    aics = np.asarray(aics, float)
    amin = np.nanmin(aics)
    d = aics - amin
    w = np.exp(-0.5 * d)
    w /= np.nansum(w)
    return list(w)

def _analyze_block(v: np.ndarray,
                   kvec: np.ndarray,
                   B_boot: int = 200,
                   min_tail: int = 100,
                   seed: int = 0) -> Dict[str, float]:
    """
    v: vector of finite β values in the block (A×B)
    kvec: degrees from the same block (row+col sums)
    """
    p_boot, fit = _bootstrap_plausibility_pvalue(v, B=B_boot, min_tail=min_tail, seed=seed)
    alpha, xmin, ksD, n_tail = fit["alpha"], fit["xmin"], fit["ks_D"], fit["n_tail"]
    # compute alternative model fits on the FITTED tail
    vv = np.asarray(v, float)
    tail = vv[np.isfinite(vv) & (vv >= xmin)] if np.isfinite(xmin) else np.array([])
    if tail.size >= min_tail and np.isfinite(alpha) and np.isfinite(xmin):
        ll_pl = _ll_powerlaw_continuous(tail, xmin, alpha)
        ll_exp, lam = _ll_trunc_exponential(tail, xmin)
        ll_logn, mu, sigma = _ll_trunc_lognormal(tail, xmin)
        aic_pl  = _aic(ll_pl,  k=2)
        aic_exp = _aic(ll_exp, k=2)
        aic_lgn = _aic(ll_logn, k=3)
        dAIC_exp  = aic_exp - aic_pl
        dAIC_lgn  = aic_lgn - aic_pl
        w_pl, w_exp, w_lgn = _aic_weights([aic_pl, aic_exp, aic_lgn])
    else:
        ll_pl = ll_exp = ll_logn = -np.inf
        aic_pl = aic_exp = aic_lgn = np.inf
        dAIC_exp = dAIC_lgn = np.nan
        w_pl = w_exp = w_lgn = np.nan

    R2s = _binned_R2s_for_k(kvec)

    return dict(
        n_finite=int(np.isfinite(v).sum()),
        tail_n=int(n_tail),
        alpha=float(alpha),
        xmin=float(xmin),
        ks_D_tail=float(ksD),
        p_plaus=float(p_boot),
        ll_pl=float(ll_pl),
        aic_pl=float(aic_pl),
        ll_exp=float(ll_exp),
        aic_exp=float(aic_exp),
        ll_logn=float(ll_logn),
        aic_logn=float(aic_lgn),
        dAIC_exp=float(dAIC_exp),
        dAIC_logn=float(dAIC_lgn),
        aicw_pl=float(w_pl),
        aicw_exp=float(w_exp),
        aicw_logn=float(w_lgn),
        R2_k_linear=float(R2s["R2_linear"]),
        R2_k_loglog=float(R2s["R2_loglog"]),
        R2_k_semilog=float(R2s["R2_semilog"]),
    )

# -------------------------
# Report builder
# -------------------------
def _image_flowable(png_path: str, max_width: float) -> RLImage:
    ir = ImageReader(png_path)
    iw, ih = ir.getSize()
    scale = min(1.0, max_width / float(iw))
    return RLImage(png_path, width=iw * scale, height=ih * scale)

def _display_name_for_csv(path: str,
                          display_names: Optional[Dict[str, str]] = None,
                          baseline=False) -> str:
    if baseline:
        return "baseline (normalization #1)"
    if not display_names:
        return os.path.basename(path)
    base = os.path.basename(path)
    if path in display_names: return display_names[path]
    if base in display_names: return display_names[base]
    return base

def make_powerlaw_gof_report(
    csv_paths: Sequence[str],
    out_pdf_path: str,
    out_html_path: Optional[str] = None,
    pairs: Optional[Sequence[Tuple[str, str]]] = None,  # None => discover all prefixes; pairs = upper-tri + diagonal
    max_rows: int = 400,
    max_cols: int = 400,
    seed: int = 42,
    B_bootstrap: int = 200,
    min_tail: int = 100,
    file_display_names: Optional[Dict[str, str]] = None,
    overlay_max_pairs_html: int = 24,
) -> List[BlockGoF]:
    """
    Main entry. See module docstring for details.
    """
    os.makedirs(os.path.dirname(out_pdf_path) or ".", exist_ok=True)
    rng = random.Random(seed)

    # discover or validate pairs; build canonical sampling across files
    if pairs is None:
        norm_sets = []
        for p in csv_paths:
            _, pref2cols = _scan_prefix_to_cols(p)
            norm_sets.append(set(_normalize(x) for x in pref2cols.keys()))
        common_norm = set.intersection(*norm_sets) if norm_sets else set()
        prefixes_global = sorted(common_norm)
        requested_pairs_global = [(a, b) for i, a in enumerate(prefixes_global) for b in prefixes_global[i:]]
    else:
        requested_pairs_global = [( _normalize(a), _normalize(b) ) for (a,b) in pairs]

    canonical = _build_canonical_samples(
        csv_paths=csv_paths,
        requested_pairs=requested_pairs_global,
        max_rows=max_rows,
        max_cols=max_cols,
        seed=seed,
    )

    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle(name="Small", fontSize=9, leading=11))
    title_style = styles["Title"]
    h1 = styles["Heading1"]; h2 = styles["Heading2"]
    body = styles["BodyText"]; small = styles["Small"]

    page_w, page_h = A4
    margins = 36
    max_img_width = page_w - 2*margins

    story = []
    html_sections: List[str] = []

    report_title = "Power-law Goodness-of-Fit Report (TS/CT, tissues & pairs)"
    story.append(Paragraph(report_title, title_style))
    story.append(Spacer(1, 0.15*inch))
    story.append(Paragraph(
        "For each normalization file and each tissue/tissue-pair block, we fit a continuous power-law tail "
        "and compute KS-D with bootstrap p-values (Clauset). We compare alternatives (trunc. exponential, "
        "trunc. lognormal) via AIC and give naive Barabási-style R² diagnostics on degree vectors. "
        "Across files, we test rank-order preservation of GoF and preservation of the pairwise similarity "
        "structure via Mantel tests on distance matrices.",
        body
    ))
    story.append(PageBreak())

    if out_html_path:
        html_sections.append(f"""
<section class="section">
  <h2>Overview</h2>
  <p>Files analyzed: <b>{len(csv_paths)}</b>. Canonical sampling made all block comparisons fair.</p>
  <p class="muted">Bootstrap B={B_bootstrap}, min tail size={min_tail}. AIC weights compare power law vs truncated exponential vs truncated lognormal on the fitted tail.</p>
</section>
<div class="divider"></div>
""")

    results: List[BlockGoF] = []
    with tempfile.TemporaryDirectory() as tmpdir:
        # ---------- per-file ----------
        for idx_file, csv_path in enumerate(csv_paths):
            label = _display_name_for_csv(csv_path, file_display_names, baseline=(idx_file==0))
            story.append(Paragraph(label, h1))
            story.append(Spacer(1, 0.06*inch))
            html_sections.append(f'<section class="section"><h2>{_html_escape(label)}</h2>')

            # pairs for this file that exist in canonical
            pairs_this = [(a,b) for (a,b) in requested_pairs_global if (a,b) in canonical]

            # per-file summary rows
            table_rows = [[
                "A","B","n_finite","tail_n","alpha","xmin","KS-D","p(PL)","ΔAIC(exp)","ΔAIC(logn)","AICw(PL)",
                "R²_lin","R²_loglog","R²_semilog"
            ]]

            # analyze each block
            for (a_n, b_n) in pairs_this:
                rows, cols = canonical[(a_n, b_n)]
                df_block = _fetch_subblock(csv_path, rows, cols).to_numpy(float)
                v = df_block[np.isfinite(df_block)]
                kvec = _degrees_from_block(df_block)
                ana = _analyze_block(v, kvec, B_boot=B_bootstrap, min_tail=min_tail, seed=rng.randint(0, 2**31-1))

                results.append(BlockGoF(
                    file_label=label, file_path=csv_path,
                    region_a=a_n, region_b=b_n,
                    n_finite=ana["n_finite"],
                    tail_n=ana["tail_n"],
                    alpha=ana["alpha"], xmin=ana["xmin"],
                    ks_D_tail=ana["ks_D_tail"], p_plaus=ana["p_plaus"],
                    ll_pl=ana["ll_pl"], aic_pl=ana["aic_pl"],
                    ll_exp=ana["ll_exp"], aic_exp=ana["aic_exp"],
                    ll_logn=ana["ll_logn"], aic_logn=ana["aic_logn"],
                    dAIC_exp=ana["dAIC_exp"], dAIC_logn=ana["dAIC_logn"],
                    aicw_pl=ana["aicw_pl"], aicw_exp=ana["aicw_exp"], aicw_logn=ana["aicw_logn"],
                    R2_k_linear=ana["R2_k_linear"], R2_k_loglog=ana["R2_k_loglog"], R2_k_semilog=ana["R2_k_semilog"],
                ))

                table_rows.append([
                    a_n, b_n,
                    f"{ana['n_finite']:,}", f"{ana['tail_n']:,}",
                    f"{ana['alpha']:.3f}", f"{ana['xmin']:.4f}",
                    f"{ana['ks_D_tail']:.3f}", f"{ana['p_plaus']:.3f}",
                    f"{ana['dAIC_exp']:+.2f}", f"{ana['dAIC_logn']:+.2f}",
                    f"{ana['aicw_pl']:.2f}",
                    f"{ana['R2_k_linear']:.3f}", f"{ana['R2_k_loglog']:.3f}", f"{ana['R2_k_semilog']:.3f}",
                ])

                # Figures: subplots for a handful (to keep PDF short); HTML gets all
                fig = _subplots_block(v, kvec, f"{a_n}×{b_n} — GoF diagnostics")
                if out_html_path:
                    html_sections.append(_wrap_fig_html(fig, title=f"{a_n}×{b_n} — diagnostics"))

            # Summary table (this file)
            tbl = Table(table_rows, hAlign="LEFT")
            tbl.setStyle(TableStyle([
                ("FONTNAME", (0,0), (-1,0), "Helvetica-Bold"),
                ("BACKGROUND", (0,0), (-1,0), colors.lightgrey),
                ("GRID", (0,0), (-1,-1), 0.3, colors.grey),
                ("ALIGN", (2,1), (-1,-1), "RIGHT"),
            ]))
            story.append(tbl)
            story.append(PageBreak())
            html_sections.append(_table_html(table_rows, caption="Per-block GoF (this file)"))
            html_sections.append("</section><div class='divider'></div>")

        # ---------- Across files: rank order preservation ----------
        if len(csv_paths) >= 2:
            story.append(Paragraph("Across-files: Order preservation of GoF", h1))
            story.append(Spacer(1, 0.06*inch))

            # build per-file ranked lists by KS-D (ascending) and by AIC weight (descending)
            df_res = pd.DataFrame([r.__dict__ for r in results])
            df_res["pair"] = df_res["region_a"] + "×" + df_res["region_b"]

            # Only keep pairs present in all files
            counts = df_res.groupby("pair")["file_label"].nunique()
            keep_pairs = set(counts[counts == len(csv_paths)].index)
            df_shared = df_res[df_res["pair"].isin(keep_pairs)].copy()

            # ranks per file
            rank_tables = []
            for metric, asc in [("ks_D_tail", True), ("aicw_pl", False)]:
                # Spearman between all file pairs
                labels = [ _display_name_for_csv(p, file_display_names, baseline=(i==0)) for i,p in enumerate(csv_paths) ]
                R = []
                for i in range(len(labels)):
                    for j in range(i+1, len(labels)):
                        a = labels[i]; b = labels[j]
                        A = df_shared[df_shared.file_label==a].set_index("pair")[metric].rank(ascending=asc, method="average")
                        B = df_shared[df_shared.file_label==b].set_index("pair")[metric].rank(ascending=asc, method="average")
                        idx = A.index.intersection(B.index)
                        rho, pval = sstats.spearmanr(A.loc[idx], B.loc[idx])
                        R.append([metric, a, b, f"{rho:.3f}", f"{pval:.3g}", len(idx)])
                rank_tables.extend(R)

            rows = [["Metric","File A","File B","Spearman ρ","p","n pairs"]] + rank_tables
            tbl = Table(rows, hAlign="LEFT")
            tbl.setStyle(TableStyle([
                ("FONTNAME", (0,0), (-1,0), "Helvetica-Bold"),
                ("BACKGROUND", (0,0), (-1,0), colors.HexColor("#EFEFEF")),
                ("GRID", (0,0), (-1,-1), 0.3, colors.grey),
                ("ALIGN", (3,1), (-1,-1), "RIGHT"),
            ]))
            story.append(tbl); story.append(PageBreak())
            if out_html_path:
                html_sections.append("<section class='section'><h2>Across-files: Order preservation</h2>")
                html_sections.append(_table_html(rows, caption="Rank-order preservation (KS-D ascending; AICw power-law descending)"))
                html_sections.append("</section><div class='divider'></div>")

            # ---------- Across files: similarity preservation via Mantel ----------
            story.append(Paragraph("Across-files: Similarity structure preservation (Mantel)", h1))
            story.append(Spacer(1, 0.06*inch))

            # Build per-file distance matrices between pairs (using the distributions v of blocks)
            # We'll use Wasserstein and KS-D; distance between PAIRs computed on their value vectors.
            pair_list = sorted(keep_pairs)
            label_to_index = {lab:i for i,lab in enumerate([ _display_name_for_csv(p, file_display_names, baseline=(i==0)) for i,p in enumerate(csv_paths) ])}

            # gather block value vectors per (file, pair)
            vals_map: Dict[Tuple[str,str], np.ndarray] = {}
            for p in csv_paths:
                lab = _display_name_for_csv(p, file_display_names, baseline=(p==csv_paths[0]))
                for pr in pair_list:
                    a_n, b_n = pr.split("×")
                    rows, cols = canonical[(a_n, b_n)]
                    A = _fetch_subblock(p, rows, cols).to_numpy(float)
                    v = A[np.isfinite(A)]
                    vals_map[(lab, pr)] = v

            def _distance_matrix(labels_pairs: List[str], file_label: str, which: str) -> np.ndarray:
                n = len(labels_pairs)
                D = np.zeros((n, n), float)
                for i in range(n):
                    for j in range(i+1, n):
                        x = vals_map[(file_label, labels_pairs[i])]
                        y = vals_map[(file_label, labels_pairs[j])]
                        m = _distance_between_vectors(x, y)
                        D[i,j] = D[j,i] = m["emd"] if which=="emd" else m["ks_D"]
                return D

            labels_files = [ _display_name_for_csv(p, file_display_names, baseline=(i==0)) for i,p in enumerate(csv_paths) ]
            mantel_rows = [["Distance","File A","File B","Mantel ρ (Spearman)","p","n_pairs"]]
            if len(pair_list) >= 4:
                for dist_name in ["emd","ks_D"]:
                    for i in range(len(labels_files)):
                        for j in range(i+1, len(labels_files)):
                            A_lab, B_lab = labels_files[i], labels_files[j]
                            DA = _distance_matrix(pair_list, A_lab, dist_name)
                            DB = _distance_matrix(pair_list, B_lab, dist_name)
                            rho, p = _mantel_test(DA, DB, perms=999, seed=seed+123)
                            mantel_rows.append([dist_name, A_lab, B_lab, f"{rho:.3f}", f"{p:.3g}", len(pair_list)])
            tbl = Table(mantel_rows, hAlign="LEFT")
            tbl.setStyle(TableStyle([
                ("FONTNAME", (0,0), (-1,0), "Helvetica-Bold"),
                ("BACKGROUND", (0,0), (-1,0), colors.HexColor("#F3F3F3")),
                ("GRID", (0,0), (-1,-1), 0.3, colors.grey),
                ("ALIGN", (3,1), (-1,-1), "RIGHT"),
            ]))
            story.append(tbl); story.append(PageBreak())
            if out_html_path:
                html_sections.append("<section class='section'><h2>Across-files: Similarity preservation (Mantel)</h2>")
                html_sections.append(_table_html(mantel_rows, caption="Mantel tests on pairwise distance matrices of pairs"))
                html_sections.append("</section><div class='divider'></div>")

        # write HTML if requested
        if out_html_path:
            _write_html_report(out_html_path, html_sections, report_title)

        # build PDF
        doc = SimpleDocTemplate(
            out_pdf_path, pagesize=A4,
            leftMargin=margins, rightMargin=margins,
            topMargin=margins, bottomMargin=margins
        )
        doc.build(story)

    return results

# -------------------------
# Example usage
# -------------------------
if __name__ == "__main__":
    """
    Example:
    python powerlaw_gof_report.py

    Edit the CSV paths and output locations below.
    """
    from datetime import datetime

    # Example CSVs (replace with your paths)
    # Each CSV is a different normalization of the same adjacency/correlation matrix schema.
    CSV_CT2_TS3 = "/media/psylab-6028/DATA/Eden/CoExpression_ReProduction/nbs/xwgcna_rosmap_constBeta_CT2_TS3_adjacency.csv"
    CSV_CT1_TS1 = "/media/psylab-6028/DATA/Eden/CoExpression_ReProduction/nbs/xwgcna_rosmap_constBeta_CT1_TS1_adjacency.csv"
    CSV_CT3_TS2 = "/media/psylab-6028/DATA/Eden/CoExpression_ReProduction/notebooks/WGCNA_new_ROSMAP_AdjFunc_adjacency_no_NEWFUNC_ISONEW_zerDiag2025-09-08.csv"


    # Optional friendly names
    file_names = {
        os.path.basename(CSV_CT1_TS1): "Norm A (CT=1,TS=1)",
        os.path.basename(CSV_CT2_TS3): "Norm B (CT=2,TS=3)"
    }

    out_dir = "./outputs_powerlaw_gof"
    os.makedirs(out_dir, exist_ok=True)

    _ = make_powerlaw_gof_report(
        csv_paths=[CSV_CT1_TS1, CSV_CT2_TS3],
        out_pdf_path=os.path.join(out_dir, f"powerlaw_gof_{datetime.now().strftime('%Y%m%d_%H%M%S')}.pdf"),
        out_html_path=os.path.join(out_dir, "powerlaw_gof.html"),
        pairs=None,                # None => auto-discover all prefixes, use all diagonal+upper-triangle pairs shared by all files
        max_rows=1000,             # canonical cap per block
        max_cols=1000,
        seed=123,
        B_bootstrap=200,           # increase for more stable p-values
        min_tail=100,              # min tail points for a valid fit & bootstrap
        file_display_names=file_names,
    )
