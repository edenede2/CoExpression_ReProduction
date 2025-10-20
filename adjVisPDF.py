# -*- coding: utf-8 -*-
"""
Region-pair submatrix visualizer and PDF report builder.

Author: Eden + ChatGPT
"""

from __future__ import annotations
import os
import re
import io
import math
import random
import tempfile
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple, Iterable

import html

import numpy as np
import pandas as pd
import polars as pl
from scipy import stats as sstats
import plotly.express as px
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
# Helpers & data structures
# -------------------------
# ==== NEW: naming helper ====
def _display_name_for_csv(path: str,
                          display_names: Optional[Dict[str, str]] = None,
                          baseline=False) -> str:
    """
    Prefer user-provided display names by matching on full path or basename.
    Baseline is always shown as 'baseline TS=1, CT=1' for clarity.
    """
    if baseline:
        return "baseline TS=1, CT=1"
    if not display_names:
        return _label_for_csv(path)
    base = os.path.basename(path)
    if path in display_names:
        return display_names[path]
    if base in display_names:
        return display_names[base]
    return _label_for_csv(path)

# ==== NEW: ECDF/CCDF utilities ====
def _ecdf(x: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    x = np.asarray(x)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return np.array([]), np.array([])
    xs = np.sort(x)
    F = np.arange(1, xs.size + 1, dtype=float) / xs.size
    return xs, F

def _ccdf(x: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    xs, F = _ecdf(x)
    return xs, 1.0 - F

# ==== NEW: CCDF figures (linear & log-log) ====
def _ccdf_figs(vals: np.ndarray, title_prefix: str):
    xs, cc = _ccdf(vals)
    # linear CCDF
    import plotly.graph_objects as go
    fig1 = go.Figure()
    fig1.add_trace(go.Scatter(x=xs, y=cc, mode="lines", name="Empirical CCDF"))
    fig1.update_layout(
        title=f"{title_prefix} — CCDF",
        xaxis_title="K (value)",
        yaxis_title="P(X ≥ x)",
        width=900, height=420, margin=dict(l=60, r=20, t=60, b=50),
    )
    # log-log CCDF (positive x only)
    x_pos = xs[xs > 0]
    cc_pos = cc[xs > 0]
    fig2 = go.Figure()
    if x_pos.size:
        fig2.add_trace(go.Scatter(x=x_pos, y=cc_pos, mode="lines", name="Empirical CCDF"))
        fig2.update_xaxes(type="log")
        fig2.update_yaxes(type="log")
    fig2.update_layout(
        title=f"{title_prefix} — CCDF (log–log)",
        xaxis_title="K (value, log)",
        yaxis_title="P(X ≥ x) (log)",
        width=900, height=420, margin=dict(l=60, r=20, t=60, b=50),
    )
    return fig1, fig2

# ==== NEW: simple continuous power-law tail fit (Clauset-style) ====
def _fit_powerlaw(vals: np.ndarray,
                  xmin_grid: Optional[Sequence[float]] = None,
                  min_tail: int = 100) -> Dict[str, float]:
    """
    Fit a continuous power-law tail: x ~ x^{-alpha} for x >= xmin.
    Returns dict(alpha, xmin, ks_D, n_tail). No p-values are computed.
    """
    x = np.asarray(vals, float)
    x = x[np.isfinite(x) & (x > 0)]
    if x.size < min_tail:
        return dict(alpha=np.nan, xmin=np.nan, ks_D=np.nan, n_tail=int(x.size))

    if xmin_grid is None:
        # candidate xmin over upper quantiles (focus on tail)
        qs = np.linspace(0.70, 0.98, 15)
        xmin_grid = np.unique(np.quantile(x, qs))

    best = None
    for xmin in xmin_grid:
        tail = x[x >= xmin]
        n = tail.size
        if n < min_tail:
            continue
        # MLE for continuous power-law exponent
        alpha = 1.0 + n / np.sum(np.log(tail / float(xmin)))
        # KS distance between empirical CDF of tail and model CDF
        xs, F_emp = _ecdf(tail)
        # model CDF for x>=xmin: F(x) = 1 - (x/xmin)^{1-alpha}
        F_model = 1.0 - (xs / float(xmin)) ** (1.0 - alpha)
        D = float(np.max(np.abs(F_emp - F_model)))
        cand = dict(alpha=float(alpha), xmin=float(xmin), ks_D=D, n_tail=int(n))
        if best is None or D < best["ks_D"]:
            best = cand

    return best or dict(alpha=np.nan, xmin=np.nan, ks_D=np.nan, n_tail=0)

# ==== NEW: overlay power-law line on CCDF ====
def _add_powerlaw_to_ccdf(fig, alpha: float, xmin: float, label: str = "Power-law fit"):
    if not np.isfinite(alpha) or not np.isfinite(xmin):
        return fig
    import plotly.graph_objects as go
    # model CCDF: S(x) = (x/xmin)^{1-alpha} for x>=xmin
    xmax = float(np.nanmax(fig.data[0].x)) if fig.data and len(fig.data[0].x) else 1.0
    x_line = np.linspace(xmin, xmax, 200)
    S = (x_line / xmin) ** (1.0 - alpha)
    fig.add_trace(go.Scatter(x=x_line, y=S, mode="lines", name=label, line=dict(dash="dash")))
    return fig

# ==== NEW: compact HTML helpers ====
def _chip_list(values: Sequence[str]) -> str:
    return " ".join(f'<span class="chip">{_html_escape(v)}</span>' for v in values)

def _h3(s: str) -> str:
    return f"<h3>{_html_escape(s)}</h3>"

def _scan_prefix_to_cols(csv_path: str) -> Tuple[str, Dict[str, List[str]]]:
    lf = pl.scan_csv(csv_path)
    cols = lf.columns
    if not cols:
        raise ValueError(f"{csv_path} has no columns.")
    row_id_col = cols[0]
    data_cols = [c for c in cols if c != row_id_col]
    m: Dict[str, List[str]] = {}
    for c in data_cols:
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
    For each (ra,rb), choose the SAME row_ids (ra-genes) and SAME col_names (rb-genes)
    across ALL files by sampling from the intersection across files.
    Returns: {(NORMA(ra),NORM(rb)) -> (row_ids, col_names)}
    """
    rng = random.Random(seed)
    # Pre-scan per file
    per_file = []
    for p in csv_paths:
        row_id_col, pref2cols = _scan_prefix_to_cols(p)
        norm2raw = { _normalize(k): k for k in pref2cols.keys() }
        rows_by_pref = _scan_rows_by_prefix(p, list(pref2cols.keys()), row_id_col)
        per_file.append(dict(path=p, row_id_col=row_id_col, pref2cols=pref2cols,
                             rows_by_pref=rows_by_pref, norm2raw=norm2raw))

    canon: Dict[Tuple[str, str], Tuple[List[str], List[str]]] = {}
    for (ra, rb) in requested_pairs:
        nra, nrb = _normalize(ra), _normalize(rb)
        # gather intersections
        inter_rows = None
        inter_cols = None
        ok = True
        for meta in per_file:
            raw_a = meta["norm2raw"].get(nra)
            raw_b = meta["norm2raw"].get(nrb)
            if raw_a is None or raw_b is None:
                ok = False
                break
            rows = set(meta["rows_by_pref"].get(raw_a, []))
            cols = set(meta["pref2cols"].get(raw_b, []))
            inter_rows = rows if inter_rows is None else (inter_rows & rows)
            inter_cols = cols if inter_cols is None else (inter_cols & cols)
        if not ok or not inter_rows or not inter_cols:
            continue  # skip pairs that aren't shared across all files

        # deterministic per-pair RNG
        prng = random.Random((hash((nra, nrb)) ^ seed) & 0x7FFFFFFF)
        rows_list = list(inter_rows); cols_list = list(inter_cols)
        if max_rows and len(rows_list) > max_rows:
            rows_list = prng.sample(rows_list, max_rows)
        if max_cols and len(cols_list) > max_cols:
            cols_list = prng.sample(cols_list, max_cols)
        canon[(nra, nrb)] = (sorted(rows_list), sorted(cols_list))
    return canon

def _html_escape(s: str) -> str:
    return html.escape(s, quote=True)

def _table_html(rows: List[List[str]], caption: Optional[str] = None) -> str:
    """Render a simple HTML table from a list of rows."""
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
        <tbody>
          {''.join(body)}
        </tbody>
      </table>
    </figure>
    """

def _wrap_fig_html(fig, title: Optional[str] = None) -> str:
    """Return Plotly figure as an embeddable HTML snippet (no duplicate JS)."""
    t = f"<h3>{_html_escape(title)}</h3>" if title else ""
    # include_plotlyjs=False because we'll load it once in the page head
    return t + fig.to_html(full_html=False, include_plotlyjs=False, config={"responsive": True})

def _pair_key(ra: str, rb: str) -> Tuple[str, str]:
    # Keep orientation (A,B) as-is; key is normalized region names
    return (_normalize(ra), _normalize(rb))

def _fetch_subblock(csv_path: str, row_ids: List[str], col_names: List[str]) -> pd.DataFrame:
    """Fetch a tiny aligned submatrix for given row_ids × col_names."""
    lf = pl.scan_csv(csv_path)
    cols = lf.columns
    if not cols:
        raise ValueError(f"{csv_path} has no columns.")
    row_id_col = cols[0]
    avail_cols = [c for c in col_names if c in cols]
    if not avail_cols:
        return pd.DataFrame(index=[], columns=[])
    df = (
        lf.select([row_id_col] + avail_cols)
          .filter(pl.col(row_id_col).is_in(row_ids))
          .collect(streaming=True)
          .to_pandas()
    )
    df = df.set_index(row_id_col)
    # intersect then reindex to desired order
    idx = [r for r in row_ids if r in df.index]
    cols_keep = [c for c in col_names if c in df.columns]
    return df.loc[idx, cols_keep]


def _get_ct_ts_for_pair(
    csv_path: str,
    ra: str,
    rb: str,
    overrides: Optional[dict],
    default_ct: Optional[int],
    default_ts: Optional[int],
) -> Tuple[Optional[int], Optional[int]]:
    """
    Resolve (CT, TS) for (ra, rb) given:
      - overrides[csv_key][(ra, rb)] = (ct, ts)
      - csv_key can be full path or basename
      - (ra, rb) match ignores underscores & case
      - supports a file-level default via keys '*' or '__default__'
      - falls back to (default_ct, default_ts) if nothing matches
    """
    def _norm(s: str) -> str:
        return s.replace("_", "").upper()

    csv_keys = (csv_path, os.path.basename(csv_path))
    ct, ts = default_ct, default_ts

    if not overrides:
        return ct, ts

    # look up by either full path or basename
    file_map = None
    for k in csv_keys:
        if k in overrides:
            file_map = overrides[k]
            break
    if not file_map:
        return ct, ts

    # exact pair (directional) first, then reversed, then default
    ra_n, rb_n = _norm(ra), _norm(rb)

    # try tuple keys
    for key in file_map.keys():
        if isinstance(key, tuple) and len(key) == 2:
            a, b = key
            if _norm(a) == ra_n and _norm(b) == rb_n:
                cto, tso = file_map[key]
                return (cto if cto is not None else ct, tso if tso is not None else ts)

    # try reversed direction
    for key in file_map.keys():
        if isinstance(key, tuple) and len(key) == 2:
            a, b = key
            if _norm(a) == rb_n and _norm(b) == ra_n:
                cto, tso = file_map[key]
                return (cto if cto is not None else ct, tso if tso is not None else ts)

    # file-level default
    for default_key in ("*", "__default__"):
        if default_key in file_map:
            cto, tso = file_map[default_key]
            return (cto if cto is not None else ct, tso if tso is not None else ts)

    return ct, ts

def _normalize(s: str) -> str:
    return s.replace("_", "").upper()

def _extract_region_prefix(col: str) -> Optional[str]:
    # Region prefix is the substring before the first occurrence of "_ENSG"
    if "_ENSG" not in col:
        return None
    return col.split("_ENSG", 1)[0]

def _parse_ct_ts_from_path(p: str) -> Tuple[Optional[int], Optional[int]]:
    """
    Try to parse CT/TS integers from filename like ...CT2_TS3...
    """
    m = re.search(r"CT(\d+).*?TS(\d+)", os.path.basename(p), flags=re.IGNORECASE)
    if not m:
        return None, None
    return int(m.group(1)), int(m.group(2))

def _label_for_csv(path: str) -> str:
    base = os.path.basename(path)
    ct, ts = _parse_ct_ts_from_path(base)
    if ct is not None and ts is not None:
        return f"{base} (CT={ct}, TS={ts})"
    return base

def _image_flowable(png_path: str, max_width: float) -> RLImage:
    ir = ImageReader(png_path)
    iw, ih = ir.getSize()
    scale = min(1.0, max_width / float(iw))
    return RLImage(png_path, width=iw * scale, height=ih * scale)

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
    --fg: #111; --muted: #666; --bg: #fff; --acc: #0b7285;
    --border: #ddd;
  }}
  body {{ font-family: system-ui, -apple-system, Segoe UI, Roboto, sans-serif; margin: 24px; color: var(--fg); background: var(--bg); }}
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
  .toc a {{ text-decoration: none; color: var(--acc); }}
  .toc li {{ margin: 2px 0; }}
  .kbd {{ font-family: ui-monospace, SFMono-Regular, Menlo, monospace; background:#f6f6f6; border:1px solid #eee; padding:0 4px; border-radius:4px; }}
</style>
</head>
<body>
<header class="section">
  <h1>{_html_escape(title_text)}</h1>
  <p class="muted">Interactive version of the PDF report. Hover, zoom, and toggle legend entries on charts.</p>
</header>
<div class="divider"></div>
"""
    tail = "</body></html>"
    with open(out_html_path, "w", encoding="utf-8") as f:
        f.write(head)
        for sec in html_sections:
            f.write(sec)
        f.write(tail)
def _metrics_pair(x: np.ndarray, y: np.ndarray) -> Dict[str, float]:
    m = {}
    if x.size < 2 or y.size < 2:
        return dict(pearson=np.nan, ks_D=np.nan, emd=np.nan,
            d_mean=np.nan, d_median=np.nan, n=min(x.size, y.size))

    m["pearson"] = float(np.corrcoef(x, y)[0,1])
    ks = sstats.ks_2samp(x, y, alternative="two-sided", method="auto")
    m["ks_D"] = float(ks.statistic)
    m["emd"] = float(sstats.wasserstein_distance(x, y))
    m["d_mean"] = float(np.mean(y) - np.mean(x))
    m["d_median"] = float(np.median(y) - np.median(x))
    m["n"] = int(min(x.size, y.size))
    return m

def _concat_finite(arrs: Iterable[np.ndarray]) -> np.ndarray:
    if not arrs:
        return np.array([], dtype=float)
    xs = []
    for a in arrs:
        if a.size:
            xs.append(a[np.isfinite(a)])
    return np.concatenate(xs) if xs else np.array([], dtype=float)
def _paired_TS_CT_vectors(
    baseline_path: str,
    other_path: str,
    canonical: Dict[Tuple[str,str], Tuple[List[str], List[str]]],
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return (x_TS, y_TS, x_CT, y_CT) where x=*baseline*, y=*other* concatenated over all pairs."""
    x_TS = []; y_TS = []; x_CT = []; y_CT = []
    for (a_n, b_n), (rows, cols) in canonical.items():
        A0 = _fetch_subblock(baseline_path, rows, cols).to_numpy(dtype=float)
        A1 = _fetch_subblock(other_path,    rows, cols).to_numpy(dtype=float)
        if A0.size == 0 or A1.size == 0: 
            continue
        mask = np.isfinite(A0) & np.isfinite(A1)
        x = A0[mask]; y = A1[mask]
        if a_n == b_n:
            x_TS.append(x); y_TS.append(y)
        else:
            x_CT.append(x); y_CT.append(y)
    return (_concat_finite(x_TS), _concat_finite(y_TS),
            _concat_finite(x_CT), _concat_finite(y_CT))

def _within_TS_CT_vectors(
    file_path: str,
    canonical: Dict[Tuple[str,str], Tuple[List[str], List[str]]],
) -> Tuple[np.ndarray, np.ndarray]:
    """Return (TS_vals, CT_vals) within one file, concatenated across pairs."""
    TS = []; CT = []
    for (a_n, b_n), (rows, cols) in canonical.items():
        A = _fetch_subblock(file_path, rows, cols).to_numpy(dtype=float)
        if A.size == 0: 
            continue
        v = A[np.isfinite(A)]
        if a_n == b_n:
            TS.append(v)
        else:
            CT.append(v)
    return _concat_finite(TS), _concat_finite(CT)


@dataclass
class PairResult:
    csv_label: str
    csv_path: str
    region_a: str
    region_b: str
    ct: Optional[int]
    ts: Optional[int]
    stats: Dict[str, float]
    heatmap_png: str
    hist_png: str
    ccdf_lin_png: str          # NEW
    ccdf_log_png: str          # NEW
    sampled_row_ids: List[str]
    sampled_col_names: List[str]



# -------------------------
# Core submatrix extractor
# -------------------------

def plot_region_pair_submatrix(
    csv_path: str,
    region_a: str,
    region_b: str,
    max_rows: int = 300,
    max_cols: int = 300,
    seed: int = 42,
    value_range: Optional[Tuple[float, float]] = None,
    color_scale: str = "Viridis",
    quiet: bool = False,
    ct: Optional[int] = None,
    ts: Optional[int] = None,
    # caches to reduce repeated IO:
    _prefix_to_cols_cache: Optional[Dict[str, List[str]]] = None,
    _row_ids_cache: Optional[Dict[str, List[str]]] = None,
    # NEW:
    fixed_row_ids: Optional[List[str]] = None,
    fixed_col_names: Optional[List[str]] = None,
):
    """
    Stream a large correlation/adjacency CSV and visualize a sampled submatrix between two regions.

    CSV assumptions:
      - The first column is a row-id with region prefix (e.g., "AC_ENSG...").
      - Data columns are genes with region prefix before "_ENSG".
      - Values are numeric (0..1 adjacency/corr OK) and may contain NaNs.

    Returns
    -------
    dict with:
        heatmap_fig, hist_fig,
        stats (enriched),
        sampled_row_ids, sampled_col_names,
        resolved_prefix_a, resolved_prefix_b
    """
    random.seed(seed)

    # --------- scan header only ---------
    lf = pl.scan_csv(csv_path)
    cols = lf.columns
    if not cols:
        raise ValueError("CSV appears to have no columns.")
    row_id_col = cols[0]
    data_cols = [c for c in cols if c != row_id_col]

    # Map region prefix -> columns
    if _prefix_to_cols_cache is None:
        prefix_to_cols: Dict[str, List[str]] = {}
        for c in data_cols:
            pre = _extract_region_prefix(c)
            if pre:
                prefix_to_cols.setdefault(pre, []).append(c)
        _prefix_to_cols_cache = prefix_to_cols
    else:
        prefix_to_cols = _prefix_to_cols_cache

    if not prefix_to_cols:
        raise ValueError("No data columns with '_ENSG' found. Check your CSV schema.")

    # Resolve raw prefixes from user inputs (ignore underscores & case)
    def _resolve_prefix(user_region: str) -> str:
        want = _normalize(user_region)
        cands = [p for p in prefix_to_cols if _normalize(p) == want]
        if not cands:
            available = ", ".join(sorted(prefix_to_cols.keys())[:30])
            raise ValueError(
                f"Region '{user_region}' not found. "
                f"Available examples (first 30): {available} ..."
            )
        if len(cands) > 1 and not quiet:
            print(f"WARNING: multiple raw prefixes matched for '{user_region}': {cands}")
        return cands[0]

    raw_a = _resolve_prefix(region_a)
    raw_b = _resolve_prefix(region_b)
    cols_b_full = prefix_to_cols[raw_b]

    if not quiet:
        print(f"[Header] Region A (rows): {region_a!r} -> {raw_a} | "
              f"Region B (cols): {region_b!r} -> {raw_b} | "
              f"Cols for {raw_b}: {len(cols_b_full)}")

    # --------- collect row IDs for region A ---------
    if _row_ids_cache is not None and raw_a in _row_ids_cache:
        a_rows = _row_ids_cache[raw_a]
    else:
        row_prefix = f"{raw_a}_"
        a_rows = (
            lf.select(row_id_col)
              .filter(pl.col(row_id_col).str.starts_with(row_prefix))
              .collect(streaming=True)[row_id_col]
              .to_list()
        )
        if _row_ids_cache is not None:
            _row_ids_cache[raw_a] = a_rows

    if not a_rows:
        raise ValueError(f"No rows found with prefix '{raw_a}_' in {row_id_col}.")

    # --------- sampling (canonical if provided) ---------
    if fixed_row_ids is not None:
        avail = set(a_rows)
        a_rows_sampled = [r for r in fixed_row_ids if r in avail]
    else:
        if max_rows and len(a_rows) > max_rows:
            a_rows_sampled = sorted(random.sample(a_rows, max_rows))
        else:
            a_rows_sampled = a_rows

    if fixed_col_names is not None:
        availc = set(cols_b_full)
        b_cols_sampled = [c for c in fixed_col_names if c in availc]
    else:
        if max_cols and len(cols_b_full) > max_cols:
            b_cols_sampled = sorted(random.sample(cols_b_full, max_cols))
        else:
            b_cols_sampled = cols_b_full


    if not quiet:
        print(f"[Sample] Rows A: {len(a_rows_sampled)} / {len(a_rows)}  |  "
              f"Cols B: {len(b_cols_sampled)} / {len(cols_b_full)}")

    # --------- fetch submatrix lazily ---------
    lf_sub = (
        lf.select([row_id_col] + b_cols_sampled)
          .filter(pl.col(row_id_col).is_in(a_rows_sampled))
    )
    sub_df = lf_sub.collect(streaming=True).to_pandas()
    sub_df = sub_df.set_index(row_id_col).loc[a_rows_sampled]  # stable row order
    if fixed_row_ids is not None:
        sub_df = sub_df.reindex(index=a_rows_sampled)
    if fixed_col_names is not None:
        sub_df = sub_df.reindex(columns=b_cols_sampled)
    # --------- compute stats ---------
    sub_np = sub_df.to_numpy(dtype=float)
    finite_mask = np.isfinite(sub_np)
    finite_vals = sub_np[finite_mask]

    # Degrees (treat NaN as 0)
    z = sub_np.copy()
    z[~np.isfinite(z)] = 0.0
    row_k = np.sum(z, axis=1)
    col_k = np.sum(z, axis=0)

    def _pct(x: np.ndarray, q: float) -> float:
        return float(np.nanpercentile(x, q)) if x.size else float("nan")

    strong_thresholds = (0.2, 0.5, 0.8)
    def _prop_ge(x: np.ndarray, thr: float) -> float:
        if x.size == 0:
            return float("nan")
        return float(np.mean(x >= thr))

    stats = {
        "n_rows": int(sub_df.shape[0]),
        "n_cols": int(sub_df.shape[1]),
        "nnz_finite": int(finite_vals.size),
        "nan_count": int(np.size(sub_np) - finite_vals.size),
        "sparsity_prop_zero": float(np.mean(sub_np == 0.0)) if sub_np.size else float("nan"),

        # value distribution
        "val_mean": float(np.nanmean(finite_vals)) if finite_vals.size else float("nan"),
        "val_std":  float(np.nanstd(finite_vals))  if finite_vals.size else float("nan"),
        "val_min":  float(np.nanmin(finite_vals))  if finite_vals.size else float("nan"),
        "val_p1":   _pct(finite_vals, 1),
        "val_p5":   _pct(finite_vals, 5),
        "val_p25":  _pct(finite_vals, 25),
        "val_median": _pct(finite_vals, 50),
        "val_p75":  _pct(finite_vals, 75),
        "val_p95":  _pct(finite_vals, 95),
        "val_p99":  _pct(finite_vals, 99),
        "val_max":  float(np.nanmax(finite_vals))  if finite_vals.size else float("nan"),
        "val_skew": float(sstats.skew(finite_vals, nan_policy="omit")) if finite_vals.size else float("nan"),
        "val_kurtosis": float(sstats.kurtosis(finite_vals, nan_policy="omit")) if finite_vals.size else float("nan"),
        "finite_vals": finite_vals,

        # degree summaries
        "row_k_sum": float(np.sum(row_k)),
        "row_k_mean": float(np.mean(row_k)),
        "row_k_median": float(np.median(row_k)),
        "col_k_sum": float(np.sum(col_k)),
        "col_k_mean": float(np.mean(col_k)),
        "col_k_median": float(np.median(col_k)),

        # strong edges (over finite values only)
        "prop_ge_0.2": _prop_ge(finite_vals, strong_thresholds[0]),
        "prop_ge_0.5": _prop_ge(finite_vals, strong_thresholds[1]),
        "prop_ge_0.8": _prop_ge(finite_vals, strong_thresholds[2]),
    }

    # Top-degree genes (IDs)
    k = min(5, sub_df.shape[0])  # top 5 max
    top_rows_idx = np.argsort(-row_k)[:k] if sub_df.shape[0] else []
    top_cols_idx = np.argsort(-col_k)[:k] if sub_df.shape[1] else []
    stats["top_rows_by_degree"] = ", ".join([sub_df.index[i] for i in top_rows_idx]) if k else ""
    stats["top_cols_by_degree"] = ", ".join([sub_df.columns[i] for i in top_cols_idx]) if k else ""

    # --------- figures ---------
    ct_lab = f"{ct}" if ct is not None else "?"
    ts_lab = f"{ts}" if ts is not None else "?"
    title_main = f"Adjacency: {region_a} (rows) × {region_b} (cols) | CT={ct_lab}, TS={ts_lab}"
    s = stats
    subtitle = (
        f"Rows k (A→B): sum={s['row_k_sum']:.3f}, mean={s['row_k_mean']:.3f}, median={s['row_k_median']:.3f} | "
        f"Cols k (B→A): sum={s['col_k_sum']:.3f}, mean={s['col_k_mean']:.3f}, median={s['col_k_median']:.3f} | "
        f"Finite cells: {finite_vals.size:,} / {sub_np.size:,} | "
        f"Strong edges ≥0.5: {s.get('prop_ge_0.5', float('nan')):.2%}"
    )

    heatmap_fig = px.imshow(
        sub_np,
        x=sub_df.columns,
        y=sub_df.index,
        color_continuous_scale=color_scale,
        zmin=value_range[0] if value_range else None,
        zmax=value_range[1] if value_range else None,
        aspect="auto",
    )
    heatmap_fig.update_layout(
        title=f"{title_main}<br><sup>{subtitle}</sup>",
        xaxis_title=f"{region_b} genes (sampled {sub_df.shape[1]})",
        yaxis_title=f"{region_a} genes (sampled {sub_df.shape[0]})",
        width=1100,
        height=880,
        coloraxis_colorbar=dict(title="Adj."),
        margin=dict(l=80, r=20, t=100, b=80),
    )

    hist_fig = px.histogram(
        x=finite_vals,
        nbins=60,
        title=f"Value distribution: {region_a}×{region_b} (finite cells {finite_vals.size:,})",
    )
    hist_fig.update_layout(
        xaxis_title="Correlation / adjacency",
        yaxis_title="Count",
        bargap=0.02,
        width=1000,
        height=460,
        margin=dict(l=60, r=20, t=70, b=60),
        annotations=[
            dict(
                x=0.01, y=1.10, xref="paper", yref="paper", xanchor="left",
                text=(f"μ={stats['val_mean']:.3f}, σ={stats['val_std']:.3f}, "
                      f"median={stats['val_median']:.3f}, "
                      f"[p5,p95]=[{stats['val_p5']:.3f},{stats['val_p95']:.3f}]"),
                showarrow=False, font=dict(size=11)
            )
        ]
    )

    return {
        "heatmap_fig": heatmap_fig,
        "hist_fig": hist_fig,
        "stats": stats,
        "sampled_row_ids": a_rows_sampled,
        "sampled_col_names": b_cols_sampled,
        "resolved_prefix_a": raw_a,
        "resolved_prefix_b": raw_b,
    }


# -------------------------
# Report builder
# -------------------------
def make_region_pairs_pdf_report(
    csv_paths: Sequence[str],
    out_pdf_path: str,
    pairs: Optional[Sequence[Tuple[str, str]]] = None,
    max_rows: int = 300,
    max_cols: int = 300,
    seed: int = 42,
    global_value_range: Optional[Tuple[float, float]] = None,
    color_scale: str = "Viridis",
    quiet: bool = False,
    ct_ts_overrides: Optional[dict] = None,
    out_html_path: Optional[str] = None,
    # ==== NEW ====
    file_display_names: Optional[Dict[str, str]] = None,
    overlay_max_pairs_html: int = 24,
) -> List[PairResult]:


    """
    Build a comprehensive multi-page PDF covering requested (or all) region pairs
    across one or more CSV adjacency/correlation files.

    - Each CSV gets a summary table page.
    - Each pair gets two figures (heatmap + histogram) and an expanded stats panel.
    - Returns a list of PairResult (also persisted in the PDF).

    NOTE: Uses Plotly+Kaleido to create PNGs, and ReportLab to compose the PDF.
    """
    os.makedirs(os.path.dirname(out_pdf_path) or ".", exist_ok=True)
    rng = random.Random(seed)

    # Prepare PDF doc
    page_w, page_h = A4
    left_margin = right_margin = top_margin = bottom_margin = 36
    max_img_width = page_w - left_margin - right_margin

    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle(name="Small", fontSize=9, leading=11))
    title_style = styles["Title"]
    h1 = styles["Heading1"]
    h2 = styles["Heading2"]
    body = styles["BodyText"]
    small = styles["Small"]

    story = []
    results: List[PairResult] = []

    html_sections: List[str] = []

    # ... Overview copy (more explicit)  ==== UPDATED ====
    html_sections.append(f"""
<section class="section">
  <h2>Overview</h2>
  <p>This report compiles region-pair submatrices from <b>{len(csv_paths)}</b> adjacency/correlation file(s).</p>
  <p>For each file we show: (1) tissue-pair heatmaps and K-distributions, (2) overlays of all tissues, (3) CCDFs (linear & log–log),
     (4) power-law tail fits (α, xmin, KS-D), (5) TS vs CT within-file summaries, and (6) compact delta tables.</p>
  <p class="muted">Color scale: {"global ["+str(global_value_range[0])+", "+str(global_value_range[1])+"]" if global_value_range else "auto per submatrix"}.
  KS p-values and scatter plots were intentionally omitted for clarity.</p>
</section>
<div class="divider"></div>
""")
    report_title = "Adjacency Matrix Region-Pair Report"
    # === label each file (baseline name overridden) ==== UPDATED ====
    baseline_path = csv_paths[0]
    baseline_label = _display_name_for_csv(baseline_path, file_display_names, baseline=True)
    # Cover page
    story.append(Paragraph("Dynamic vs Constant Beta (values)", title_style))

    story.append(Spacer(1, 0.2 * inch))
    story.append(Paragraph(
        f"This report compiles region-pair submatrices from {len(csv_paths)} adjacency matrix(es) into a single document. "
        f"For each pair, we show a heatmap of top {max_rows} rows and {max_cols} cols and the correlations distribution of 90,000 values, "
        f"along with descriptive statistics, degree summaries, sparsity and strong-edge proportions.",
        body
    ))
    story.append(Spacer(1, 0.1 * inch))

    story.append(PageBreak())
    # Decide global pairs (must exist in ALL files so sampling is identical)
    # If user gave pairs -> use those; else intersect prefixes across files
    if pairs is None:
        # get intersection of normalized prefixes
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

    with tempfile.TemporaryDirectory() as tmpdir:
        for csv_path in csv_paths:
            lf = pl.scan_csv(csv_path)
            cols = lf.columns
            if not cols:
                raise ValueError(f"{csv_path} has no columns.")
            row_id_col = cols[0]
            data_cols = [c for c in cols if c != row_id_col]

            # Build prefix-to-cols map and discover all region prefixes
            prefix_to_cols: Dict[str, List[str]] = {}
            for c in data_cols:
                pre = _extract_region_prefix(c)
                if pre:
                    prefix_to_cols.setdefault(pre, []).append(c)
            prefixes = sorted(prefix_to_cols.keys())
            if not prefixes:
                raise ValueError(f"{csv_path}: no columns matched '*_ENSG...'")

            csv_label = _display_name_for_csv(csv_path, file_display_names, baseline=(csv_path == baseline_path))
            ct, ts = _parse_ct_ts_from_path(csv_path)

            # Decide which pairs to run
            if pairs is not None:
                requested_pairs = pairs
            else:
                # all pairs including within-tissue (diagonal)
                requested_pairs = [(a, b) for i, a in enumerate(prefixes) for b in prefixes[i:]]  # upper triangle including diagonal
            
            # Use global requested pairs that survived canonical sampling
            requested_pairs = []
            for (a_n, b_n) in requested_pairs_global:
                if (a_n, b_n) in canonical:
                    # map back to this file's raw labels for printing
                    # but keep normalized left/right names for keys
                    requested_pairs.append((a_n, b_n))

            # Fast row ids cache per region (avoid re-scanning)
            row_ids_cache: Dict[str, List[str]] = {}

            # Summary table rows for this CSV
            table_rows = [["Region A", "Region B", "Rows", "Cols", "μ", "σ", "Median", "p5", "p95", "≥0.5"]]

            # CSV section header
            story.append(Paragraph(csv_label, h1))
            story.append(Spacer(1, 0.08 * inch))
            story.append(Paragraph(
                f"Detected regions (prefixes): {', '.join(prefixes)}",
                small
            ))
            story.append(Spacer(1, 0.12 * inch))
            html_sections.append(f"""
<section class="section" id="{_html_escape(csv_label)}">
  <h2>{_html_escape(csv_label)}</h2>
  <p class="muted">Detected regions: {" ".join(f'<span class="chip">{_html_escape(p)}</span>' for p in prefixes)}</p>
""")
            

            #for (ra, rb) in requested_pairs:
            for (a_n, b_n) in requested_pairs:
                ra, rb = a_n, b_n  # keep normalized for titles; display is fine
                # resolve CT/TS for this specific pair (ra, rb)
                pair_ct, pair_ts = _get_ct_ts_for_pair(
                    csv_path=csv_path,
                    ra=ra,
                    rb=rb,
                    overrides=ct_ts_overrides,
                    default_ct=ct,
                    default_ts=ts,
                )
                fixed_rows, fixed_cols = canonical[(a_n, b_n)]

                try:
                    out = plot_region_pair_submatrix(
                        csv_path=csv_path,
                        region_a=ra,
                        region_b=rb,
                        max_rows=max_rows,
                        max_cols=max_cols,
                        seed=rng.randint(0, 2**31-1),
                        value_range=global_value_range,
                        color_scale=color_scale,
                        quiet=quiet,
                        ct=pair_ct,    # <— use per-pair CT/TS
                        ts=pair_ts,
                        _prefix_to_cols_cache=prefix_to_cols,
                        _row_ids_cache=row_ids_cache,
                        fixed_row_ids=fixed_rows,           # NEW
                        fixed_col_names=fixed_cols,         # NEW
                    )

                except Exception as e:
                    story.append(Paragraph(f"<b>{ra} × {rb}</b> — ERROR: {e}", body))
                    story.append(Spacer(1, 0.05 * inch))
                    continue
                heatmap_html = _wrap_fig_html(out["heatmap_fig"], title=f"{ra} × {rb} — Heatmap")
                hist_html    = _wrap_fig_html(out["hist_fig"],    title=f"{ra} × {rb} — Distribution")
                # ---- (3) power-law tail fit + CCDF overlays for this pair ----
                # Re-extract the exact finite K values for this pair using the canonical rows/cols
                finite_vals = out['stats']["finite_vals"]

                # Fit a continuous power-law tail (α, xmin) and KS-D on the tail only
                fit = _fit_powerlaw(finite_vals, xmin_grid=None, min_tail=100)

                # Build CCDF figures (linear, log-log) and overlay the theoretical power-law
                ccdf_lin, ccdf_log = _ccdf_figs(finite_vals, f"{ra}×{rb}")
                _add_powerlaw_to_ccdf(ccdf_lin, fit["alpha"], fit["xmin"], "Power-law fit")
                _add_powerlaw_to_ccdf(ccdf_log, fit["alpha"], fit["xmin"], "Power-law fit")
                ccdf_lin_png = os.path.join(tmpdir, f"{os.path.basename(csv_path)}_{ra}_x_{rb}_ccdf_lin.png")
                ccdf_log_png = os.path.join(tmpdir, f"{os.path.basename(csv_path)}_{ra}_x_{rb}_ccdf_log.png")
                ccdf_lin.write_image(ccdf_lin_png, scale=2)
                ccdf_log.write_image(ccdf_log_png, scale=2)

                # A small tail-fit summary (no KS p-value)
                tail_html = (
                    f"<div class='desc'>Tail fit: "
                    f"α={fit['alpha']:.3f} (xmin={fit['xmin']:.4f}), "
                    f"KS-D={fit['ks_D']:.4f}, n_tail={fit['n_tail']}</div>"
                )

                # Save figs to temp PNGs
                hm_png = os.path.join(tmpdir, f"{os.path.basename(csv_path)}_{ra}_x_{rb}_hm.png")
                hist_png = os.path.join(tmpdir, f"{os.path.basename(csv_path)}_{ra}_x_{rb}_hist.png")
                out["heatmap_fig"].write_image(hm_png, scale=2)
                out["hist_fig"].write_image(hist_png, scale=2)

                # Append summary row
                s = out["stats"]
                table_rows.append([
                    ra, rb,
                    f"{s['n_rows']}", f"{s['n_cols']}",
                    f"{s['val_mean']:.3f}", f"{s['val_std']:.3f}",
                    f"{s['val_median']:.3f}",
                    f"{s['val_p5']:.3f}", f"{s['val_p95']:.3f}",
                    f"{s['prop_ge_0.5']:.2%}" if not math.isnan(s['prop_ge_0.5']) else "NA",
                ])

                # results.append(PairResult(
                #     csv_label=csv_label, csv_path=csv_path,   # <— add csv_path
                #     region_a=ra, region_b=rb,
                #     ct=ct, ts=ts, stats=s,
                #     heatmap_png=hm_png, hist_png=hist_png,
                #     sampled_row_ids=out["sampled_row_ids"],
                #     sampled_col_names=out["sampled_col_names"],
                # ))
                s_small = {k: v for k, v in out["stats"].items() if k != "finite_vals"}

                results.append(PairResult(
                    csv_label=csv_label, csv_path=csv_path,
                    region_a=ra, region_b=rb,
                    ct=pair_ct, ts=pair_ts,
                    stats=s_small,                      # use s_small if you applied the memory tip
                    heatmap_png=hm_png, hist_png=hist_png,
                    ccdf_lin_png=ccdf_lin_png,          # NEW
                    ccdf_log_png=ccdf_log_png,          # NEW
                    sampled_row_ids=out["sampled_row_ids"],
                    sampled_col_names=out["sampled_col_names"],
                ))


                desc_html = (
                    f"<div class='desc'>"
                    f"<b>Rows</b>={s['n_rows']}, <b>Cols</b>={s['n_cols']}; "
                    f"<b>finite</b>={s['nnz_finite']:,}, <b>NaN</b>={s['nan_count']:,}; "
                    f"<b>μ</b>={s['val_mean']:.3f}, <b>σ</b>={s['val_std']:.3f}, "
                    f"<b>median</b>={s['val_median']:.3f}, "
                    f"<b>[p5,p95]</b>=[{s['val_p5']:.3f},{s['val_p95']:.3f}], "
                    f"<b>min</b>={s['val_min']:.3f}, <b>max</b>={s['val_max']:.3f}; "
                    f"<b>skew</b>={s['val_skew']:.3f}, <b>kurtosis</b>={s['val_kurtosis']:.3f}. "
                    f"<br/>A→B degree: sum={s['row_k_sum']:.3f}, μ={s['row_k_mean']:.3f}, median={s['row_k_median']:.3f}; "
                    f"B→A degree: sum={s['col_k_sum']:.3f}, μ={s['col_k_mean']:.3f}, median={s['col_k_median']:.3f}. "
                    f"<br/>Strong edges: ≥0.2={s['prop_ge_0.2']:.2%}, ≥0.5={s['prop_ge_0.5']:.2%}, ≥0.8={s['prop_ge_0.8']:.2%}."
                    f"</div>"
                )
                html_sections.append(f"<h3>{_html_escape(ra)} × {_html_escape(rb)}</h3>")
                html_sections.append(desc_html)
                html_sections.append(tail_html)
                html_sections.append(_wrap_fig_html(ccdf_lin, title=f"{ra}×{rb} — CCDF + power-law"))
                html_sections.append(_wrap_fig_html(ccdf_log, title=f"{ra}×{rb} — CCDF (log–log) + power-law"))

                html_sections.append(heatmap_html)
                html_sections.append(hist_html)


            # Add summary table for this CSV
            tbl = Table(table_rows, hAlign="LEFT")
            tbl.setStyle(TableStyle([
                ("FONTNAME", (0,0), (-1,0), "Helvetica-Bold"),
                ("BACKGROUND", (0,0), (-1,0), colors.lightgrey),
                ("GRID", (0,0), (-1,-1), 0.3, colors.grey),
                ("ALIGN", (2,1), (-1,-1), "RIGHT"),
            ]))
            story.append(tbl)
            story.append(Spacer(1, 0.15 * inch))
            html_sections.append(_table_html(table_rows, caption="Per-pair summary (this file)"))

                        # ---- (4) WITHIN-FILE compacts & overlays ------------------------------------
            # Build a tidy frame of K summaries per tissue pair in this file
            per_pair_records = []
            for pr in [r for r in results if r.csv_label == csv_label]:
                s = pr.stats
                per_pair_records.append(dict(
                    pair=f"{pr.region_a}|{pr.region_b}",
                    A=pr.region_a, B=pr.region_b,
                    mean=s["val_mean"], median=s["val_median"], p_ge_05=s["prop_ge_0.5"],
                    n=s["nnz_finite"]
                ))
            df_pairs = pd.DataFrame(per_pair_records)

            # (4a) Table: K summaries per tissue pair
            html_sections.append(_h3("Within-file: tissue pair summaries"))
            html_sections.append(_table_html(
                [["Pair","K mean","L median","prop≥0.5","n"]] + 
                [[r["pair"], f"{r['mean']:.4f}", f"{r['median']:.4f}", f"{r['p_ge_05']:.3f}", f"{r['n']:,}"] for r in per_pair_records],
                caption="K summaries per tissue pair"
            ))

            # (4b) Tables: Pair-of-pairs ΔK
            rows_delta = [["Comparison", "Δ mean K", "Δ median K", "Δ prop≥0.5"]]
            if not df_pairs.empty:
                idx = {(r.A, r.B): r for r in df_pairs.itertuples(index=False)}
                speciesA = sorted(set(df_pairs.A))
                speciesB = sorted(set(df_pairs.B))
                # same A: (A|B) vs (A|C)
                for A in speciesA:
                    Bs = [b for b in speciesB if (A,b) in idx]
                    for i in range(len(Bs)):
                        for j in range(i+1, len(Bs)):
                            r1 = idx[(A, Bs[i])]; r2 = idx[(A, Bs[j])]
                            rows_delta.append([f"{A}|{Bs[i]}  vs  {A}|{Bs[j]}",
                                f"{(r1.mean - r2.mean):+.4f}", f"{(r1.median - r2.median):+.4f}", f"{(r1.p_ge_05 - r2.p_ge_05):+.3f}"])
                # same B: (C|B) vs (A|B)
                for B in speciesB:
                    As = [a for a in speciesA if (a,B) in idx]
                    for i in range(len(As)):
                        for j in range(i+1, len(As)):
                            r1 = idx[(As[i], B)]; r2 = idx[(As[j], B)]
                            rows_delta.append([f"{As[i]}|{B}  vs  {As[j]}|{B}",
                                f"{(r1.mean - r2.mean):+.4f}", f"{(r1.median - r2.median):+.4f}", f"{(r1.p_ge_05 - r2.p_ge_05):+.3f}"])
                # diagonal-style: (C|B) vs (A|C)
                for A in speciesA:
                    for B in speciesB:
                        for C in speciesA:
                            if (C,B) in idx and (A,C) in idx and not (A==C and B==C):
                                r1 = idx[(C,B)]; r2 = idx[(A,C)]
                                rows_delta.append([f"{C}|{B}  vs  {A}|{C}",
                                    f"{(r1.mean - r2.mean):+.4f}", f"{(r1.median - r2.median):+.4f}", f"{(r1.p_ge_05 - r2.p_ge_05):+.3f}"])

            if len(rows_delta) > 1:
                html_sections.append(_h3("Within-file: pair-of-pairs ΔK"))
                html_sections.append(_table_html(rows_delta, caption="All relevant pair-of-pairs comparisons"))

            # (4c) Overlays: all tissue-pair histograms + CCDFs (cap number of traces)
            show_pairs = per_pair_records[:overlay_max_pairs_html]
            overlay_df = []
            pair_vals_cache: Dict[str, np.ndarray] = {}
            for r in show_pairs:
                pr = next(x for x in results if x.csv_label==csv_label and f"{x.region_a}|{x.region_b}" == r["pair"])
                Ablock = _fetch_subblock(pr.csv_path, pr.sampled_row_ids, pr.sampled_col_names).to_numpy(float)
                v = Ablock[np.isfinite(Ablock)].ravel()
                pair_vals_cache[r["pair"]] = v
                overlay_df.append(pd.DataFrame({"value": v, "pair": r["pair"]}))

            if overlay_df:
                df_overlay = pd.concat(overlay_df, ignore_index=True)
                figH = px.histogram(df_overlay, x="value", color="pair", barmode="overlay", nbins=60, opacity=0.45,
                                    title=f"All tissue pairs — histogram (max {overlay_max_pairs_html} traces)")
                figH.update_layout(xaxis_title="K (value)", yaxis_title="Count", width=1000, height=520, legend=dict(orientation="h"))
                html_sections.append(_wrap_fig_html(figH, title=figH.layout.title.text))

                import plotly.graph_objects as go
                figC = go.Figure(); figClog = go.Figure()
                for r in show_pairs:
                    xs, cc = _ccdf(pair_vals_cache[r["pair"]])
                    figC.add_trace(go.Scatter(x=xs, y=cc, mode="lines", name=r["pair"]))
                    m = xs > 0
                    figClog.add_trace(go.Scatter(x=xs[m], y=cc[m], mode="lines", name=r["pair"]))
                figC.update_layout(title="All tissue pairs — CCDF", xaxis_title="K", yaxis_title="P(X ≥ x)",
                                width=1000, height=520, legend=dict(orientation="h"))
                figClog.update_layout(title="All tissue pairs — CCDF (log–log)", xaxis_type="log", yaxis_type="log",
                                    xaxis_title="K (log)", yaxis_title="P(X ≥ x) (log)",
                                    width=1000, height=520, legend=dict(orientation="h"))
                html_sections.append(_wrap_fig_html(figC, title=figC.layout.title.text))
                html_sections.append(_wrap_fig_html(figClog, title=figClog.layout.title.text))

            # (4d) TS vs CT in this file (CCDFs)
            TS_vals, CT_vals = _within_TS_CT_vectors(csv_path, canonical)
            for fig in (_ccdf_figs(TS_vals, f"{csv_label}: TS") + _ccdf_figs(CT_vals, f"{csv_label}: CT")):
                html_sections.append(_wrap_fig_html(fig, title=fig.layout.title.text))
            html_sections.append("</section><div class='divider'></div>")
            # ---- end (4) ---------------------------------------------------------------


            # Add pages per pair
            for pr in [r for r in results if r.csv_label == csv_label]:
                story.append(Paragraph(f"{pr.region_a} × {pr.region_b}", h2))
                s = pr.stats
                # concise description block
                desc = (
                    f"Rows={s['n_rows']}, Cols={s['n_cols']}; "
                    f"finite={s['nnz_finite']:,}, NaN={s['nan_count']:,}, "
                    f"Value μ={s['val_mean']:.3f}, σ={s['val_std']:.3f}, "
                    f"median={s['val_median']:.3f}, "
                    f"[p5,p95]=[{s['val_p5']:.3f},{s['val_p95']:.3f}], "
                    f"min={s['val_min']:.3f}, max={s['val_max']:.3f}, "
                    f"skew={s['val_skew']:.3f}, kurtosis={s['val_kurtosis']:.3f}<br/>"
                    f"\n"
                    f"A→B degree: sum={s['row_k_sum']:.3f}, μ={s['row_k_mean']:.3f}, median={s['row_k_median']:.3f}; "
                    f"\n"
                    f"B→A degree: sum={s['col_k_sum']:.3f}, μ={s['col_k_mean']:.3f}, median={s['col_k_median']:.3f}<br/>"
                    f"\n"
                    f"Strong edges: ≥0.2={s['prop_ge_0.2']:.2%}, ≥0.5={s['prop_ge_0.5']:.2%}, ≥0.8={s['prop_ge_0.8']:.2%}<br/>"
                )
                story.append(Paragraph(desc, small))
                story.append(Spacer(1, 0.08 * inch))

                # figures
                story.append(_image_flowable(pr.heatmap_png, max_img_width))
                story.append(Spacer(1, 0.08 * inch))
                story.append(_image_flowable(pr.hist_png, max_img_width))
                story.append(Spacer(1, 0.06 * inch))
                story.append(_image_flowable(pr.ccdf_lin_png, max_img_width))  # NEW
                story.append(Spacer(1, 0.06 * inch))
                story.append(_image_flowable(pr.ccdf_log_png, max_img_width))  # NEW
                story.append(Spacer(1, 0.06 * inch))


                story.append(PageBreak())
        
        
        if len(csv_paths) >= 2:
            story.append(Paragraph("Cross-file comparison", h1))
            story.append(Spacer(1, 0.12 * inch))

            # Map: label -> results by pair key
            # Keep orientation as in stored results; handle reversed if needed.
            by_label: Dict[str, Dict[Tuple[str,str], PairResult]] = {}
            for r in results:
                pk = _pair_key(r.region_a, r.region_b)
                by_label.setdefault(r.csv_label, {})[pk] = r

            baseline_path = csv_paths[0]
            baseline_label = _display_name_for_csv(baseline_path, file_display_names, baseline=True)
            story.append(Paragraph(f"Baseline file: <b>{baseline_label}</b>", body))
            story.append(Spacer(1, 0.08 * inch))
            ks_global: List[Dict[str, float]] = []  # <-- move here, outside the loop

            # For each other file, compare per-pair vs baseline
            for other_path in csv_paths[1:]:
                other_label = _display_name_for_csv(other_path, file_display_names, baseline=False)

                story.append(Paragraph(f"Compared file: {other_label}", h2))
                story.append(Spacer(1, 0.06 * inch))

                base_pairs = by_label.get(baseline_label, {})
                other_pairs = by_label.get(other_label, {})

                # If orientation differs (B×A in other), we'll try reversed key on demand
                def _get_other_pair(ra_norm: str, rb_norm: str) -> Tuple[Optional[PairResult], bool]:
                    # exact orientation
                    p = other_pairs.get((ra_norm, rb_norm))
                    if p is not None:
                        return p, False
                    # reversed (needs transpose)
                    p_rev = other_pairs.get((rb_norm, ra_norm))
                    if p_rev is not None:
                        return p_rev, True
                    return None, False

                # Collect summary rows
                comp_rows = [["Pair", "Aligned cells", "Pearson r", "KS D", "Δ mean", "Δ median", "Δ prop≥0.5"]]
                ks_records: List[Dict[str, float]] = []


                for (ra_n, rb_n), base_pr in base_pairs.items():
                    other_pr, need_transpose = _get_other_pair(ra_n, rb_n)
                    if other_pr is None:
                        continue

                    # Align on intersection of sampled ids
                    base_rows = base_pr.sampled_row_ids
                    base_cols = base_pr.sampled_col_names
                    # if reversed, swap row/col lists from 'other'
                    other_rows = other_pr.sampled_row_ids if not need_transpose else other_pr.sampled_col_names
                    other_cols = other_pr.sampled_col_names if not need_transpose else other_pr.sampled_row_ids

                    rows_inter = [r for r in base_rows if r in set(other_rows)]
                    cols_inter = [c for c in base_cols if c in set(other_cols)]
                    if len(rows_inter) == 0 or len(cols_inter) == 0:
                        continue

                    # Fetch aligned blocks
                    A0 = _fetch_subblock(baseline_path, rows_inter, cols_inter).to_numpy(dtype=float)
                    A1 = _fetch_subblock(other_path,    rows_inter, cols_inter).to_numpy(dtype=float)
                    if A0.size == 0 or A1.size == 0:
                        continue

                    mask = np.isfinite(A0) & np.isfinite(A1)
                    x = A0[mask]; y = A1[mask]
                    n_aligned = int(x.size)
                    if n_aligned < 3:
                        continue

                    # Metrics
                    pearson = float(np.corrcoef(x, y)[0, 1]) if n_aligned >= 2 else float("nan")
                    ks_res = sstats.ks_2samp(x, y, alternative="two-sided", method="auto") if n_aligned >= 2 else None
                    ks_D = float(ks_res.statistic) if ks_res else float("nan")
                    ks_p = float(ks_res.pvalue) if ks_res else float("nan")
                    d_mean = float(np.mean(y) - np.mean(x))
                    d_median = float(np.median(y) - np.median(x))
                    d_prop05 = float(np.mean(y >= 0.5) - np.mean(x >= 0.5))

                    comp_rows.append([
                        f"{base_pr.region_a}×{base_pr.region_b}",
                        f"{n_aligned:,}",
                        f"{pearson:.3f}",
                        f"{ks_D:.3f}",
                        f"{d_mean:+.4f}",
                        f"{d_median:+.4f}",
                        f"{d_prop05:+.3f}",
                    ])


                    ks_records.append({
                        "pair": f"{base_pr.region_a}×{base_pr.region_b}",
                        "ks_D": ks_D,
                        "ks_p": ks_p,
                        "n_aligned": n_aligned
                    })


                    # Plots: overlay histogram + Δ heatmap
                    # Histogram overlay
                    df_overlay = pd.DataFrame({
                        "value": np.concatenate([x, y]),
                        "which": ([baseline_label] * x.size) + ([other_label] * y.size),
                    })
                    hist_overlay = px.histogram(
                        df_overlay, x="value", color="which",
                        barmode="overlay", nbins=60, opacity=0.6,
                        title=f"Distribution comparison: {base_pr.region_a}×{base_pr.region_b}"
                    )
                    hist_overlay.update_layout(
                        xaxis_title="Adjacency / correlation",
                        yaxis_title="Count",
                        width=900, height=420, margin=dict(l=50, r=20, t=60, b=50),
                    )
                    hist_png = os.path.join(tmpdir, f"overlay_{baseline_label}_{other_label}_{base_pr.region_a}x{base_pr.region_b}.png")
                    hist_overlay.write_image(hist_png, scale=2)
                    
                    if out_html_path:
                        html_sections.append(_wrap_fig_html(
                            hist_overlay, title=f"Distribution comparison: {base_pr.region_a}×{base_pr.region_b} ({_html_escape(other_label)} vs baseline)"
                        ))

                    story.append(_image_flowable(hist_png, max_img_width))
                    story.append(Spacer(1, 0.05 * inch))

                    # Δ heatmap (other - baseline), center at 0 with symmetric limits
                    D = A1 - A0
                    if D.size >= 100:  # only show if reasonably sized
                        vmax = float(np.nanmax(np.abs(D))) if np.isfinite(D).any() else 1.0
                        delta_fig = px.imshow(
                            D, color_continuous_scale="RdBu", zmin=-vmax, zmax=+vmax,
                            title=f"Δ heatmap (other − baseline): {base_pr.region_a}×{base_pr.region_b}"
                        )
                        delta_fig.update_layout(
                            width=900, height=620, margin=dict(l=60, r=20, t=60, b=60),
                            coloraxis_colorbar=dict(title="Δ")
                        )
                        delta_png = os.path.join(tmpdir, f"delta_{baseline_label}_{other_label}_{base_pr.region_a}x{base_pr.region_b}.png")
                        delta_fig.write_image(delta_png, scale=2)
                        story.append(_image_flowable(delta_png, max_img_width))
                        story.append(Spacer(1, 0.08 * inch))
                        if out_html_path:
                            html_sections.append(_wrap_fig_html(
                                delta_fig, title=f"Δ heatmap (other − baseline): {base_pr.region_a}×{base_pr.region_b}"
                            ))
                # After finishing the for (ra_n, rb_n) ... loop:
                for rec in ks_records:
                    ks_global.append({"file": other_label, "pair": rec["pair"], "ks_D": rec["ks_D"]})


                # Per-file summary table
                tbl = Table(comp_rows, hAlign="LEFT")
                tbl.setStyle(TableStyle([
                    ("FONTNAME", (0,0), (-1,0), "Helvetica-Bold"),
                    ("BACKGROUND", (0,0), (-1,0), colors.HexColor("#EFEFEF")),
                    ("GRID", (0,0), (-1,-1), 0.3, colors.grey),
                    ("ALIGN", (1,1), (-1,-1), "RIGHT"),
                ]))
                story.append(tbl)
                story.append(Spacer(1, 0.15 * inch))
                # ===== TS/CT similarity preservation & convergence =====
                # Cross-file: TS vs baseline, CT vs baseline (preservation)
                x_TS, y_TS, x_CT, y_CT = _paired_TS_CT_vectors(
                    baseline_path=baseline_path,
                    other_path=other_path,
                    canonical=canonical,
                )
                m_TS = _metrics_pair(x_TS, y_TS)
                m_CT = _metrics_pair(x_CT, y_CT)

                # Within-file: TS vs CT (convergence) for baseline and other
                base_TS, base_CT = _within_TS_CT_vectors(baseline_path, canonical)
                oth_TS,  oth_CT  = _within_TS_CT_vectors(other_path,    canonical)

                w_base = _metrics_pair(base_TS, base_CT)
                w_oth  = _metrics_pair(oth_TS,  oth_CT)

                conv_rows = [
                    ["Scope", "n", "Pearson r", "KS D", "EMD", "Δ mean", "Δ median"],
                    [f"TS vs baseline ({other_label})", f"{m_TS['n']:,}", f"{m_TS['pearson']:.3f}", f"{m_TS['ks_D']:.3f}", f"{m_TS['emd']:.4f}", f"{m_TS['d_mean']:+.4f}", f"{m_TS['d_median']:+.4f}"],
                    [f"CT vs baseline ({other_label})", f"{m_CT['n']:,}", f"{m_CT['pearson']:.3f}", f"{m_CT['ks_D']:.3f}", f"{m_CT['emd']:.4f}", f"{m_CT['d_mean']:+.4f}", f"{m_CT['d_median']:+.4f}"],
                    ["Baseline TS vs CT", f"{w_base['n']:,}", "—", f"{w_base['ks_D']:.3f}", f"{w_base['emd']:.4f}", "—", "—"],
                    [f"{other_label} TS vs CT", f"{w_oth['n']:,}", "—", f"{w_oth['ks_D']:.3f}", f"{w_oth['emd']:.4f}", "—", "—"],
                    ["Δ(EMD) TS↔CT (other - base)", "", "—", "—", f"{(w_oth['emd']-w_base['emd']):+.4f}", "—", "—"],
                ]

                tblc = Table(conv_rows, hAlign="LEFT")
                tblc.setStyle(TableStyle([
                    ("FONTNAME", (0,0), (-1,0), "Helvetica-Bold"),
                    ("BACKGROUND", (0,0), (-1,0), colors.HexColor("#F3F3F3")),
                    ("GRID", (0,0), (-1,-1), 0.3, colors.grey),
                    ("ALIGN", (1,1), (-1,-1), "RIGHT"),
                ]))
                story.append(Spacer(1, 0.1 * inch))
                story.append(Paragraph("TS/CT similarity (preservation & convergence)", h2))
                story.append(tblc)
                story.append(Spacer(1, 0.1 * inch))

                # Visuals: overlay histograms
                def _overlay_hist(valsA, labA, valsB, labB, title):
                    df = pd.DataFrame({"value": np.concatenate([valsA, valsB]),
                                       "which": ([labA]*valsA.size)+([labB]*valsB.size)})
                    fig = px.histogram(df, x="value", color="which",
                                       barmode="overlay", nbins=60, opacity=0.6,
                                       title=title)
                    fig.update_layout(xaxis_title="Adjacency / correlation",
                                      yaxis_title="Count",
                                      width=900, height=420,
                                      margin=dict(l=50, r=20, t=60, b=50))
                    return fig

                figs = [
                    _overlay_hist(base_TS, "baseline TS", base_CT, "baseline CT",
                                  "Within baseline: TS vs CT"),
                    _overlay_hist(oth_TS,  f"{other_label} TS",  oth_CT,  f"{other_label} CT",
                                  f"Within {other_label}: TS vs CT"),
                    _overlay_hist(x_TS, "baseline TS", y_TS, f"{other_label} TS",
                                  "Cross-file: TS baseline vs TS other"),
                    _overlay_hist(x_CT, "baseline CT", y_CT, f"{other_label} CT",
                                  "Cross-file: CT baseline vs CT other"),
                ]
                # Save static for PDF, embed interactive for HTML
                for i, fig in enumerate(figs, 1):
                    png = os.path.join(tmpdir, f"tsct_{i}_{baseline_label}_{other_label}.png")
                    fig.write_image(png, scale=2)
                    story.append(_image_flowable(png, max_img_width))
                    story.append(Spacer(1, 0.06 * inch))
                    if out_html_path:
                        html_sections.append(_wrap_fig_html(fig, title=fig.layout.title.text))

                if out_html_path:
                    html_sections.append(_table_html(conv_rows, caption=f"TS/CT similarity — {other_label} vs baseline"))
                    html_sections.append("<div class='divider'></div>")
                # === KS-D summary diagrams for this compared file ===
                if ks_records:
                    dfk = pd.DataFrame(ks_records)

                    # 2a) Horizontal bar: KS D per pair (sorted)
                    dfk_sort = dfk.sort_values("ks_D", ascending=False)
                    bar_fig = px.bar(
                        dfk_sort, y="pair", x="ks_D", orientation="h",
                        hover_data={"pair": True, "n_aligned": ":,"},
                        title=f"KS D per pair — {other_label} vs baseline"
                    )
                    bar_fig.update_layout(
                        xaxis_title="KS D (0–1; higher = more different)",
                        yaxis_title="Pair (A×B)",
                        width=900,
                        height=max(420, min(1400, 24 * len(dfk_sort))),  # auto-height
                        margin=dict(l=180, r=20, t=60, b=40),
                    )

                    # Save static for PDF; embed interactive for HTML
                    png1 = os.path.join(tmpdir, f"ksD_bar_{baseline_label}_{other_label}.png")
                    bar_fig.write_image(png1, scale=2)
                    story.append(_image_flowable(png1, max_img_width)); story.append(Spacer(1, 0.06 * inch))
                    if out_html_path:
                        html_sections.append(_wrap_fig_html(bar_fig, title=bar_fig.layout.title.text))

                story.append(PageBreak())
                if out_html_path:
                    html_sections.append(_table_html(comp_rows, caption=f"Comparison vs baseline — {other_label}"))
                    html_sections.append("<div class='divider'></div>")

            # Global median deltas across all compared files vs baseline
            # Aggregate over comp_rows captured above for all files
            # Re-scan results to recompute global deltas succinctly
            global_rows = [["Compared file", "Median Δmean", "Median Δmedian", "Median Δprop≥0.5", "Median Pearson r"]]

            for other_path in csv_paths[1:]:
                other_label = _display_name_for_csv(other_path, file_display_names, baseline=False)
                base_pairs = by_label.get(baseline_label, {})
                other_pairs = by_label.get(other_label, {})
                deltas = []
                for (ra_n, rb_n), base_pr in base_pairs.items():
                    other_pr = other_pairs.get((ra_n, rb_n)) or other_pairs.get((rb_n, ra_n))
                    if not other_pr:
                        continue
                    # Quick scalar deltas from pair stats (approx; not aligned-cell exact)
                    # Using summary stats: mean, median, prop≥0.5
                    dm = other_pr.stats["val_mean"] - base_pr.stats["val_mean"]
                    dmed = other_pr.stats["val_median"] - base_pr.stats["val_median"]
                    dprop = other_pr.stats["prop_ge_0.5"] - base_pr.stats["prop_ge_0.5"]
                    deltas.append([dm, dmed, dprop])
                if deltas:
                    arr = np.array(deltas)
                    global_rows.append([
                        other_label,
                        f"{np.nanmedian(arr[:,0]):+.4f}",
                        f"{np.nanmedian(arr[:,1]):+.4f}",
                        f"{np.nanmedian(arr[:,2]):+.3f}",
                        "—",
                    ])

            tblg = Table(global_rows, hAlign="LEFT")
            tblg.setStyle(TableStyle([
                ("FONTNAME", (0,0), (-1,0), "Helvetica-Bold"),
                ("BACKGROUND", (0,0), (-1,0), colors.lightgrey),
                ("GRID", (0,0), (-1,-1), 0.3, colors.grey),
                ("ALIGN", (1,1), (-1,-1), "RIGHT"),
            ]))
            story.append(Paragraph("Global summary (median deltas over shared pairs):", body))
            # === Global KS-D heatmap (other vs baseline across all files) ===
            if ks_global:
                dfG = pd.DataFrame(ks_global)
                if not dfG.empty:
                    pivot = dfG.pivot(index="pair", columns="file", values="ks_D")
                    # Stable ordering: highest median KS D first
                    row_order = pivot.median(axis=1, skipna=True).sort_values(ascending=False).index
                    pivot = pivot.loc[row_order]

                    heat_fig = px.imshow(
                        pivot.to_numpy(),
                        x=list(pivot.columns),
                        y=list(pivot.index),
                        color_continuous_scale="Viridis",
                        aspect="auto",
                        title="Global KS D heatmap — compared files vs baseline"
                    )
                    heat_fig.update_layout(
                        width=1000,
                        height=max(420, min(1600, 22 * len(pivot.index))),
                        coloraxis_colorbar=dict(title="KS D"),
                        margin=dict(l=220, r=20, t=60, b=60),
                    )
                    # Save/Embed
                    png = os.path.join(tmpdir, "ksD_global_heatmap.png")
                    heat_fig.write_image(png, scale=2)
                    story.append(_image_flowable(png, max_img_width))
                    story.append(Spacer(1, 0.12 * inch))
                    if out_html_path:
                        html_sections.append(_wrap_fig_html(heat_fig, title=heat_fig.layout.title.text))
                        html_sections.append("<div class='divider'></div>")

            story.append(Spacer(1, 0.06 * inch))
            story.append(tblg)
            story.append(PageBreak())
            if out_html_path:
                html_sections.append("<section class='section'>")
                html_sections.append("<h2>Global summary (median deltas over shared pairs)</h2>")
                html_sections.append(_table_html(global_rows))
                html_sections.append("</section><div class='divider'></div>")
                # ---- (5D) Across files: pooled EMD + within-file TS↔CT EMD ------------------
                html_sections.append("<section class='section'>")
                html_sections.append("<h2>Across files: pooled distances & within-file TS↔CT</h2>")

                # pooled K per file (TS+CT) over canonical pairs
                pooled_vals: Dict[str, np.ndarray] = {}
                for p in csv_paths:
                    label = _display_name_for_csv(p, file_display_names, baseline=(p == baseline_path))
                    TSv, CTv = _within_TS_CT_vectors(p, canonical)
                    pooled_vals[label] = _concat_finite([TSv, CTv])

                labels = list(pooled_vals.keys())
                M = np.empty((len(labels), len(labels)), float)
                for i,a in enumerate(labels):
                    for j,b in enumerate(labels):
                        xa, xb = pooled_vals[a], pooled_vals[b]
                        M[i,j] = sstats.wasserstein_distance(xa, xb) if xa.size and xb.size else np.nan

                rows = [[""] + labels]
                for i,a in enumerate(labels):
                    rows.append([a] + [f"{M[i,j]:.4f}" if np.isfinite(M[i,j]) else "NA" for j in range(len(labels))])
                html_sections.append(_table_html(rows, caption="Pairwise pooled EMD (lower = more similar)"))

                # Within-file TS↔CT EMD + rough correlations vs baseline TS/CT
                rows2 = [["File", "EMD(TS,CT)", "ρ(TS vs baseline TS)", "ρ(CT vs baseline CT)"]]
                base_TS, base_CT = _within_TS_CT_vectors(baseline_path, canonical)
                def _rho(a,b):
                    a = a[np.isfinite(a)]; b = b[np.isfinite(b)]
                    if a.size < 2 or b.size < 2:
                        return np.nan
                    n = min(a.size, b.size, 100000)
                    ai = np.linspace(0, a.size-1, n).astype(int)
                    bi = np.linspace(0, b.size-1, n).astype(int)
                    return float(np.corrcoef(np.sort(a)[ai], np.sort(b)[bi])[0,1])

                for p in csv_paths:
                    label = _display_name_for_csv(p, file_display_names, baseline=(p == baseline_path))
                    TSv, CTv = _within_TS_CT_vectors(p, canonical)
                    emd_tsct = sstats.wasserstein_distance(TSv, CTv) if TSv.size and CTv.size else np.nan
                    rows2.append([label, f"{emd_tsct:.4f}", f"{_rho(base_TS, TSv):.3f}", f"{_rho(base_CT, CTv):.3f}"])

                html_sections.append(_table_html(rows2, caption="Within-file TS↔CT & cross-file correlations vs baseline"))
                html_sections.append("</section><div class='divider'></div>")
                # ---- end (5D) ----------------------------------------------------------------

        # Write interactive HTML if requested
        if out_html_path:
            _write_html_report(out_html_path, html_sections, title_text=report_title)

        # finally build
        doc = SimpleDocTemplate(
            out_pdf_path, pagesize=A4,
            leftMargin=left_margin, rightMargin=right_margin,
            topMargin=top_margin, bottomMargin=bottom_margin
        )
        doc.build(story)

    return results
