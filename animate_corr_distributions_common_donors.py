
"""
Animate within-tissue (TS) and cross-tissue (CT) correlation DISTRIBUTIONS
vs WGCNA power (beta), **forcing a common set of donors across tissues**.

Key changes vs your original:
- Adds robust CSV loaders that (a) detect orientation, (b) coerce numeric, and
  (c) optionally aggregate technical replicates into donor-level rows.
- Aligns **all** tissues to the same intersection of donor IDs before computing
  TS/CT correlations — guaranteeing comparability across tissues/pairs.
- Provides a `build_ts_ct_correlations_common_donors(...)` that mirrors your
  builder, and a convenience `animate_from_csvs_common_donors(...)` that plugs
  the alignment + animation together for quick runs.

Assumptions (configurable):
- Each CSV is an expression matrix (donors × genes) or (genes × donors).
- Donor ID is in the index or in the first column (if present).
- Gene columns are numeric; any non-numeric are coerced to NaN and then dropped.
- Correlations are computed on donor-aligned matrices; absolute Pearson by default.

Outputs
- Two interactive HTML files saved to disk (same as your original):
  1) <out_html_prefix>__TS.html   (one facet per tissue)
  2) <out_html_prefix>__CT.html   (one facet per tissue‑pair)

TIP
- For very large matrices, pass `sample_per_group=...` to downsample values before plotting.
- If your correlations are signed, set ABS_INPUT=False below so we take abs() only once here.

Author: ChatGPT (common-donors edition)
"""
from __future__ import annotations

from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple, Callable
from pathlib import Path
import re
import numpy as np
import pandas as pd

import plotly.graph_objects as go
from plotly.subplots import make_subplots

# ---------------------------- Config ----------------------------
ABS_INPUT = True  # If True, we expect to work with |corr| in [0,1] for visualization.
                  # Set to False to take absolute values inside the flatteners.

# --------------------- Helpers: RNG & sampling ------------------

def _rng(seed: Optional[int]) -> np.random.Generator:
    return np.random.default_rng(seed) if seed is not None else np.random.default_rng()


def sample_vec(v: np.ndarray, max_n: Optional[int], seed: Optional[int]) -> np.ndarray:
    """Downsample vector `v` without replacement to length <= max_n."""
    if (max_n is None) or (v.size <= max_n):
        return v
    rng = _rng(seed)
    idx = rng.choice(v.size, size=max_n, replace=False)
    return v[idx]


# --------------------- CSV loading & alignment ------------------

_gene_like = re.compile(r'^(ENSG|ENST|Gene|G\d+|AC\d+|LOC\d+)', re.IGNORECASE)

def _maybe_set_index(df: pd.DataFrame) -> pd.DataFrame:
    """If the first column looks like an ID column (non-numeric, unique),
    set it as index. Otherwise return df unchanged.
    """
    if df.shape[1] == 0:
        return df
    first_col = df.columns[0]
    col = df[first_col]
    # A candidate ID column: mostly non-numeric and unique
    non_numeric_frac = pd.to_numeric(col, errors='coerce').isna().mean()
    if non_numeric_frac > 0.8 and col.is_unique:
        df = df.set_index(first_col)
    return df


def _coerce_numeric(df: pd.DataFrame) -> pd.DataFrame:
    """Try to coerce everything (except index) to numeric; drop all-NaN cols/rows."""
    # Preserve index, but ensure it's string for donor-id ops
    df = df.copy()
    df.index = df.index.astype(str)
    for c in df.columns:
        df[c] = pd.to_numeric(df[c], errors='coerce')
    # Drop empty / constant columns
    df = df.dropna(axis=1, how='all')
    # Drop rows with all NaNs
    df = df.dropna(axis=0, how='all')
    return df


def _detect_orientation(df: pd.DataFrame) -> str:
    """Return 'donors_by_genes' if donors are rows and genes are columns.
    Otherwise return 'genes_by_donors'.
    Heuristics:
      - If many column names look gene-like, assume donors are rows.
      - Else if many index values look gene-like, assume donors are columns.
      - Else fallback: donors are rows if n_rows < n_cols (common in bulk RNA-seq).
    """
    cols_gene_like = np.mean([bool(_gene_like.match(str(c))) for c in df.columns]) if df.columns.size else 0.0
    idx_gene_like  = np.mean([bool(_gene_like.match(str(i))) for i in df.index]) if df.index.size else 0.0
    if cols_gene_like >= 0.1 and cols_gene_like > idx_gene_like:
        return 'donors_by_genes'
    if idx_gene_like >= 0.1 and idx_gene_like > cols_gene_like:
        return 'genes_by_donors'
    # shape fallback
    return 'donors_by_genes' if df.shape[0] < df.shape[1] else 'genes_by_donors'


def _to_donors_by_genes(df: pd.DataFrame) -> pd.DataFrame:
    """Ensure table is (donors × genes)."""
    orient = _detect_orientation(df)
    if orient == 'genes_by_donors':
        df = df.transpose()
    # Drop any columns with zero variance (cannot be correlated)
    nunq = df.nunique(dropna=True)
    keep = nunq[nunq > 1].index
    if len(keep) < df.shape[1]:
        df = df.loc[:, keep]
    return df


def _aggregate_technical_reps(df: pd.DataFrame, donor_id_func: Optional[Callable[[str], str]] = None) -> pd.DataFrame:
    """If donor IDs include technical replicate suffixes, aggregate to donor-level by mean.
    donor_id_func maps each index string to a 'base' donor ID. If None, identity is used.
    """
    if donor_id_func is None:
        donor_id_func = lambda x: x
    donors = pd.Index([donor_id_func(idx) for idx in df.index], name='donor')
    df = df.copy()
    df.index = donors
    # Average replicates
    df = df.groupby(level=0).mean(numeric_only=True)
    return df


def read_expression_csv(
    path: str | Path,
    donor_id_func: Optional[Callable[[str], str]] = None,
) -> pd.DataFrame:
    """Robust CSV reader that returns donors × genes (float), donor-aggregated if needed."""
    df = pd.read_csv(path)
    df = _maybe_set_index(df)
    df = _coerce_numeric(df)
    df = _to_donors_by_genes(df)
    df = _aggregate_technical_reps(df, donor_id_func=donor_id_func)
    # Final cleanup: drop columns with any NaNs (post-aggregation) to standardize Ns
    df = df.dropna(axis=1, how='any')
    return df


def align_to_common_donors(expr_by_tissue: Mapping[str, pd.DataFrame]) -> Tuple[Dict[str, pd.DataFrame], List[str]]:
    """Intersect donor sets across all tissues and subset each matrix accordingly."""
    donor_sets = [set(df.index) for df in expr_by_tissue.values() if df is not None and len(df)]
    if not donor_sets:
        return {k: v for k, v in expr_by_tissue.items()}, []
    common = sorted(set.intersection(*donor_sets))
    aligned = {t: df.loc[common].sort_index() for t, df in expr_by_tissue.items()}
    return aligned, common


# --------------------- Correlation builders ---------------------

def _corrcoef_abs(A: np.ndarray, B: Optional[np.ndarray] = None) -> np.ndarray:
    """Absolute Pearson correlation. If B is None, return corr(A) square.
    A and B are (n_samples × n_features)."""
    if B is None:
        C = np.corrcoef(A, rowvar=False)
        return np.abs(C)
    # Cross-correlation block using corrcoef
    AB = np.corrcoef(A, B, rowvar=False)
    p = A.shape[1]
    q = B.shape[1]
    block = AB[:p, p:]
    return np.abs(block)


def build_ts_ct_correlations_common_donors(
    tissue_names: Sequence[str],
    tissue_files: Sequence[str | Path],
    *,
    donor_id_func: Optional[Callable[[str], str]] = None,
    cor_method: str = "pearson",
    verbose: bool = True,
) -> Tuple[Dict[str, pd.DataFrame], Dict[str, pd.DataFrame]]:
    """Compute TS and CT correlation matrices using the **same donor set** across tissues.
    Returns (TS_corrs, CT_corrs) where keys:
      - TS_corrs: tissue -> (genes×genes) DataFrame
      - CT_corrs: "Ti||Tj" -> (Gi×Gj) DataFrame
    Notes
    -----
    - Currently supports Pearson correlations. Spearman support could be added if needed.
    - All correlations are absolute (|r|) for consistency with WGCNA adjacency (use abs).
    """
    if cor_method.lower() not in ("pearson",):
        raise NotImplementedError("Only Pearson is supported here.")

    # Load all expression matrices (donors × genes)
    expr_by_tissue: Dict[str, pd.DataFrame] = {}
    for t, f in zip(tissue_names, tissue_files):
        df = read_expression_csv(f, donor_id_func=donor_id_func)
        if verbose:
            print(f"[load] {t}: donors={df.shape[0]}, genes={df.shape[1]} (from {Path(f).name})")
        expr_by_tissue[t] = df

    # Align to common donors
    aligned, common = align_to_common_donors(expr_by_tissue)
    if verbose:
        print(f"[align] Common donors across {len(expr_by_tissue)} tissues: {len(common)}")
        if len(common) == 0:
            print("WARNING: No common donors — results will be empty.")

    # Compute TS correlations per tissue
    TS_corrs: Dict[str, pd.DataFrame] = {}
    for t, df in aligned.items():
        X = df.to_numpy(copy=False)  # (n_donors × n_genes)
        C = _corrcoef_abs(X)  # square
        TS_corrs[t] = pd.DataFrame(C, index=df.columns, columns=df.columns)

    # Compute CT correlations for each pair
    tissues = list(aligned.keys())
    CT_corrs: Dict[str, pd.DataFrame] = {}
    for i in range(len(tissues)):
        for j in range(i + 1, len(tissues)):
            ti, tj = tissues[i], tissues[j]
            Ai = aligned[ti].to_numpy(copy=False)
            Bj = aligned[tj].to_numpy(copy=False)
            block = _corrcoef_abs(Ai, Bj)  # (Gi × Gj)
            CT_corrs[f"{ti}||{tj}"] = pd.DataFrame(block, index=aligned[ti].columns, columns=aligned[tj].columns)

    return TS_corrs, CT_corrs


# ----------------- Flatteners for histogramming -----------------

def flatten_upper_triangle(C: pd.DataFrame) -> np.ndarray:
    """Return upper‑triangle (i<j) vector of a square matrix. Clips to [0,1]."""
    m = C.to_numpy(copy=False)
    n = m.shape[0]
    iu = np.triu_indices(n, k=1)
    v = m[iu]
    if not ABS_INPUT:
        v = np.abs(v)
    v = v[np.isfinite(v)]
    v = v[(v >= 0) & (v <= 1)]
    return v.astype(np.float64, copy=False)


def flatten_rect(C: pd.DataFrame) -> np.ndarray:
    """Return all values from a rectangular CT matrix. Clips to [0,1]."""
    v = C.to_numpy(copy=False).ravel()
    if not ABS_INPUT:
        v = np.abs(v)
    v = v[np.isfinite(v)]
    v = v[(v >= 0) & (v <= 1)]
    return v.astype(np.float64, copy=False)


# ----------------------- Histograms per beta ---------------------

def compute_histograms_for_betas(
    values_by_group: Mapping[str, np.ndarray],
    betas: Sequence[float],
    bins: int = 40,
    density: bool = True,
) -> Tuple[Dict[float, Dict[str, np.ndarray]], np.ndarray]:
    """For each beta, compute histogram counts for each group.
    Returns (frame_counts, bin_edges)
    """
    bin_edges = np.linspace(0.0, 1.0, bins + 1)
    frame_counts: Dict[float, Dict[str, np.ndarray]] = {}
    bin_widths = np.diff(bin_edges)

    for b in betas:
        frame_counts[b] = {}
        for g, v in values_by_group.items():
            if v.size == 0:
                frame_counts[b][g] = np.zeros(bins, dtype=float)
                continue
            w = np.power(v, b)
            counts, _ = np.histogram(w, bins=bin_edges, range=(0.0, 1.0))
            if density:
                n = w.size
                counts = counts.astype(float) / (n * bin_widths)
            frame_counts[b][g] = counts.astype(float)
    return frame_counts, bin_edges


# ------------------- Animated subplot figure --------------------

def make_animated_hist_subplots(
    frame_counts: Mapping[float, Mapping[str, np.ndarray]],
    bin_edges: np.ndarray,
    subplot_titles: Sequence[str],
    title: str,
    height: int = 420,
) -> go.Figure:
    """Create a subplot figure (1 row × N cols), one histogram per group,
    animated over frames keyed by beta.
    """
    groups = list(subplot_titles)
    n = len(groups)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.0

    palette = ["#1f77b4", "#d62728", "#2ca02c", "#9467bd", "#8c564b",
               "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"]
    edge_color = "#2b2b2b"
    bar_opacity = 0.95

    fig = make_subplots(rows=1, cols=n, shared_yaxes=True, horizontal_spacing=0.06,
                        subplot_titles=tuple(groups))

    # Initial traces at the first beta
    betas_sorted = sorted(frame_counts.keys(), key=float)
    b0 = betas_sorted[0] if betas_sorted else 1.0

    for j, g in enumerate(groups, start=1):
        y0 = frame_counts.get(b0, {}).get(g, np.zeros_like(bin_centers))
        color = palette[(j - 1) % len(palette)]
        fig.add_trace(
            go.Bar(x=bin_centers, y=y0, name=g, showlegend=False,
                   marker=dict(color=color, line=dict(color=edge_color, width=1.2)),
                   opacity=bar_opacity),
            row=1, col=j
        )

    # Frames: one per beta, update each subplot trace in order
    frames = []
    for b in betas_sorted:
        data = []
        for j, g in enumerate(groups, start=1):
            y = frame_counts[b].get(g, np.zeros_like(bin_centers))
            color = palette[(j - 1) % len(palette)]
            data.append(
                go.Bar(
                    x=bin_centers,
                    y=y,
                    marker=dict(color=color, line=dict(color=edge_color, width=1.2)),
                    opacity=bar_opacity,
                    showlegend=False,
                )
            )
        frames.append(go.Frame(name=f"beta={b}", data=data))

    fig.update(frames=frames)

    # Slider & play controls
    steps = []
    for k, b in enumerate(betas_sorted):
        steps.append({
            "args": [[f"beta={b}"], {"frame": {"duration": 0, "redraw": True},
                                      "mode": "immediate", "transition": {"duration": 0}}],
            "label": str(b),
            "method": "animate",
        })

    fig.update_layout(
        title=title,
        height=height,
        template="plotly_white",
        plot_bgcolor="white",
        bargap=0.05,
        xaxis_title="adjacency = |corr|^beta",
        yaxis_title="density",
        updatemenus=[{
            "type": "buttons",
            "showactive": True,
            "x": 1.05, "y": 1.15, "xanchor": "right", "yanchor": "top",
            "buttons": [
                {"label": "▶ Play", "method": "animate",
                 "args": [None, {"frame": {"duration": 300, "redraw": True},
                                   "transition": {"duration": 0},
                                   "fromcurrent": True}]},
                {"label": "⏸ Pause", "method": "animate",
                 "args": [[None], {"frame": {"duration": 0, "redraw": False},
                                    "mode": "immediate"}]}
            ]
        }],
        sliders=[{
            "active": 0, "y": -0.08, "x": 0.5, "len": 0.9,
            "xanchor": "center", "yanchor": "top",
            "pad": {"b": 10, "t": 30},
            "steps": steps
        }]
    )

    # Lock axes to [0,1] on x for all subplots
    for i in range(n):
        fig.update_xaxes(range=[0, 1], showline=True, linewidth=1, linecolor="#2b2b2b",
                         gridcolor="#dddddd", row=1, col=i+1)
        fig.update_yaxes(showline=True, linewidth=1, linecolor="#2b2b2b",
                         gridcolor="#dddddd", row=1, col=i+1)
    return fig


# ------------------- High-level animate function ----------------

def animate_ts_ct_distributions(
    TS_corrs: Mapping[str, pd.DataFrame],
    CT_corrs: Mapping[str, pd.DataFrame],
    betas: Sequence[float] = tuple(range(1, 21)),
    sample_per_group: Optional[int] = 200_000,
    bins: int = 40,
    density: bool = True,
    seed: Optional[int] = 0,
    out_html_prefix: str = "corr_beta_anim",
    height: int = 420,
) -> Tuple[str, str]:
    """Create two animated Plotly histograms (TS and CT) over beta."""
    rng = _rng(seed)

    # --- Prepare TS values ---
    def _flatten_ts(C: pd.DataFrame) -> np.ndarray:
        m = C.to_numpy(copy=False)
        n = m.shape[0]
        iu = np.triu_indices(n, k=1)
        v = m[iu]
        if not ABS_INPUT: v = np.abs(v)
        v = v[np.isfinite(v)]
        v = v[(v >= 0) & (v <= 1)]
        return v.astype(np.float64, copy=False)

    ts_values: Dict[str, np.ndarray] = {}
    for tname, C in TS_corrs.items():
        if C is None or C.size == 0:
            ts_values[tname] = np.array([], dtype=float)
            continue
        v = _flatten_ts(C)
        v = sample_vec(v, sample_per_group, seed=rng.integers(0, 2**31 - 1))
        ts_values[tname] = v

    # --- Prepare CT values ---
    ct_values: Dict[str, np.ndarray] = {}
    for pair, C in CT_corrs.items():
        if C is None or C.size == 0:
            ct_values[pair] = np.array([], dtype=float)
            continue
        v = flatten_rect(C)
        v = sample_vec(v, sample_per_group, seed=rng.integers(0, 2**31 - 1))
        ct_values[pair] = v

    # Histograms over betas (TS)
    ts_frame_counts, bin_edges = compute_histograms_for_betas(
        ts_values, betas=betas, bins=bins, density=density
    )
    ts_groups = list(ts_values.keys()) or ["No TS"]
    if not ts_values:
        ts_frame_counts = {b: {"No TS": np.zeros(bins)} for b in betas}
    fig_TS = make_animated_hist_subplots(
        ts_frame_counts, bin_edges, subplot_titles=ts_groups,
        title="Within‑tissue adjacency distribution vs β (common donors)", height=height
    )

    # Histograms over betas (CT)
    ct_frame_counts, bin_edges_ct = compute_histograms_for_betas(
        ct_values, betas=betas, bins=bins, density=density
    )
    ct_groups = list(ct_values.keys()) or ["No CT"]
    if not ct_values:
        ct_frame_counts = {b: {"No CT": np.zeros(bins)} for b in betas}
    fig_CT = make_animated_hist_subplots(
        ct_frame_counts, bin_edges_ct, subplot_titles=ct_groups,
        title="Cross‑tissue adjacency distribution vs β (common donors)", height=height
    )

    # Save HTML files
    out_TS = f"{out_html_prefix}__TS.html"
    out_CT = f"{out_html_prefix}__CT.html"
    fig_TS.write_html(out_TS, include_plotlyjs="cdn", auto_play=False)
    fig_CT.write_html(out_CT, include_plotlyjs="cdn", auto_play=False)

    return out_TS, out_CT


# -------- Convenience: end-to-end from CSVs (common donors) -----

def animate_from_csvs_common_donors(
    tissue_names: Sequence[str],
    tissue_files: Sequence[str | Path],
    *,
    donor_id_func: Optional[Callable[[str], str]] = None,
    betas: Sequence[float] = tuple(range(1, 21)),
    sample_per_group: Optional[int] = 200_000,
    bins: int = 40,
    density: bool = True,
    seed: Optional[int] = 42,
    out_html_prefix: str = "corr_beta_anim_COMMON_DONORS",
    height: int = 420,
    verbose: bool = True,
) -> Tuple[str, str]:
    """One-call pipeline: load CSVs → align common donors → TS/CT corr → animate."""
    TS_corrs, CT_corrs = build_ts_ct_correlations_common_donors(
        tissue_names=tissue_names,
        tissue_files=tissue_files,
        donor_id_func=donor_id_func,
        verbose=verbose,
    )
    return animate_ts_ct_distributions(
        TS_corrs, CT_corrs,
        betas=betas,
        sample_per_group=sample_per_group,
        bins=bins,
        density=density,
        seed=seed,
        out_html_prefix=out_html_prefix,
        height=height,
    )


# ------------------- If run as a script -------------------
if __name__ == "__main__":
    # Example using the filenames you uploaded (OLD cohort).
    # Adjust paths as needed; they must exist on your machine.
    files_old = [
        "/mnt/data/Adipose - Subcutaneous_old.csv",
        "/mnt/data/Muscle - Skeletal_old.csv",
        "/mnt/data/Brain - Cortex_old.csv",
    ]
    tissues = ["Adipose - Subcutaneous", "Muscle - Skeletal", "Brain - Cortex"]

    def donor_base(x: str) -> str:
        """
        Optional donor normalizer. Edit to match your donor naming convention.
        Example mirrors your earlier R helper: keeps the first two hyphen groups.
        """
        m = re.match(r"^([^-]+-[^-]+)", str(x))
        return m.group(1) if m else str(x)

    try:
        out_TS, out_CT = animate_from_csvs_common_donors(
            tissue_names=tissues,
            tissue_files=files_old,
            donor_id_func=donor_base,
            betas=list(range(1, 21)),
            sample_per_group=200_000,
            bins=40,
            density=True,
            seed=42,
            out_html_prefix="/mnt/data/corr_beta_anim_COMMON_DONORS_OLD",
            height=420,
            verbose=True,
        )
        print("Saved:", out_TS, out_CT)
    except Exception as e:
        print("Self-test run failed (this is OK if paths are missing):", e)
