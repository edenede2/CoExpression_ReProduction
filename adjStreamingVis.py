# -*- coding: utf-8 -*-
import random
from typing import Dict, List, Optional, Tuple

import numpy as np
import polars as pl
import pandas as pd
import plotly.express as px


def plot_region_pair_submatrix(
    csv_path: str,
    region_a: str,
    region_b: str,
    max_rows: Optional[int] = 300,
    max_cols: Optional[int] = 300,
    seed: int = 42,
    value_range: Optional[Tuple[float, float]] = None,
    color_scale: str = "Viridis",
    quiet: bool = False,
    ct: int = 3,
    ts: int = 2,
    chunk_size: int = 2000,   # NEW: for full-block streaming sum
):
    """
    Stream a large correlation/adjacency CSV and visualize a sampled submatrix between two regions.
    Also computes a full (unsampled) A×B block sum equivalent to your streaming script.
    """
    random.seed(seed)

    # --------- helpers ---------
    def _normalize(s: str) -> str:
        return s.replace("_", "").upper()

    def _extract_region_prefix(col: str) -> Optional[str]:
        # Region prefix is the substring before the first occurrence of "_ENSG"
        if "_ENSG" not in col:
            return None
        return col.split("_ENSG", 1)[0]

    # Full-block streaming sum (unsampled), same policy as other script
    def _full_block_sum(rows_prefix: str, cols_all: List[str]) -> float:
        total = 0.0
        for k in range(0, len(cols_all), chunk_size):
            chunk = cols_all[k : k + chunk_size]
            if not chunk:
                continue
            lf = (
                pl.scan_csv(csv_path)
                  .select([row_id_col] + chunk)
                  .filter(pl.col(row_id_col).str.starts_with(f"{rows_prefix}_"))
                  # sum each selected column; treat nulls as 0, cast to float
                  .select([pl.col(c).cast(pl.Float64).fill_null(0.0).sum().alias(c) for c in chunk])
                  # now sum horizontally to get total over this chunk
                  .select(pl.sum_horizontal(pl.all()).alias("tot"))
            )
            # collect streaming and accumulate
            total += float(lf.collect(streaming=True)["tot"][0] or 0.0)
        return total

    # --------- scan header only ---------
    cols = pl.scan_csv(csv_path).columns
    if not cols:
        raise ValueError("CSV appears to have no columns.")
    row_id_col = cols[0]
    data_cols = [c for c in cols if c != row_id_col]

    # Map region prefix -> columns
    prefix_to_cols: Dict[str, List[str]] = {}
    for c in data_cols:
        pre = _extract_region_prefix(c)
        if pre:
            prefix_to_cols.setdefault(pre, []).append(c)

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
              f"Region B (cols): {region_b!r} -> {raw_b}")
        print(f"[Header] Columns for {raw_b}: {len(cols_b_full)}")

    # --------- collect row IDs for region A (all, unsampled) ---------
    row_prefix = f"{raw_a}_"
    a_rows_all = (
        pl.scan_csv(csv_path)
          .select(row_id_col)
          .filter(pl.col(row_id_col).str.starts_with(row_prefix))
          .collect(streaming=True)[row_id_col]
          .to_list()
    )
    if not a_rows_all:
        raise ValueError(f"No rows found with prefix '{row_prefix}' in {row_id_col}.")

    # --------- sampling for plots ---------
    if max_rows and len(a_rows_all) > max_rows:
        a_rows_sampled = sorted(random.sample(a_rows_all, max_rows))
    else:
        a_rows_sampled = a_rows_all

    if max_cols and len(cols_b_full) > max_cols:
        b_cols_sampled = sorted(random.sample(cols_b_full, max_cols))
    else:
        b_cols_sampled = cols_b_full

    if not quiet:
        print(f"[Sample] Rows A: {len(a_rows_sampled)} / {len(a_rows_all)}  |  "
              f"Cols B: {len(b_cols_sampled)} / {len(cols_b_full)}")

    # --------- fetch sampled submatrix lazily (for plots) ---------
    lf_sub = (
        pl.scan_csv(csv_path)
          .select([row_id_col] + b_cols_sampled)
          .filter(pl.col(row_id_col).is_in(a_rows_sampled))
    )
    sub_df = lf_sub.collect(streaming=True).to_pandas()
    sub_df = sub_df.set_index(row_id_col).loc[a_rows_sampled]

    # --------- sampled stats (current behavior) ---------
    sub_np = sub_df.to_numpy(dtype=float)
    row_k = np.nansum(sub_np, axis=1)
    col_k = np.nansum(sub_np, axis=0)

    stats = {
        "rows": sub_df.shape[0],
        "cols": sub_df.shape[1],
        "row_k_mean": float(np.nanmean(row_k)) if row_k.size else float("nan"),
        "row_k_median": float(np.nanmedian(row_k)) if row_k.size else float("nan"),
        "row_k_sum": float(np.nansum(row_k)) if row_k.size else float("nan"),
        "col_k_mean": float(np.nanmean(col_k)) if col_k.size else float("nan"),
        "col_k_median": float(np.nanmedian(col_k)) if col_k.size else float("nan"),
        "col_k_sum": float(np.nansum(col_k)) if col_k.size else float("nan"),
    }

    # --------- FULL (unsampled) A×B block totals — matches your other script ---------
    full_sum = _full_block_sum(raw_a, cols_b_full)
    n_rows_full = len(a_rows_all)
    n_cols_full = len(cols_b_full)
    full_block = {
        "rows": n_rows_full,
        "cols": n_cols_full,
        "sum": full_sum,                         # == block_total_sum(raw_a, raw_b)
        "mean_k_rows": full_sum / max(1, n_rows_full),  # == mean_k_CT[(A,B)] in your other script
        "mean_k_cols": full_sum / max(1, n_cols_full),  # optional: mean degree per column
    }

    # --------- distribution (sampled) ---------
    flat_vals = sub_np.ravel()
    flat_vals = flat_vals[~np.isnan(flat_vals)]

    # --------- plots ---------
    title_main = f"Adjacency: {raw_a} (rows) × {raw_b} (cols) | CTbeta {ct}, TSbeta {ts}"
    subtitle = (
        f"Rows k (A→B): sum={stats['row_k_sum']:.3f}, mean={stats['row_k_mean']:.3f}, median={stats['row_k_median']:.3f} | "
        f"Cols k (B→A): sum={stats['col_k_sum']:.3f}, mean={stats['col_k_mean']:.3f}, median={stats['col_k_median']:.3f} | "
        f"sum row wise (full A×B block)={full_block['sum']:.3f}, mean_k_rows (full)={full_block['mean_k_rows']:.3f}"
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
        xaxis_title=f"{raw_b} genes (sampled {sub_df.shape[1]})",
        yaxis_title=f"{raw_a} genes (sampled {sub_df.shape[0]})",
        width=1000,
        height=800,
        coloraxis_colorbar=dict(title="Adj"),
        margin=dict(l=80, r=20, t=90, b=80),
    )

    hist_fig = px.histogram(
        x=flat_vals,
        nbins=60,
        title=f"Correlation/Adjacency Distribution: {raw_a}×{raw_b} (sampled {flat_vals.size:,} cells)",
    )
    hist_fig.update_layout(
        xaxis_title="Correlation / adjacency",
        yaxis_title="Count",
        bargap=0.02,
        width=900,
        height=450,
        margin=dict(l=60, r=20, t=70, b=60),
    )

    if not quiet:
        print(f"[Done] submatrix shape: {sub_df.shape} | histogram cells: {flat_vals.size:,}")
        print(f"[Full] rows={n_rows_full}, cols={n_cols_full}, sum={full_sum:.6g}, "
              f"mean_k_rows={full_block['mean_k_rows']:.6g}")

    return {
        "heatmap_fig": heatmap_fig,
        "hist_fig": hist_fig,
        "stats": stats,                 # sampled stats (unchanged)
        "full_block": full_block,       # NEW: matches streaming script values
        "sampled_row_ids": a_rows_sampled,
        "sampled_col_names": b_cols_sampled,
        "resolved_prefix_a": raw_a,
        "resolved_prefix_b": raw_b,
    }
