# -*- coding: utf-8 -*-
"""
FAST Multi-report Power-law GoF analyzer for adjacency/correlation matrices.

Speed-oriented changes vs previous version:
- Fast bootstrap: fixed xmin, re-estimate alpha only (no xmin grid per bootstrap)
- Smaller xmin grid (by default)
- Downsample per-block values to max_cells_per_block for fitting/KS
- Read each pair's submatrix ONCE per file; reuse for integrated TS/CT/ALL
- Parallel per-pair analysis (n_jobs > 1)
- Only integrated figures (TS, CT, ALL) saved; per-block visuals summarized as tables

Author: Eden + ChatGPT
"""

from __future__ import annotations

import os
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
from scipy.stats import mannwhitneyu, kruskal
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.io as pio

# PDF
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, Image as RLImage, Table, TableStyle, PageBreak
)
from reportlab.lib.utils import ImageReader

# parallel
from concurrent.futures import ProcessPoolExecutor, as_completed

# ensure plotly static export works
import plotly.io as pio
pio.defaults.default_format = "png"

# Import Barabási Scale-Free Analysis Module
import warnings
import sys
from pathlib import Path

try:
    from barabasi_scale_free_analysis import (
        analyze_scale_free_properties, 
        save_scale_free_analysis_report,
        ScaleFreeAnalysis,
        BarabasiScaleFreeAnalyzer
    )
    BARABASI_AVAILABLE = True
    print("✓ Barabási Scale-Free Analysis module loaded successfully")
except ImportError as e:
    BARABASI_AVAILABLE = False
    print(f"⚠ Barabási module not available: {e}")
    print("  Continuing with standard analysis only")

# -------------------- Enhanced Statistical Metrics --------------------
def _calculate_enhanced_metrics(values: np.ndarray) -> Dict[str, float]:
    """Calculate enhanced statistical metrics for a set of values."""
    if len(values) == 0:
        return {
            'mean': np.nan, 'median': np.nan, 'p5': np.nan, 'p95': np.nan,
            'frac_02': np.nan, 'frac_04': np.nan, 'frac_06': np.nan, 'frac_08': np.nan,
            'std': np.nan, 'iqr': np.nan, 'cv': np.nan, 'skewness': np.nan, 'kurtosis': np.nan
        }
    
    finite_vals = values[np.isfinite(values)]
    if len(finite_vals) == 0:
        return {
            'mean': np.nan, 'median': np.nan, 'p5': np.nan, 'p95': np.nan,
            'frac_02': np.nan, 'frac_04': np.nan, 'frac_06': np.nan, 'frac_08': np.nan,
            'std': np.nan, 'iqr': np.nan, 'cv': np.nan, 'skewness': np.nan, 'kurtosis': np.nan
        }
    
    from scipy.stats import skew, kurtosis
    
    mean_val = np.mean(finite_vals)
    std_val = np.std(finite_vals)
    
    return {
        'mean': float(mean_val),
        'median': float(np.median(finite_vals)),
        'p5': float(np.percentile(finite_vals, 5)),
        'p95': float(np.percentile(finite_vals, 95)),
        'frac_02': float(np.mean(finite_vals > 0.2)),
        'frac_04': float(np.mean(finite_vals > 0.4)),
        'frac_06': float(np.mean(finite_vals > 0.6)),
        'frac_08': float(np.mean(finite_vals > 0.8)),
        'std': float(std_val),
        'iqr': float(np.percentile(finite_vals, 75) - np.percentile(finite_vals, 25)),
        'cv': float(std_val / mean_val if mean_val != 0 else np.nan),
        'skewness': float(skew(finite_vals)) if len(finite_vals) > 2 else np.nan,
        'kurtosis': float(kurtosis(finite_vals)) if len(finite_vals) > 3 else np.nan
    }

def _calculate_emd(values1: np.ndarray, values2: np.ndarray) -> float:
    """Calculate Earth Mover's Distance (Wasserstein-1) between two distributions."""
    try:
        from scipy.stats import wasserstein_distance
        finite1 = values1[np.isfinite(values1)]
        finite2 = values2[np.isfinite(values2)]
        if len(finite1) == 0 or len(finite2) == 0:
            return np.nan
        return float(wasserstein_distance(finite1, finite2))
    except ImportError:
        return np.nan

def _calculate_delta_metrics(metrics1: Dict[str, float], metrics2: Dict[str, float], label1: str, label2: str) -> Dict[str, float]:
    """Calculate delta (difference) metrics between two metric dictionaries."""
    deltas = {}
    for key in metrics1.keys():
        if key in metrics2:
            deltas[f'delta_{key}_{label1}_{label2}'] = metrics1[key] - metrics2[key]
    return deltas

def _visualize_ks_statistic(x_data: np.ndarray, xmin: float, alpha: float, title: str = "KS Test Visualization") -> go.Figure:
    """Create a detailed visualization of the KS statistic to explain the test."""
    if len(x_data) == 0:
        fig = go.Figure()
        fig.add_annotation(text="No data available", xref="paper", yref="paper", x=0.5, y=0.5)
        return fig
    
    # Sort data and filter tail
    x_sorted = np.sort(x_data[x_data >= xmin])
    if len(x_sorted) == 0:
        fig = go.Figure()
        fig.add_annotation(text="No data in tail", xref="paper", yref="paper", x=0.5, y=0.5)
        return fig
    
    n = len(x_sorted)
    
    # Empirical CDF
    empirical_cdf = np.arange(1, n + 1) / n
    
    # Theoretical power-law CDF
    theoretical_cdf = 1 - (x_sorted / xmin) ** (1 - alpha)
    
    # KS statistic (maximum difference)
    ks_differences = np.abs(empirical_cdf - theoretical_cdf)
    ks_stat = np.max(ks_differences)
    ks_idx = np.argmax(ks_differences)
    
    fig = make_subplots(
        rows=2, cols=2,
        subplot_titles=("Empirical vs Theoretical CDF", "KS Differences", "Power-law Fit (Log-Log)", "Data Distribution"),
        specs=[[{"secondary_y": False}, {"secondary_y": False}],
               [{"secondary_y": False}, {"secondary_y": False}]]
    )
    
    # Plot 1: CDFs comparison
    fig.add_trace(go.Scatter(x=x_sorted, y=empirical_cdf, name="Empirical CDF", 
                            line=dict(color='blue', width=2)), row=1, col=1)
    fig.add_trace(go.Scatter(x=x_sorted, y=theoretical_cdf, name="Theoretical CDF", 
                            line=dict(color='red', width=2)), row=1, col=1)
    fig.add_trace(go.Scatter(x=[x_sorted[ks_idx]], y=[empirical_cdf[ks_idx]], 
                            name=f"Max Diff (KS={ks_stat:.4f})", mode='markers',
                            marker=dict(color='orange', size=10)), row=1, col=1)
    
    # Plot 2: KS differences
    fig.add_trace(go.Scatter(x=x_sorted, y=ks_differences, name="KS Differences", 
                            line=dict(color='green', width=2), showlegend=False), row=1, col=2)
    fig.add_hline(y=ks_stat, line_dash="dash", line_color="orange", row=1, col=2,
                  annotation_text=f"KS Statistic = {ks_stat:.4f}")
    
    # Plot 3: Log-log power-law fit
    log_x = np.log10(x_sorted)
    log_ccdf = np.log10(1 - theoretical_cdf + 1e-10)  # Add small value to avoid log(0)
    empirical_ccdf = 1 - empirical_cdf
    log_emp_ccdf = np.log10(empirical_ccdf + 1e-10)
    
    fig.add_trace(go.Scatter(x=log_x, y=log_emp_ccdf, name="Empirical CCDF", 
                            mode='markers', marker=dict(color='blue', size=4), showlegend=False), row=2, col=1)
    fig.add_trace(go.Scatter(x=log_x, y=log_ccdf, name="Power-law Fit", 
                            line=dict(color='red', width=2), showlegend=False), row=2, col=1)
    
    # Plot 4: Data histogram
    fig.add_trace(go.Histogram(x=x_sorted, name="Data Distribution", 
                              marker=dict(color='lightblue', opacity=0.7), showlegend=False), row=2, col=2)
    fig.add_vline(x=xmin, line_dash="dash", line_color="red", row=2, col=2,
                  annotation_text=f"xmin = {xmin:.3f}")
    
    # Update layout
    fig.update_xaxes(title_text="Value", row=1, col=1)
    fig.update_yaxes(title_text="CDF", row=1, col=1)
    fig.update_xaxes(title_text="Value", row=1, col=2)
    fig.update_yaxes(title_text="KS Difference", row=1, col=2)
    fig.update_xaxes(title_text="log₁₀(Value)", row=2, col=1)
    fig.update_yaxes(title_text="log₁₀(CCDF)", row=2, col=1)
    fig.update_xaxes(title_text="Value", row=2, col=2)
    fig.update_yaxes(title_text="Count", row=2, col=2)
    
    fig.update_layout(
        title=f"{title}<br><sub>α={alpha:.3f}, xmin={xmin:.3f}, KS={ks_stat:.4f}</sub>",
        width=1000, height=800,
        showlegend=True
    )
    
    return fig

# -------------------- P-value Formatting Functions --------------------

def _format_pvalue(p: float, scientific_threshold: float = 0.001) -> str:
    """Format p-value with appropriate notation - scientific for small values."""
    if np.isnan(p) or p is None:
        return "N/A"
    
    if p == 0:
        return "< 1e-16"
    
    if p < scientific_threshold:
        return f"{p:.2e}"
    else:
        return f"{p:.4f}"

def _format_pvalue_hover(p: float) -> str:
    """Format p-value for hover text with consistent scientific notation."""
    if np.isnan(p) or p is None:
        return "N/A"
    
    if p == 0:
        return "< 1e-16"
    
    if p < 0.001:
        return f"{p:.2e}"
    else:
        return f"{p:.4f}"

# -------------------- HTML helpers --------------------
def _html_escape(s: str) -> str:
    return html.escape(s, quote=True)

def _table_html(rows: List[List[str]], caption: Optional[str] = None) -> str:
    if not rows: return ""
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
    """Enhanced HTML report writer with better styling and navigation."""
    os.makedirs(os.path.dirname(out_html_path) or ".", exist_ok=True)
    
    # Count sections for navigation
    section_count = len(html_sections)
    has_plots = any('plotly-div' in section for section in html_sections)
    
    head = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8"/>
<meta name="viewport" content="width=device-width, initial-scale=1"/>
<title>{_html_escape(title_text)} - Power-Law Analysis Report</title>
<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
<link href="https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&display=swap" rel="stylesheet">
<style>
  :root {{ 
    --primary: #2563eb; --primary-dark: #1d4ed8; --primary-light: #3b82f6;
    --secondary: #64748b; --success: #10b981; --warning: #f59e0b; --danger: #ef4444;
    --bg-primary: #ffffff; --bg-secondary: #f8fafc; --bg-accent: #f1f5f9;
    --text-primary: #0f172a; --text-secondary: #475569; --text-muted: #64748b;
    --border: #e2e8f0; --border-light: #f1f5f9;
    --shadow-sm: 0 1px 2px 0 rgb(0 0 0 / 0.05); --shadow-md: 0 4px 6px -1px rgb(0 0 0 / 0.1);
    --radius: 8px; --radius-lg: 12px;
  }}
  
  * {{ box-sizing: border-box; }}
  
  body {{ 
    font-family: 'Inter', -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
    margin: 0; padding: 0; color: var(--text-primary); background: var(--bg-secondary);
    line-height: 1.6; font-size: 14px;
  }}
  
  .container {{ max-width: 1400px; margin: 0 auto; padding: 0 20px; }}
  
  .header {{ 
    background: linear-gradient(135deg, var(--primary) 0%, var(--primary-dark) 100%);
    color: white; padding: 2rem 0; margin-bottom: 2rem;
    box-shadow: var(--shadow-md);
  }}
  
  .header h1 {{ 
    margin: 0; font-size: 2.25rem; font-weight: 700; line-height: 1.2;
    text-shadow: 0 2px 4px rgb(0 0 0 / 0.1);
  }}
  
  .header .subtitle {{ 
    margin: 0.5rem 0 0 0; font-size: 1.1rem; opacity: 0.9; font-weight: 400;
  }}
  
  .header .meta {{ 
    margin-top: 1rem; display: flex; gap: 2rem; flex-wrap: wrap;
    font-size: 0.9rem; opacity: 0.8;
  }}
  
  .header .meta-item {{ display: flex; align-items: center; gap: 0.5rem; }}
  
  .nav-toc {{ 
    background: var(--bg-primary); border-radius: var(--radius-lg);
    padding: 1.5rem; margin-bottom: 2rem; box-shadow: var(--shadow-sm);
    border: 1px solid var(--border);
  }}
  
  .nav-toc h2 {{ 
    margin: 0 0 1rem 0; font-size: 1.25rem; color: var(--primary);
    border-bottom: 2px solid var(--primary-light); padding-bottom: 0.5rem;
  }}
  
  .nav-list {{ 
    display: grid; grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
    gap: 0.5rem; list-style: none; padding: 0; margin: 0;
  }}
  
  .nav-item {{ 
    padding: 0.5rem 0.75rem; border-radius: var(--radius);
    border: 1px solid var(--border-light); transition: all 0.2s;
  }}
  
  .nav-item:hover {{ 
    background: var(--bg-accent); border-color: var(--primary);
    transform: translateY(-1px); box-shadow: var(--shadow-sm);
  }}
  
  .nav-link {{ 
    text-decoration: none; color: var(--text-secondary); font-weight: 500;
    display: block; line-height: 1.4;
  }}
  
  .nav-link:hover {{ color: var(--primary); }}
  
  .main-content {{ background: var(--bg-primary); border-radius: var(--radius-lg); }}
  
  .section {{ 
    margin-bottom: 3rem; padding: 2rem; border-radius: var(--radius-lg);
    background: var(--bg-primary); box-shadow: var(--shadow-sm);
    border: 1px solid var(--border);
  }}
  
  .section:target {{ 
    border-color: var(--primary); box-shadow: 0 0 0 3px rgb(37 99 235 / 0.1);
  }}
  
  .section h2 {{ 
    margin: 0 0 1.5rem 0; font-size: 1.5rem; font-weight: 600;
    color: var(--primary); display: flex; align-items: center; gap: 0.5rem;
  }}
  
  .section h3 {{ 
    margin: 1.5rem 0 1rem 0; font-size: 1.25rem; color: var(--text-primary);
    border-left: 4px solid var(--primary-light); padding-left: 1rem;
  }}
  
  .section h4 {{ 
    margin: 1rem 0 0.5rem 0; font-size: 1rem; color: var(--text-secondary);
    text-transform: uppercase; letter-spacing: 0.05em; font-weight: 500;
  }}
  
  .plot-container {{ 
    margin: 1.5rem 0; border-radius: var(--radius); overflow: hidden;
    border: 1px solid var(--border); background: white;
  }}
  
  .tbl {{ margin: 1.5rem 0; }}
  .tbl table {{ 
    border-collapse: collapse; width: 100%; font-size: 13px;
    border-radius: var(--radius); overflow: hidden; box-shadow: var(--shadow-sm);
  }}
  .tbl th {{ 
    background: var(--bg-accent); color: var(--text-primary); font-weight: 600;
    border: 1px solid var(--border); padding: 12px 16px; text-align: left;
  }}
  .tbl td {{ 
    border: 1px solid var(--border); padding: 10px 16px;
    transition: background-color 0.2s;
  }}
  .tbl tr:hover td {{ background: var(--bg-accent); }}
  .tbl th:first-child, .tbl td:first-child {{ text-align: left; }}
  .tbl th:not(:first-child), .tbl td:not(:first-child) {{ text-align: right; }}
  .tbl caption {{ 
    text-align: left; margin-bottom: 0.5rem; color: var(--text-muted);
    font-weight: 500; font-size: 0.9rem;
  }}
  
  .stats-grid {{ 
    display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
    gap: 1rem; margin: 1rem 0;
  }}
  
  .stat-card {{ 
    background: var(--bg-accent); padding: 1rem; border-radius: var(--radius);
    text-align: center; border: 1px solid var(--border-light);
  }}
  
  .stat-value {{ 
    font-size: 1.75rem; font-weight: 700; color: var(--primary);
    display: block; line-height: 1;
  }}
  
  .stat-label {{ 
    font-size: 0.875rem; color: var(--text-muted); text-transform: uppercase;
    letter-spacing: 0.05em; margin-top: 0.25rem;
  }}
  
  .divider {{ 
    height: 1px; background: var(--border); margin: 2rem 0;
    border-radius: 1px;
  }}
  
  .footer {{ 
    margin-top: 3rem; padding: 2rem 0; text-align: center;
    border-top: 1px solid var(--border); color: var(--text-muted);
    font-size: 0.875rem;
  }}
  
  .badge {{ 
    display: inline-block; padding: 0.25rem 0.5rem; border-radius: 4px;
    font-size: 0.75rem; font-weight: 600; text-transform: uppercase;
    letter-spacing: 0.05em;
  }}
  
  .badge-primary {{ background: var(--primary); color: white; }}
  .badge-success {{ background: var(--success); color: white; }}
  .badge-warning {{ background: var(--warning); color: white; }}
  .badge-secondary {{ background: var(--secondary); color: white; }}
  
  .alert {{ 
    padding: 1rem; border-radius: var(--radius); margin: 1rem 0;
    border-left: 4px solid;
  }}
  
  .alert-info {{ 
    background: rgb(59 130 246 / 0.1); border-color: var(--primary);
    color: var(--primary-dark);
  }}
  
  .alert-success {{ 
    background: rgb(16 185 129 / 0.1); border-color: var(--success);
    color: #065f46;
  }}
  
  .alert-warning {{ 
    background: rgb(245 158 11 / 0.1); border-color: var(--warning);
    color: #92400e;
  }}
  
  @media (max-width: 768px) {{
    .header h1 {{ font-size: 1.75rem; }}
    .header .meta {{ flex-direction: column; gap: 0.5rem; }}
    .nav-list {{ grid-template-columns: 1fr; }}
    .stats-grid {{ grid-template-columns: repeat(2, 1fr); }}
    .section {{ padding: 1.5rem 1rem; }}
  }}
  
  .plotly-graph-div {{ border-radius: var(--radius) !important; }}
</style>
</head>
<body>
<header class="header">
  <div class="container">
    <h1>{_html_escape(title_text)}</h1>
    <p class="subtitle">Comprehensive Power-Law Goodness-of-Fit Analysis Report</p>
    <div class="meta">
      <div class="meta-item">
        <span>📊</span>
        <span>{section_count} Analysis Section{'s' if section_count != 1 else ''}</span>
      </div>
      <div class="meta-item">
        <span>📈</span>
        <span>{'Interactive Plots' if has_plots else 'Statistical Tables'}</span>
      </div>
      <div class="meta-item">
        <span>🔬</span>
        <span>Network Science Analysis</span>
      </div>
      <div class="meta-item">
        <span>⏱️</span>
        <span>Generated: {html.escape(str(pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')))}</span>
      </div>
    </div>
  </div>
</header>

<div class="container">
  <div class="nav-toc">
    <h2>📋 Report Navigation</h2>
    <p style="margin-bottom: 1rem; color: var(--text-muted);">
      Jump to specific sections of this comprehensive power-law analysis report.
    </p>
    <ul class="nav-list">"""
    
    # Add navigation items based on detected sections
    nav_items = []
    for i, section in enumerate(html_sections):
        if '<h2' in section or '<h3' in section:
            # Extract section title (simplified)
            section_title = f"Analysis Section {i+1}"
            if 'Distribution' in section:
                section_title = "📊 Distribution Analysis"
            elif 'Comparison' in section:
                section_title = "📈 Statistical Comparison"
            elif 'Goodness' in section:
                section_title = "🎯 Goodness-of-Fit Tests"
            elif 'table' in section.lower():
                section_title = "📋 Summary Statistics"
                
            nav_items.append(f'<li class="nav-item"><a href="#section-{i+1}" class="nav-link">{section_title}</a></li>')
    
    if not nav_items:
        nav_items = [f'<li class="nav-item"><a href="#section-{i+1}" class="nav-link">📊 Analysis Section {i+1}</a></li>' 
                    for i in range(section_count)]
    
    head += '\n'.join(nav_items[:12])  # Limit to 12 navigation items
    
    head += """
    </ul>
  </div>
  
  <div class="main-content">"""
    
    tail = """
  </div>
  
  <footer class="footer">
    <div class="container">
      <p><strong>Power-Law Analysis Report</strong> | 
         Generated by Enhanced Statistical Analysis Pipeline | 
         <span class="badge badge-primary">Network Science</span>
         <span class="badge badge-secondary">Statistical Modeling</span>
      </p>
      <p style="margin-top: 0.5rem; font-size: 0.8rem;">
        For questions about this analysis, refer to the methodology documentation 
        or the associated research publications.
      </p>
    </div>
  </footer>
</body>
</html>"""

    with open(out_html_path, "w", encoding="utf-8") as f:
        f.write(head)
        for i, sec in enumerate(html_sections):
            # Wrap each section with enhanced styling
            section_wrapper = f"""
    <div class="section" id="section-{i+1}">
      {sec}
    </div>"""
            f.write(section_wrapper)
        f.write(tail)

# -------------------- CSV structure --------------------
def _normalize(s: str) -> str:
    return s.replace("_", "").upper()

def _extract_region_prefix(col: str) -> Optional[str]:
    if "_ENSG" not in col:
        return None
    return col.split("_ENSG", 1)[0]

def _scan_prefix_to_cols(csv_path: str) -> Tuple[str, Dict[str, List[str]]]:
    lf = pl.scan_csv(csv_path)
    cols = lf.collect_schema().names()
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
    """Ultra-optimized scanning for massive files with size limits."""
    out: Dict[str, List[str]] = {pre: [] for pre in prefixes}
    
    # Check file size and add warnings
    file_size_gb = os.path.getsize(csv_path) / (1024**3)
    print(f"      🔍 Scanning {os.path.basename(csv_path)} ({file_size_gb:.1f}GB) for prefixes: {prefixes}")
    
    if file_size_gb > 10:  # For files > 10GB, use aggressive optimization
        print(f"         ⚡ Large file detected - using ultra-fast sampling approach")
        
        try:
            # For massive files, just sample the first N rows to find prefixes
            lf = pl.scan_csv(csv_path)
            
            # Build filter for all prefixes
            prefix_filters = [pl.col(row_id_col).str.starts_with(f"{pre}_") for pre in prefixes]
            combined_filter = prefix_filters[0]
            for pf in prefix_filters[1:]:
                combined_filter = combined_filter | pf
            
            # Sample first 50000 rows only for massive files
            sample_size = 50000 if file_size_gb > 20 else 100000
            print(f"         📊 Sampling first {sample_size:,} rows for pattern detection")
            
            rows_df = (
                lf.select(row_id_col)
                  .filter(combined_filter)
                  .limit(sample_size)  # Critical: limit rows processed
                  .collect(streaming=True)
            )
            
            # Distribute sampled rows to prefixes
            all_rows = rows_df[row_id_col].to_list()
            for row in all_rows:
                for pre in prefixes:
                    if row.startswith(f"{pre}_"):
                        out[pre].append(row)
                        break
                        
            # For massive files, artificially expand the sample if patterns found
            for pre in prefixes:
                if out[pre]:  # If we found some examples
                    # Estimate total count and create synthetic entries if needed
                    sample_count = len(out[pre])
                    if sample_count > 0:
                        # Keep only first 1000 per prefix to avoid memory issues
                        out[pre] = out[pre][:1000]
                        print(f"         📝 {pre}: found {sample_count} examples (limited to {len(out[pre])})")
                    
            print(f"         ✅ Fast scan completed - {sum(len(rows) for rows in out.values())} total rows sampled")
            
        except Exception as e:
            print(f"         ⚠️  Fast scan failed: {e}, using minimal fallback")
            # Ultra-minimal fallback for massive files
            for pre in prefixes:
                out[pre] = [f"{pre}_dummy_0", f"{pre}_dummy_1"]  # Dummy entries
                
    else:
        # Original method for smaller files
        try:
            lf = pl.scan_csv(csv_path)
            prefix_filters = [pl.col(row_id_col).str.starts_with(f"{pre}_") for pre in prefixes]
            combined_filter = prefix_filters[0]
            for pf in prefix_filters[1:]:
                combined_filter = combined_filter | pf
            
            rows_df = (
                lf.select(row_id_col)
                  .filter(combined_filter)
                  .collect(streaming=True)
            )
            
            all_rows = rows_df[row_id_col].to_list()
            for row in all_rows:
                for pre in prefixes:
                    if row.startswith(f"{pre}_"):
                        out[pre].append(row)
                        break
                        
            print(f"         ✅ Found {sum(len(rows) for rows in out.values())} total rows")
            
        except Exception as e:
            print(f"         ⚠️  Fallback to minimal scan due to: {e}")
            for pre in prefixes:
                out[pre] = [f"{pre}_minimal_0"]
    
    return out

def _build_canonical_samples(
    csv_paths: Sequence[str],
    requested_pairs: Sequence[Tuple[str, str]],
    max_rows: int,
    max_cols: int,
    seed: int = 42,
) -> Dict[Tuple[str, str], Tuple[List[str], List[str]]]:
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
    """Optimized subblock fetching for huge files using chunked streaming."""
    import os
    
    # Get file size for progress tracking
    file_size_gb = os.path.getsize(csv_path) / (1024**3)
    print(f"         📊 Reading subblock from {file_size_gb:.1f}GB file...")
    
    try:
        # Use streaming with row filtering for memory efficiency
        lf = pl.scan_csv(csv_path)
        cols = lf.collect_schema().names()
        if not cols:
            return pd.DataFrame(index=[], columns=[])
        
        row_id_col = cols[0]
        keep = [c for c in col_names if c in cols]
        if not keep:
            return pd.DataFrame(index=[], columns=[])
        
        # Optimize query with early filtering and limited columns
        row_ids_set = set(row_ids)  # Convert to set for faster lookup
        
        print(f"         🔍 Filtering {len(row_ids)} rows from {len(keep)} columns...")
        
        # Use streaming with optimized query
        df = (
            lf.select([row_id_col] + keep)
              .filter(pl.col(row_id_col).is_in(row_ids))
              .collect(streaming=True)  # Use streaming to reduce memory usage
              .to_pandas()
              .set_index(row_id_col)
        )
        
        # Preserve order and ensure all requested data
        idx = [r for r in row_ids if r in df.index]
        keep2 = [c for c in keep if c in df.columns]
        result = df.loc[idx, keep2]
        
        print(f"         ✅ Retrieved subblock: {result.shape[0]} rows × {result.shape[1]} cols")
        return result
        
    except Exception as e:
        print(f"         ⚠️ Error reading subblock: {e}")
        print(f"         🔄 Attempting fallback method...")
        
        # Fallback: try pandas chunked reading if polars fails
        try:
            import pandas as pd
            chunk_size = 10000
            found_rows = []
            
            for chunk in pd.read_csv(csv_path, chunksize=chunk_size):
                if len(chunk.columns) == 0:
                    continue
                    
                row_id_col = chunk.columns[0]
                chunk_match = chunk[chunk[row_id_col].isin(row_ids)]
                if not chunk_match.empty:
                    found_rows.append(chunk_match)
                    
                # Early exit if we found all rows
                if len(found_rows) > 0:
                    combined = pd.concat(found_rows, ignore_index=True)
                    if len(combined) >= len(row_ids):
                        break
            
            if found_rows:
                result = pd.concat(found_rows, ignore_index=True).set_index(row_id_col)
                keep2 = [c for c in col_names if c in result.columns]
                idx = [r for r in row_ids if r in result.index]
                return result.loc[idx, keep2]
            else:
                return pd.DataFrame(index=[], columns=[])
                
        except Exception as e2:
            print(f"         ❌ Fallback failed: {e2}")
            return pd.DataFrame(index=[], columns=[])

# -------------------- ECDF/CCDF --------------------
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

# -------------------- Power-law tools --------------------
def _alpha_mle_continuous(x: np.ndarray, xmin: float) -> float:
    n = x.size
    return 1.0 + n / np.sum(np.log(x / xmin))

def _ksD_tail(x: np.ndarray, xmin: float, alpha: float) -> float:
    xs = np.sort(x)
    F_emp = np.arange(1, xs.size + 1, dtype=float) / xs.size
    F_model = 1.0 - (xs / float(xmin)) ** (1.0 - alpha)
    return float(np.max(np.abs(F_emp - F_model)))

def _ks_pvalue_from_statistic(n_tail: int, ks_D: float) -> float:
    """Compute KS p-value from KS statistic and sample size using asymptotic distribution."""
    from scipy import stats
    if n_tail < 30:
        # For small samples, use scipy's exact distribution
        return float(stats.ksone.sf(ks_D, n_tail))
    else:
        # For large samples, use asymptotic approximation
        pval = 2.0 * np.exp(-2.0 * n_tail * ks_D**2)
        return min(1.0, pval)

def _ks_pvalue_tail(x: np.ndarray, xmin: float, alpha: float) -> float:
    """Compute KS p-value for power-law fit using asymptotic distribution."""
    from scipy import stats
    n = len(x)
    D = _ksD_tail(x, xmin, alpha)
    # Use Kolmogorov distribution for p-value approximation
    # For large n, KS statistic follows asymptotic distribution
    if n > 30:
        # Asymptotic p-value (two-sided test)
        pval = 2.0 * np.exp(-2.0 * n * D**2)
        return min(1.0, pval)
    else:
        # For small samples, use scipy's ksone distribution
        return float(stats.ksone.sf(D, n))

def _ks_pvalue_from_statistic(n_tail: int, ks_d_stat: float) -> float:
    """Compute KS p-value from precomputed statistic and sample size."""
    from scipy import stats
    if n_tail > 30:
        # Asymptotic p-value (two-sided test)
        pval = 2.0 * np.exp(-2.0 * n_tail * ks_d_stat**2)
        return min(1.0, pval)
    else:
        # For small samples, use scipy's ksone distribution
        return float(stats.ksone.sf(ks_d_stat, n_tail))

def _fit_powerlaw_tail_grid(x_all: np.ndarray,
                            xmin_grid: Sequence[float],
                            min_tail: int) -> Dict[str, float]:
    x = np.asarray(x_all, float)
    x = x[np.isfinite(x) & (x > 0)]
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

def _bootstrap_plausibility_pvalue_fast(x_all: np.ndarray,
                                        B: int,
                                        xmin: float,
                                        reestimate_alpha: bool,
                                        seed: int = 0) -> Tuple[float, float, int]:
    """Fast bootstrap: fixed xmin; alpha re-estimated (optional) per bootstrap."""
    rng = np.random.default_rng(seed)
    x = np.asarray(x_all, float)
    tail = x[np.isfinite(x) & (x >= xmin)]
    if tail.size == 0:
        return np.nan, np.nan, 0
    alpha_hat = _alpha_mle_continuous(tail, xmin)
    D_obs = _ksD_tail(tail, xmin, alpha_hat)
    if not np.isfinite(D_obs):
        return np.nan, np.nan, tail.size
    ge = 0
    for _ in range(B):
        sim = xmin * (1.0 - rng.uniform(0, 1, size=tail.size)) ** (-1.0 / (alpha_hat - 1.0))
        a = _alpha_mle_continuous(sim, xmin) if reestimate_alpha else alpha_hat
        Db = _ksD_tail(sim, xmin, a)
        if np.isfinite(Db) and Db >= D_obs:
            ge += 1
    return ge / float(B), D_obs, tail.size

# -------------------- Alternative models & AIC --------------------
def _ll_powerlaw_continuous(x: np.ndarray, xmin: float, alpha: float) -> float:
    return (x.size * (np.log(alpha - 1.0) - np.log(xmin))
            - alpha * np.sum(np.log(x / xmin)))

def _ll_trunc_exponential(x: np.ndarray, xmin: float) -> Tuple[float, float]:
    lam = 1.0 / np.mean(x - xmin)
    ll = x.size * np.log(lam) - lam * np.sum(x - xmin)
    return ll, lam

def _ll_trunc_lognormal(x: np.ndarray, xmin: float) -> Tuple[float, float, float]:
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

def _aic_weights(aics: List[float]) -> List[float]:
    aics = np.asarray(aics, float)
    amin = np.nanmin(aics)
    d = aics - amin
    w = np.exp(-0.5 * d)
    w /= np.nansum(w)
    return list(w)

# -------------------- Degrees & diagnostics --------------------
def _concat_finite(arrs: Iterable[np.ndarray]) -> np.ndarray:
    xs = []
    for a in arrs:
        a = np.asarray(a, float)
        if a.size:
            xs.append(a[np.isfinite(a)])
    return np.concatenate(xs) if xs else np.array([], dtype=float)

def _degrees_from_block(A: np.ndarray) -> np.ndarray:
    if A.size == 0: return np.array([], float)
    Z = A.copy(); Z[~np.isfinite(Z)] = 0.0
    row_k = np.sum(Z, axis=1)
    col_k = np.sum(Z, axis=0)
    return _concat_finite([row_k, col_k])

def _density_from_bins(x: np.ndarray, bins: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    hist, edges = np.histogram(x, bins=bins)
    widths = np.diff(edges)
    mids = 0.5*(edges[:-1] + edges[1:])
    dens = hist / (np.sum(hist) * widths + 1e-12)
    return mids, dens

def _linear_regression_R2(x: np.ndarray, y: np.ndarray) -> float:
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
    lin_bins = np.linspace(np.min(k), np.max(k), bins_linear + 1)
    xm_lin, d_lin = _density_from_bins(k, lin_bins)
    kmin = np.min(k[k>0]) if np.any(k>0) else 1e-9
    log_bins = np.logspace(np.log10(kmin), np.log10(np.max(k)), bins_log + 1)
    xm_log, d_log = _density_from_bins(k, log_bins)
    R2_linear = _linear_regression_R2(xm_lin, d_lin)
    m = (xm_log > 0) & (d_log > 0)
    R2_loglog = _linear_regression_R2(np.log10(xm_log[m]), np.log10(d_log[m]))
    m2 = d_lin > 0
    R2_semilog = _linear_regression_R2(xm_lin[m2], -np.log10(d_lin[m2]))
    return dict(R2_linear=R2_linear, R2_loglog=R2_loglog, R2_semilog=R2_semilog)

def _subplots_block(v: np.ndarray, kvec: np.ndarray, title: str) -> go.Figure:
    fig = make_subplots(rows=2, cols=3, subplot_titles=(
        "Histogram(v)", "CCDF(v)", "CCDF(v) log–log",
        "k: linear bin (-log10 dens)", "k: log bin (R²)", "k: semilog (R²)"
    ))
    # histogram
    fig.add_trace(go.Histogram(x=v, nbinsx=60, name="v"), row=1, col=1)
    # CCDF linear
    xs, cc = _ccdf(v)
    fig.add_trace(go.Scatter(x=xs, y=cc, mode="lines", name="CCDF"), row=1, col=2)
    # CCDF log–log
    mpos = xs > 0
    fig.add_trace(go.Scatter(x=xs[mpos], y=cc[mpos], mode="lines", name="CCDF"), row=1, col=3)
    fig.update_xaxes(type="log", row=1, col=3); fig.update_yaxes(type="log", row=1, col=3)
    # k diagnostics
    R2s = _binned_R2s_for_k(kvec)
    kpos = kvec[kvec > 0]
    if kpos.size >= 12:
        # Linear bins with -log10 density on y-axis
        lin_bins = np.linspace(np.min(kpos), np.max(kpos), 12 + 1)
        xm, d = _density_from_bins(kpos, lin_bins)
        d_finite = d[d > 0]  # Only positive densities for log transform
        xm_finite = xm[d > 0]
        if len(d_finite) > 0:
            fig.add_trace(go.Scatter(x=xm_finite, y=-np.log10(d_finite), mode="markers", name="-log10 dens"), row=2, col=1)
            if np.isfinite(R2s["R2_linear"]) and len(xm_finite) > 2:
                X = np.vstack([np.ones_like(xm_finite), xm_finite]).T
                beta, *_ = np.linalg.lstsq(X, -np.log10(d_finite), rcond=None); yhat = X @ beta
                fig.add_trace(go.Scatter(x=xm_finite, y=yhat, mode="lines", name="lin fit"), row=2, col=1)
        
        # Log bins (unchanged)
        log_bins = np.logspace(np.log10(np.min(kpos)), np.log10(np.max(kpos)), 12 + 1)
        xm2, d2 = _density_from_bins(kpos, log_bins)
        m = (xm2 > 0) & (d2 > 0)
        fig.add_trace(go.Scatter(x=np.log10(xm2[m]), y=np.log10(d2[m]), mode="markers", name="log bins"), row=2, col=2)
        if np.isfinite(R2s["R2_loglog"]):
            X = np.vstack([np.ones(m.sum()), np.log10(xm2[m])]).T
            beta, *_ = np.linalg.lstsq(X, np.log10(d2[m]), rcond=None); yhat = X @ beta
            fig.add_trace(go.Scatter(x=np.log10(xm2[m]), y=yhat, mode="lines", name="log fit"), row=2, col=2)
        
        # Semilog (unchanged)
        xm3, d3 = _density_from_bins(kpos, lin_bins)
        m3 = d3 > 0
        fig.add_trace(go.Scatter(x=xm3[m3], y=-np.log10(d3[m3]), mode="markers", name="-log10 dens"), row=2, col=3)
        if np.isfinite(R2s["R2_semilog"]):
            X = np.vstack([np.ones(m3.sum()), xm3[m3]]).T
            beta, *_ = np.linalg.lstsq(X, -np.log10(d3[m3]), rcond=None); yhat = X @ beta
            fig.add_trace(go.Scatter(x=xm3[m3], y=yhat, mode="lines", name="semilog fit"), row=2, col=3)
    fig.update_layout(title=title, width=1100, height=760,
                      margin=dict(l=60, r=20, t=70, b=60), legend=dict(orientation="h"))
    return fig

def _image_flowable(png_path: str, max_width: float) -> RLImage:
    ir = ImageReader(png_path); iw, ih = ir.getSize()
    scale = min(1.0, max_width / float(iw))
    return RLImage(png_path, width=iw*scale, height=ih*scale)

def _display_name_for_csv(path: str,
                          display_names: Optional[Dict[str, str]] = None,
                          baseline=False) -> str:
    if display_names:
        base = os.path.basename(path)
        if path in display_names: return display_names[path]
        if base in display_names: return display_names[base]
    return "baseline (file #1)" if baseline else os.path.basename(path)

# -------------------- Data classes --------------------
@dataclass
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
    ks_pvalue: float  # Add KS p-value
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
    # Add raw data for visualization
    v_data: np.ndarray = None
    k_data: np.ndarray = None

@dataclass
class AggGoF:
    file_label: str
    kind: str  # "TS", "CT", "ALL"
    n_finite: int
    tail_n: int
    alpha: float
    xmin: float
    ks_D_tail: float
    ks_pvalue: float  # Add KS p-value
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

# -------------------- Fast analysis core --------------------
def _pairs_auto(csv_paths: Sequence[str], pairs: Optional[Sequence[Tuple[str,str]]]) -> List[Tuple[str,str]]:
    if pairs is None:
        norm_sets = []
        for p in csv_paths:
            _, pref2cols = _scan_prefix_to_cols(p)
            norm_sets.append(set(_normalize(x) for x in pref2cols.keys()))
        common_norm = set.intersection(*norm_sets) if norm_sets else set()
        prefixes_global = sorted(common_norm)
        requested_pairs_global = [(a, b) for i, a in enumerate(prefixes_global) for b in prefixes_global[i:]]
        return requested_pairs_global
    else:
        return [(_normalize(a), _normalize(b)) for (a, b) in pairs]

def _downsample_values(x: np.ndarray, max_cells: int, seed: int) -> np.ndarray:
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    n = x.size
    if max_cells and n > max_cells:
        rng = np.random.default_rng(seed)
        idx = rng.choice(n, size=max_cells, replace=False)
        return x[idx]
    return x

def _fit_block_fast(v: np.ndarray,
                    kvec: np.ndarray,
                    xmin_quantiles: Tuple[float, float, int],
                    B_bootstrap: int,
                    min_tail: int,
                    fast_bootstrap: bool,
                    seed: int) -> Dict[str, float]:
    # choose xmin grid
    x = v[v > 0]
    if x.size == 0:
        return dict(n_finite=0, tail_n=0, alpha=np.nan, xmin=np.nan, ks_D_tail=np.nan,
                    p_plaus=np.nan, ll_pl=-np.inf, aic_pl=np.inf,
                    ll_exp=-np.inf, aic_exp=np.inf, ll_logn=-np.inf, aic_logn=np.inf,
                    dAIC_exp=np.nan, dAIC_logn=np.nan, aicw_pl=np.nan, aicw_exp=np.nan, aicw_logn=np.nan,
                    R2_k_linear=np.nan, R2_k_loglog=np.nan, R2_k_semilog=np.nan)
    qlo, qhi, qn = xmin_quantiles
    qs = np.linspace(qlo, qhi, qn)
    xmin_grid = np.unique(np.quantile(x, qs))
    fit = _fit_powerlaw_tail_grid(x, xmin_grid, min_tail=min_tail)
    alpha, xmin = fit["alpha"], fit["xmin"]
    if not np.isfinite(alpha) or not np.isfinite(xmin) or fit["n_tail"] < min_tail:
        # still compute R² diagnostics on k
        R2s = _binned_R2s_for_k(kvec)
        return dict(n_finite=int(np.isfinite(v).size), tail_n=int(fit["n_tail"]), alpha=float(alpha),
                    xmin=float(xmin), ks_D_tail=float(fit["ks_D"]),
                    p_plaus=np.nan, ll_pl=-np.inf, aic_pl=np.inf,
                    ll_exp=-np.inf, aic_exp=np.inf, ll_logn=-np.inf, aic_logn=np.inf,
                    dAIC_exp=np.nan, dAIC_logn=np.nan, aicw_pl=np.nan, aicw_exp=np.nan, aicw_logn=np.nan,
                    R2_k_linear=float(R2s["R2_linear"]), R2_k_loglog=float(R2s["R2_loglog"]), R2_k_semilog=float(R2s["R2_semilog"]))
    # fast bootstrap p-value
    p_plaus, ksD, tail_n = _bootstrap_plausibility_pvalue_fast(x, B=B_bootstrap, xmin=xmin,
                                                               reestimate_alpha=True if fast_bootstrap else False,
                                                               seed=seed)
    tail = x[x >= xmin]
    ll_pl = _ll_powerlaw_continuous(tail, xmin, _alpha_mle_continuous(tail, xmin))
    ll_exp, _ = _ll_trunc_exponential(tail, xmin)
    ll_logn, _, _ = _ll_trunc_lognormal(tail, xmin)
    aic_pl  = _aic(ll_pl, 2)
    aic_exp = _aic(ll_exp, 2)
    aic_lgn = _aic(ll_logn, 3)
    dAIC_exp  = aic_exp - aic_pl
    dAIC_lgn  = aic_lgn - aic_pl
    w_pl, w_exp, w_lgn = _aic_weights([aic_pl, aic_exp, aic_lgn])
    R2s = _binned_R2s_for_k(kvec)
    return dict(
        n_finite=int(np.isfinite(v).sum()),
        tail_n=int(tail_n),
        alpha=float(alpha),
        xmin=float(xmin),
        ks_D_tail=float(ksD),
        p_plaus=float(p_plaus),
        ll_pl=float(ll_pl), aic_pl=float(aic_pl),
        ll_exp=float(ll_exp), aic_exp=float(aic_exp),
        ll_logn=float(ll_logn), aic_logn=float(aic_lgn),
        dAIC_exp=float(dAIC_exp), dAIC_logn=float(dAIC_lgn),
        aicw_pl=float(w_pl), aicw_exp=float(w_exp), aicw_logn=float(w_lgn),
        R2_k_linear=float(R2s["R2_linear"]),
        R2_k_loglog=float(R2s["R2_loglog"]),
        R2_k_semilog=float(R2s["R2_semilog"]),
    )

# -------------------- Rendering helpers --------------------
def _story_styles():
    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle(name="Small", fontSize=9, leading=11))
    return (styles["Title"], styles["Heading1"], styles["Heading2"],
            styles["BodyText"], styles["Small"])

def _render_per_file_report(file_label: str,
                            blocks: List[BlockGoF],
                            aggs: List[AggGoF],
                            figs: Dict[str, go.Figure],
                            out_pdf_path: Optional[str],  # Can be None to skip PDF
                            out_html_path: Optional[str]):
    title_style, h1, h2, body, small = _story_styles()
    page_w, page_h = A4; margins = 36; max_img_width = page_w - 2*margins
    story = []; html_sections: List[str] = []

    story.append(Paragraph(f"Power-law GoF Report — {file_label}", title_style))
    story.append(Spacer(1, 0.12*inch))
    story.append(Paragraph("Fast mode: fixed-xmin bootstrap, downsampling, integrated visuals only.", body))
    story.append(PageBreak())

    # Per-block table
    rows = [[
        "A","B","n_finite","tail_n","alpha","xmin","KS-D","p(PL)","ΔAIC(exp)","ΔAIC(logn)","AICw(PL)",
        "R²_lin","R²_loglog","R²_semilog"
    ]]
    for b in blocks:
        rows.append([
            b.region_a, b.region_b, f"{b.n_finite:,}", f"{b.tail_n:,}",
            f"{b.alpha:.4f}", f"{b.xmin:.4f}",
            f"{b.ks_D_tail:.3f}", f"{_format_pvalue(b.ks_pvalue)}", f"{_format_pvalue(b.p_plaus)}",
            f"{b.dAIC_exp:+.5f}", f"{b.dAIC_logn:+.5f}",
            f"{b.aicw_pl:.5f}",
            f"{b.R2_k_linear:.3f}", f"{b.R2_k_loglog:.3f}", f"{b.R2_k_semilog:.3f}",
        ])
    tbl = Table(rows, hAlign="LEFT")
    tbl.setStyle(TableStyle([
        ("FONTNAME", (0,0), (-1,0), "Helvetica-Bold"),
        ("BACKGROUND", (0,0), (-1,0), colors.lightgrey),
        ("GRID", (0,0), (-1,-1), 0.3, colors.grey),
        ("ALIGN", (2,1), (-1,-1), "RIGHT"),
    ]))
    story.append(Paragraph("Per-block GoF", h1)); story.append(tbl); story.append(PageBreak())
    if out_html_path:
        html_sections.append(f"<section class='section'><h2>{_html_escape(file_label)} — Per-block GoF</h2>")
        html_sections.append(_table_html(rows)); html_sections.append("</section><div class='divider'></div>")

    # Integrated figures - optimized for speed
    with tempfile.TemporaryDirectory() as tmpdir:
        for kind in ["TS","CT","ALL"]:
            print(f"         📊 Processing {kind} figure...")
            fig = figs[kind]
            png = os.path.join(tmpdir, f"{file_label}_{kind}_integrated.png").replace("/","_").replace("\\","_")
            
            # Optimize image generation for speed
            try:
                fig.write_image(png, scale=1.5, width=1200, height=800)  # Reduced scale and fixed size
            except Exception as e:
                print(f"         ⚠️  Image generation failed for {kind}: {e}")
                continue
                
            story.append(Paragraph(f"Integrated {kind}", h2))
            story.append(_image_flowable(png, max_img_width))
            story.append(Spacer(1, 0.08*inch))
            if out_html_path:
                html_sections.append(_wrap_fig_html(fig, title=f"{file_label} — Integrated {kind}"))

    # Aggregated table
    arows = [["Scope","n_finite","tail_n","alpha","xmin","KS-D","KS-p","p(PL)","ΔAIC(exp)","ΔAIC(logn)","AICw(PL)",
              "R²_lin","R²_loglog","R²_semilog"]]
    for a in aggs:
        arows.append([
            a.kind, f"{a.n_finite:,}", f"{a.tail_n:,}",
            f"{a.alpha:.4f}", f"{a.xmin:.4f}",
            f"{a.ks_D_tail:.3f}", f"{_format_pvalue(a.ks_pvalue)}", f"{_format_pvalue(a.p_plaus)}",
            f"{a.dAIC_exp:+.5f}", f"{a.dAIC_logn:+.5f}",
            f"{a.aicw_pl:.5f}",
            f"{a.R2_k_linear:.3f}", f"{a.R2_k_loglog:.3f}", f"{a.R2_k_semilog:.3f}",
        ])
    tbl2 = Table(arows, hAlign="LEFT")
    tbl2.setStyle(TableStyle([
        ("FONTNAME", (0,0), (-1,0), "Helvetica-Bold"),
        ("BACKGROUND", (0,0), (-1,0), colors.HexColor("#EFEFEF")),
        ("GRID", (0,0), (-1,-1), 0.3, colors.grey),
        ("ALIGN", (1,1), (-1,-1), "RIGHT"),
    ]))
    story.append(Paragraph("Integrated TS/CT/ALL — GoF", h1)); story.append(tbl2)

    # Generate PDF only if path provided
    if out_pdf_path:
        doc = SimpleDocTemplate(out_pdf_path, pagesize=A4,
                                leftMargin=36, rightMargin=36, topMargin=36, bottomMargin=36)
        doc.build(story)
    
    # Always try to generate HTML if path provided
    if out_html_path:
        _write_html_report(out_html_path, html_sections, f"GoF — {file_label}")

def _render_combined_tissues_report(all_gofs: list[BlockGoF], tissues: list[str], file_prefix: str) -> go.Figure:
    """Create combined report showing all tissues in same subplots"""
    colors = px.colors.qualitative.Plotly
    
    # Create combined plot
    fig = make_subplots(rows=2, cols=3, subplot_titles=(
        f"Histogram(v) - All {len(tissues)} Tissues", 
        f"CCDF(v) - All {len(tissues)} Tissues", 
        f"CCDF(v) log–log - All {len(tissues)} Tissues",
        f"k: linear bin (-log10 dens) - All {len(tissues)} Tissues", 
        f"k: log bin (R²) - All {len(tissues)} Tissues", 
        f"k: semilog (R²) - All {len(tissues)} Tissues"
    ))
    
    # Plot each tissue with different color
    for i, (gof, tissue) in enumerate(zip(all_gofs, tissues)):
        color = colors[i % len(colors)]
        v, kvec = gof.v_data, gof.k_data
        
        # Histograms
        fig.add_trace(go.Histogram(x=v, nbinsx=60, name=f"{tissue}", 
                                   marker_color=color, opacity=0.7), row=1, col=1)
        
        # CCDF linear
        xs, cc = _ccdf(v)
        fig.add_trace(go.Scatter(x=xs, y=cc, mode="lines", name=f"{tissue} CCDF",
                                line=dict(color=color)), row=1, col=2)
        
        # CCDF log–log
        mpos = xs > 0
        fig.add_trace(go.Scatter(x=xs[mpos], y=cc[mpos], mode="lines", name=f"{tissue} CCDF",
                                line=dict(color=color)), row=1, col=3)
        
        # k diagnostics
        kpos = kvec[kvec > 0]
        if kpos.size >= 12:
            # Linear bins with -log10 density
            lin_bins = np.linspace(np.min(kpos), np.max(kpos), 12 + 1)
            xm, d = _density_from_bins(kpos, lin_bins)
            d_finite = d[d > 0]
            xm_finite = xm[d > 0]
            if len(d_finite) > 0:
                fig.add_trace(go.Scatter(x=xm_finite, y=-np.log10(d_finite), mode="markers", 
                                        name=f"{tissue} -log10 dens", marker=dict(color=color)), row=2, col=1)
            
            # Log bins
            log_bins = np.logspace(np.log10(np.min(kpos)), np.log10(np.max(kpos)), 12 + 1)
            xm2, d2 = _density_from_bins(kpos, log_bins)
            m = (xm2 > 0) & (d2 > 0)
            fig.add_trace(go.Scatter(x=np.log10(xm2[m]), y=np.log10(d2[m]), mode="markers", 
                                    name=f"{tissue} log bins", marker=dict(color=color)), row=2, col=2)
            
            # Semilog
            xm3, d3 = _density_from_bins(kpos, lin_bins)
            m3 = d3 > 0
            fig.add_trace(go.Scatter(x=xm3[m3], y=-np.log10(d3[m3]), mode="markers", 
                                    name=f"{tissue} -log10 dens", marker=dict(color=color)), row=2, col=3)
    
    # Update log axes
    fig.update_xaxes(type="log", row=1, col=3)
    fig.update_yaxes(type="log", row=1, col=3)
    
    fig.update_layout(title=f"Combined Tissues Analysis - {file_prefix}", 
                      width=1400, height=900,
                      margin=dict(l=60, r=20, t=100, b=60), 
                      legend=dict(orientation="h", y=-0.1))
    
    return fig

def _render_ts_vs_ct_comparison(ts_gofs: list[BlockGoF], ct_gofs: list[BlockGoF], 
                               ts_tissues: list[str], ct_tissues: list[str]) -> go.Figure:
    """Create comparison report showing all TS vs all CT"""
    
    fig = make_subplots(rows=2, cols=3, subplot_titles=(
        "Histogram(v) - TS vs CT", "CCDF(v) - TS vs CT", "CCDF(v) log–log - TS vs CT",
        "k: linear bin (-log10 dens) - TS vs CT", "k: log bin (R²) - TS vs CT", "k: semilog (R²) - TS vs CT"
    ))
    
    # Plot TS data in blue tones
    for i, (gof, tissue) in enumerate(zip(ts_gofs, ts_tissues)):
        v, kvec = gof.v_data, gof.k_data
        color = f"rgb(0, {100 + i*30}, 255)"  # Blue gradient
        
        # Histograms
        fig.add_trace(go.Histogram(x=v, nbinsx=60, name=f"TS-{tissue}", 
                                   marker_color=color, opacity=0.6), row=1, col=1)
        
        # CCDF linear
        xs, cc = _ccdf(v)
        fig.add_trace(go.Scatter(x=xs, y=cc, mode="lines", name=f"TS-{tissue}",
                                line=dict(color=color, dash='solid')), row=1, col=2)
        
        # CCDF log–log
        mpos = xs > 0
        fig.add_trace(go.Scatter(x=xs[mpos], y=cc[mpos], mode="lines", name=f"TS-{tissue}",
                                line=dict(color=color, dash='solid')), row=1, col=3)
        
        # k diagnostics
        kpos = kvec[kvec > 0]
        if kpos.size >= 12:
            lin_bins = np.linspace(np.min(kpos), np.max(kpos), 12 + 1)
            xm, d = _density_from_bins(kpos, lin_bins)
            d_finite = d[d > 0]
            xm_finite = xm[d > 0]
            if len(d_finite) > 0:
                fig.add_trace(go.Scatter(x=xm_finite, y=-np.log10(d_finite), mode="markers", 
                                        name=f"TS-{tissue}", marker=dict(color=color, symbol='circle')), row=2, col=1)
            
            log_bins = np.logspace(np.log10(np.min(kpos)), np.log10(np.max(kpos)), 12 + 1)
            xm2, d2 = _density_from_bins(kpos, log_bins)
            m = (xm2 > 0) & (d2 > 0)
            fig.add_trace(go.Scatter(x=np.log10(xm2[m]), y=np.log10(d2[m]), mode="markers", 
                                    name=f"TS-{tissue}", marker=dict(color=color, symbol='circle')), row=2, col=2)
            
            xm3, d3 = _density_from_bins(kpos, lin_bins)
            m3 = d3 > 0
            fig.add_trace(go.Scatter(x=xm3[m3], y=-np.log10(d3[m3]), mode="markers", 
                                    name=f"TS-{tissue}", marker=dict(color=color, symbol='circle')), row=2, col=3)
    
    # Plot CT data in red tones
    for i, (gof, tissue) in enumerate(zip(ct_gofs, ct_tissues)):
        v, kvec = gof.v_data, gof.k_data
        color = f"rgb(255, {100 + i*30}, 0)"  # Red gradient
        
        # Histograms
        fig.add_trace(go.Histogram(x=v, nbinsx=60, name=f"CT-{tissue}", 
                                   marker_color=color, opacity=0.6), row=1, col=1)
        
        # CCDF linear
        xs, cc = _ccdf(v)
        fig.add_trace(go.Scatter(x=xs, y=cc, mode="lines", name=f"CT-{tissue}",
                                line=dict(color=color, dash='dash')), row=1, col=2)
        
        # CCDF log–log
        mpos = xs > 0
        fig.add_trace(go.Scatter(x=xs[mpos], y=cc[mpos], mode="lines", name=f"CT-{tissue}",
                                line=dict(color=color, dash='dash')), row=1, col=3)
        
        # k diagnostics
        kpos = kvec[kvec > 0]
        if kpos.size >= 12:
            lin_bins = np.linspace(np.min(kpos), np.max(kpos), 12 + 1)
            xm, d = _density_from_bins(kpos, lin_bins)
            d_finite = d[d > 0]
            xm_finite = xm[d > 0]
            if len(d_finite) > 0:
                fig.add_trace(go.Scatter(x=xm_finite, y=-np.log10(d_finite), mode="markers", 
                                        name=f"CT-{tissue}", marker=dict(color=color, symbol='diamond')), row=2, col=1)
            
            log_bins = np.logspace(np.log10(np.min(kpos)), np.log10(np.max(kpos)), 12 + 1)
            xm2, d2 = _density_from_bins(kpos, log_bins)
            m = (xm2 > 0) & (d2 > 0)
            fig.add_trace(go.Scatter(x=np.log10(xm2[m]), y=np.log10(d2[m]), mode="markers", 
                                    name=f"CT-{tissue}", marker=dict(color=color, symbol='diamond')), row=2, col=2)
            
            xm3, d3 = _density_from_bins(kpos, lin_bins)
            m3 = d3 > 0
            fig.add_trace(go.Scatter(x=xm3[m3], y=-np.log10(d3[m3]), mode="markers", 
                                    name=f"CT-{tissue}", marker=dict(color=color, symbol='diamond')), row=2, col=3)
    
    # Update log axes
    fig.update_xaxes(type="log", row=1, col=3)
    fig.update_yaxes(type="log", row=1, col=3)
    
    fig.update_layout(title="TS vs CT Tissue Comparison", 
                      width=1400, height=900,
                      margin=dict(l=60, r=20, t=100, b=60), 
                      legend=dict(orientation="h", y=-0.1))
    
    return fig

def _render_all_files_report(analyses: Dict[str, Dict], out_pdf_path: str, out_html_path: Optional[str]):
    title_style, h1, h2, body, small = _story_styles()
    story = []; html_sections: List[str] = []
    story.append(Paragraph("Power-law GoF — All Files (FAST)", title_style)); story.append(Spacer(1, 0.12*inch))
    story.append(PageBreak())
    labels = list(analyses.keys())
    for lab in labels:
        blocks: List[BlockGoF] = analyses[lab]["blocks"]
        aggs: List[AggGoF] = analyses[lab]["aggs"]
        story.append(Paragraph(lab, h1))
        rows = [[
            "A","B","n_finite","tail_n","alpha","xmin","KS-D","KS-p","p(PL)","ΔAIC(exp)","ΔAIC(logn)","AICw(PL)",
            "R²_lin","R²_loglog","R²_semilog"
        ]]
        for b in blocks:
            rows.append([
                b.region_a, b.region_b, f"{b.n_finite:,}", f"{b.tail_n:,}",
                f"{b.alpha:.4f}", f"{b.xmin:.4f}",
                f"{b.ks_D_tail:.3f}", f"{b.ks_pvalue:.3f}", f"{b.p_plaus:.3f}",
                f"{b.dAIC_exp:+.5f}", f"{b.dAIC_logn:+.5f}",
                f"{b.aicw_pl:.5f}",
                f"{b.R2_k_linear:.3f}", f"{b.R2_k_loglog:.3f}", f"{b.R2_k_semilog:.3f}",
            ])
        tbl = Table(rows, hAlign="LEFT")
        tbl.setStyle(TableStyle([
            ("FONTNAME", (0,0), (-1,0), "Helvetica-Bold"),
            ("BACKGROUND", (0,0), (-1,0), colors.lightgrey),
            ("GRID", (0,0), (-1,-1), 0.3, colors.grey),
            ("ALIGN", (2,1), (-1,-1), "RIGHT"),
        ]))
        story.append(tbl); story.append(PageBreak())

        arows = [["Scope","n_finite","tail_n","alpha","xmin","KS-D","KS-p","p(PL)","ΔAIC(exp)","ΔAIC(logn)","AICw(PL)",
                  "R²_lin","R²_loglog","R²_semilog"]]
        for a in aggs:
            arows.append([
                a.kind, f"{a.n_finite:,}", f"{a.tail_n:,}",
                f"{a.alpha:.4f}", f"{a.xmin:.4f}",
                f"{a.ks_D_tail:.3f}", f"{a.ks_pvalue:.3f}", f"{a.p_plaus:.3f}",
                f"{a.dAIC_exp:+.5f}", f"{a.dAIC_logn:+.5f}",
                f"{a.aicw_pl:.5f}",
                f"{a.R2_k_linear:.3f}", f"{a.R2_k_loglog:.3f}", f"{a.R2_k_semilog:.3f}",
            ])
        tbl2 = Table(arows, hAlign="LEFT")
        tbl2.setStyle(TableStyle([
            ("FONTNAME", (0,0), (-1,0), "Helvetica-Bold"),
            ("BACKGROUND", (0,0), (-1,0), colors.HexColor("#EFEFEF")),
            ("GRID", (0,0), (-1,-1), 0.3, colors.grey),
            ("ALIGN", (1,1), (-1,-1), "RIGHT"),
        ]))
        story.append(tbl2); story.append(PageBreak())
        
        # Generate HTML sections for this file
        html_sections.append(f'<section class="section"><h2>{_html_escape(lab)}</h2>')
        
        # Add blocks table to HTML
        html_rows = [["A","B","n_finite","tail_n","alpha","xmin","KS-D","KS-p","p(PL)","ΔAIC(exp)","ΔAIC(logn)","AICw(PL)",
                      "R²_lin","R²_loglog","R²_semilog"]]
        for b in blocks:
            html_rows.append([
                b.region_a, b.region_b, f"{b.n_finite:,}", f"{b.tail_n:,}",
                f"{b.alpha:.4f}", f"{b.xmin:.4f}",
                f"{b.ks_D_tail:.3f}", f"{b.ks_pvalue:.3f}", f"{b.p_plaus:.3f}",
                f"{b.dAIC_exp:+.5f}", f"{b.dAIC_logn:+.5f}",
                f"{b.aicw_pl:.5f}",
                f"{b.R2_k_linear:.3f}", f"{b.R2_k_loglog:.3f}", f"{b.R2_k_semilog:.3f}",
            ])
        html_sections.append(_table_html(html_rows, "Block-level analysis"))
        
        # Add aggregations table to HTML
        html_agg_rows = [["Scope","n_finite","tail_n","alpha","xmin","KS-D","KS-p","p(PL)","ΔAIC(exp)","ΔAIC(logn)","AICw(PL)",
                          "R²_lin","R²_loglog","R²_semilog"]]
        for a in aggs:
            html_agg_rows.append([
                a.kind, f"{a.n_finite:,}", f"{a.tail_n:,}",
                f"{a.alpha:.4f}", f"{a.xmin:.4f}",
                f"{a.ks_D_tail:.3f}", f"{a.ks_pvalue:.3f}", f"{a.p_plaus:.3f}",
                f"{a.dAIC_exp:+.5f}", f"{a.dAIC_logn:+.5f}",
                f"{a.aicw_pl:.5f}",
                f"{a.R2_k_linear:.3f}", f"{a.R2_k_loglog:.3f}", f"{a.R2_k_semilog:.3f}",
            ])
        html_sections.append(_table_html(html_agg_rows, "Aggregated analysis"))
        
        # Add figures if available
        if lab in analyses and "figs" in analyses[lab]:
            figs = analyses[lab]["figs"]
            for fig_name, fig in figs.items():
                html_sections.append(f'<h3>{_html_escape(fig_name)}</h3>')
                html_sections.append(_wrap_fig_html(fig))
        
        html_sections.append('</section>')

    doc = SimpleDocTemplate(out_pdf_path, pagesize=A4,
                            leftMargin=36, rightMargin=36, topMargin=36, bottomMargin=36)
    doc.build(story)
    if out_html_path:
        _write_html_report(out_html_path, html_sections, "GoF — All Files (FAST)")

# -------------------- Statistical Comparison Functions --------------------
def _compare_tissues_within_file(blocks: List[BlockGoF], metric: str = "alpha") -> Dict[str, float]:
    """Compare tissues/regions within same file using Kruskal-Wallis test"""
    # Group blocks by region pairs
    tissue_data = {}
    for block in blocks:
        tissue_pair = f"{block.region_a}×{block.region_b}"
        if tissue_pair not in tissue_data:
            tissue_data[tissue_pair] = []
        tissue_data[tissue_pair].append(getattr(block, metric))
    
    if len(tissue_data) < 2:
        return {"kruskal_p": np.nan, "n_groups": len(tissue_data)}
    
    # Perform Kruskal-Wallis test
    values = [np.array(data) for data in tissue_data.values()]
    try:
        h_stat, p_val = kruskal(*values)
        return {
            "kruskal_p": p_val,
            "kruskal_h": h_stat,
            "n_groups": len(tissue_data),
            "group_sizes": [len(v) for v in values]
        }
    except Exception as e:
        return {"kruskal_p": np.nan, "n_groups": len(tissue_data), "error": str(e)}

def _compare_files_pairwise(blocks_file1: List[BlockGoF], blocks_file2: List[BlockGoF], 
                           metric: str = "alpha") -> Dict[str, float]:
    """Compare same tissue pairs between two files using Mann-Whitney U test"""
    # Match tissue pairs between files
    pairs1 = {f"{b.region_a}×{b.region_b}": getattr(b, metric) for b in blocks_file1}
    pairs2 = {f"{b.region_a}×{b.region_b}": getattr(b, metric) for b in blocks_file2}
    
    # Find common pairs
    common_pairs = set(pairs1.keys()) & set(pairs2.keys())
    if len(common_pairs) < 3:
        return {"mannwhitney_p": np.nan, "n_pairs": len(common_pairs)}
    
    # Get values for common pairs
    values1 = [pairs1[pair] for pair in common_pairs]
    values2 = [pairs2[pair] for pair in common_pairs]
    
    try:
        stat, p_val = mannwhitneyu(values1, values2, alternative='two-sided')
        return {
            "mannwhitney_p": p_val,
            "mannwhitney_u": stat,
            "n_pairs": len(common_pairs),
            "mean1": np.mean(values1),
            "mean2": np.mean(values2),
            "median1": np.median(values1),
            "median2": np.median(values2)
        }
    except Exception as e:
        return {"mannwhitney_p": np.nan, "n_pairs": len(common_pairs), "error": str(e)}

def _compare_pairs_within_tissues(blocks: List[BlockGoF], metric: str = "alpha") -> Dict[str, Dict[str, float]]:
    """Compare different pairs within same tissues (e.g., A×B vs A×C vs B×C)"""
    # Group by individual tissues/regions
    tissue_pairs = {}
    for block in blocks:
        tissues = sorted([block.region_a, block.region_b])  # Normalize order
        tissue_key = f"{tissues[0]}+{tissues[1]}"
        if tissue_key not in tissue_pairs:
            tissue_pairs[tissue_key] = []
        tissue_pairs[tissue_key].append(getattr(block, metric))
    
    results = {}
    tissue_keys = list(tissue_pairs.keys())
    
    # Pairwise comparisons between tissue combinations
    for i in range(len(tissue_keys)):
        for j in range(i+1, len(tissue_keys)):
            key1, key2 = tissue_keys[i], tissue_keys[j]
            values1, values2 = tissue_pairs[key1], tissue_pairs[key2]
            
            if len(values1) >= 3 and len(values2) >= 3:
                try:
                    stat, p_val = mannwhitneyu(values1, values2, alternative='two-sided')
                    results[f"{key1}_vs_{key2}"] = {
                        "mannwhitney_p": p_val,
                        "mannwhitney_u": stat,
                        "n1": len(values1),
                        "n2": len(values2),
                        "mean1": np.mean(values1),
                        "mean2": np.mean(values2)
                    }
                except Exception as e:
                    results[f"{key1}_vs_{key2}"] = {
                        "mannwhitney_p": np.nan,
                        "error": str(e)
                    }
    
    return results

def _generate_comparison_tables(all_blocks: List[BlockGoF], 
                               analyses: Dict[str, Dict]) -> Tuple[List[List[str]], List[List[str]], List[List[str]]]:
    """Generate statistical comparison tables for within-file, between-file, and within-tissue comparisons"""
    
    # Within-file tissue comparisons
    within_file_rows = [["File", "Metric", "Kruskal-Wallis p", "H-statistic", "N groups", "Group sizes"]]
    
    for label, content in analyses.items():
        blocks = content["blocks"]
        for metric in ["alpha", "ks_D_tail", "aicw_pl", "ks_pvalue"]:
            comparison = _compare_tissues_within_file(blocks, metric)
            within_file_rows.append([
                label,
                metric,
                f"{comparison.get('kruskal_p', np.nan):.3g}",
                f"{comparison.get('kruskal_h', np.nan):.3f}",
                str(comparison.get('n_groups', 0)),
                str(comparison.get('group_sizes', []))
            ])
    
    # Between-file comparisons
    between_file_rows = [["File 1", "File 2", "Metric", "Mann-Whitney p", "U-statistic", "N pairs", "Mean1", "Mean2"]]
    
    file_labels = list(analyses.keys())
    for i in range(len(file_labels)):
        for j in range(i+1, len(file_labels)):
            label1, label2 = file_labels[i], file_labels[j]
            blocks1 = analyses[label1]["blocks"]
            blocks2 = analyses[label2]["blocks"]
            
            for metric in ["alpha", "ks_D_tail", "aicw_pl", "ks_pvalue"]:
                comparison = _compare_files_pairwise(blocks1, blocks2, metric)
                between_file_rows.append([
                    label1,
                    label2,
                    metric,
                    f"{comparison.get('mannwhitney_p', np.nan):.3g}",
                    f"{comparison.get('mannwhitney_u', np.nan):.3f}",
                    str(comparison.get('n_pairs', 0)),
                    f"{comparison.get('mean1', np.nan):.4f}",
                    f"{comparison.get('mean2', np.nan):.4f}"
                ])
    
    # Within-tissue pair comparisons
    within_tissue_rows = [["File", "Comparison", "Metric", "Mann-Whitney p", "U-statistic", "N1", "N2", "Mean1", "Mean2"]]
    
    for label, content in analyses.items():
        blocks = content["blocks"]
        for metric in ["alpha", "ks_D_tail", "aicw_pl", "ks_pvalue"]:
            pair_comparisons = _compare_pairs_within_tissues(blocks, metric)
            for comparison_name, stats in pair_comparisons.items():
                within_tissue_rows.append([
                    label,
                    comparison_name,
                    metric,
                    f"{stats.get('mannwhitney_p', np.nan):.3g}",
                    f"{stats.get('mannwhitney_u', np.nan):.3f}",
                    str(stats.get('n1', 0)),
                    str(stats.get('n2', 0)),
                    f"{stats.get('mean1', np.nan):.4f}",
                    f"{stats.get('mean2', np.nan):.4f}"
                ])
    
    return within_file_rows, between_file_rows, within_tissue_rows

def _render_comparisons_only(blocks_all: List[BlockGoF],
                             csv_paths: Sequence[str],
                             file_display_names: Optional[Dict[str,str]],
                             out_pdf_path: str,
                             out_html_path: Optional[str],
                             analyses: Optional[Dict[str, Dict]] = None):
    title_style, h1, h2, body, small = _story_styles()
    story=[]; html_sections=[]
    story.append(Paragraph("Comparisons Only — Rank Preservation (FAST)", title_style))
    story.append(Spacer(1, 0.12*inch))

    df = pd.DataFrame([b.__dict__ for b in blocks_all])
    df["pair"] = df["region_a"] + "×" + df["region_b"]
    labels = [ _display_name_for_csv(p, file_display_names, baseline=(i==0)) for i,p in enumerate(csv_paths) ]
    counts = df.groupby("pair")["file_label"].nunique()
    keep_pairs = set(counts[counts == len(labels)].index)
    df = df[df["pair"].isin(keep_pairs)].copy()

    rows_rank = [["Metric","File A","File B","Spearman ρ","p","n pairs"]]
    for metric, asc in [("ks_D_tail", True), ("aicw_pl", False)]:
        for i in range(len(labels)):
            for j in range(i+1, len(labels)):
                a, b = labels[i], labels[j]
                A = df[df.file_label==a].set_index("pair")[metric].rank(ascending=asc, method="average")
                B = df[df.file_label==b].set_index("pair")[metric].rank(ascending=asc, method="average")
                idx = A.index.intersection(B.index)
                rho, pval = sstats.spearmanr(A.loc[idx], B.loc[idx])
                rows_rank.append([metric, a, b, f"{rho:.3f}", f"{pval:.3g}", len(idx)])

    tblR = Table(rows_rank, hAlign="LEFT")
    tblR.setStyle(TableStyle([
        ("FONTNAME",(0,0),(-1,0),"Helvetica-Bold"),
        ("BACKGROUND",(0,0),(-1,0),colors.lightgrey),
        ("GRID",(0,0),(-1,-1),0.3,colors.grey),
        ("ALIGN",(3,1),(-1,-1),"RIGHT"),
    ]))
    story.append(Paragraph("Rank-order preservation across files", h1)); story.append(tblR)
    
    # Add new statistical comparison tables if analyses provided
    if analyses:
        story.append(PageBreak())
        story.append(Paragraph("Statistical Comparisons of Power-law Parameters", h1))
        story.append(Spacer(1, 0.1*inch))
        
        within_file_rows, between_file_rows, within_tissue_rows = _generate_comparison_tables(blocks_all, analyses)
        
        # Within-file tissue comparisons
        story.append(Paragraph("Within-file tissue comparisons (Kruskal-Wallis)", h2))
        story.append(Spacer(1, 0.05*inch))
        tbl_within_file = Table(within_file_rows, hAlign="LEFT")
        tbl_within_file.setStyle(TableStyle([
            ("FONTNAME",(0,0),(-1,0),"Helvetica-Bold"),
            ("BACKGROUND",(0,0),(-1,0),colors.lightblue),
            ("GRID",(0,0),(-1,-1),0.3,colors.grey),
            ("ALIGN",(2,1),(-1,-1),"RIGHT"),
        ]))
        story.append(tbl_within_file)
        story.append(Spacer(1, 0.1*inch))
        
        # Between-file comparisons
        story.append(Paragraph("Between-file comparisons (Mann-Whitney U)", h2))
        story.append(Spacer(1, 0.05*inch))
        tbl_between_file = Table(between_file_rows, hAlign="LEFT")
        tbl_between_file.setStyle(TableStyle([
            ("FONTNAME",(0,0),(-1,0),"Helvetica-Bold"),
            ("BACKGROUND",(0,0),(-1,0),colors.lightgreen),
            ("GRID",(0,0),(-1,-1),0.3,colors.grey),
            ("ALIGN",(3,1),(-1,-1),"RIGHT"),
        ]))
        story.append(tbl_between_file)
        story.append(Spacer(1, 0.1*inch))
        
        # Within-tissue pair comparisons
        story.append(Paragraph("Within-tissue pair comparisons (Mann-Whitney U)", h2))
        story.append(Spacer(1, 0.05*inch))
        tbl_within_tissue = Table(within_tissue_rows, hAlign="LEFT")
        tbl_within_tissue.setStyle(TableStyle([
            ("FONTNAME",(0,0),(-1,0),"Helvetica-Bold"),
            ("BACKGROUND",(0,0),(-1,0),colors.lightyellow),
            ("GRID",(0,0),(-1,-1),0.3,colors.grey),
            ("ALIGN",(3,1),(-1,-1),"RIGHT"),
        ]))
        story.append(tbl_within_tissue)
        
        # Prepare HTML sections
        html_sections.append(_table_html(rows_rank, "Rank-order preservation"))
        html_sections.append(_table_html(within_file_rows, "Within-file tissue comparisons (Kruskal-Wallis)"))
        html_sections.append(_table_html(between_file_rows, "Between-file comparisons (Mann-Whitney U)"))
        html_sections.append(_table_html(within_tissue_rows, "Within-tissue pair comparisons (Mann-Whitney U)"))
    else:
        html_sections.append(_table_html(rows_rank, "Rank-order preservation"))
    
    doc = SimpleDocTemplate(out_pdf_path, pagesize=A4,
                            leftMargin=36, rightMargin=36, topMargin=36, bottomMargin=36)
    doc.build(story)
    if out_html_path:
        _write_html_report(out_html_path, html_sections, "Comparisons Only (FAST)")

# -------------------- Module-level functions for multiprocessing --------------------
def _analyze_pair_task(args):
    """Task function for multiprocessing in make_powerlaw_gof_multi_reports_fast."""
    (a_n, b_n), (v, k), xmin_quantiles, B_bootstrap, min_tail, fast, seed0 = args
    ana = _fit_block_fast(v, k, xmin_quantiles, B_bootstrap, min_tail,
                          fast_bootstrap=fast, seed=seed0)
    # Include raw data in the analysis result
    ana["v_data"] = v
    ana["k_data"] = k
    return (a_n, b_n, ana)

# -------------------- Public API --------------------
def _create_html_index(out_dir: str, analyses: Dict[str, Dict], report_count: int):
    """Generate an HTML index page for navigating all generated reports."""
    import glob
    from datetime import datetime
    
    # Discover all HTML files in output directory
    html_files = glob.glob(os.path.join(out_dir, "*.html"))
    html_files = [os.path.basename(f) for f in html_files]
    
    # Categorize files
    per_file_reports = [f for f in html_files if f.startswith("GoF_FAST_") and not any(x in f for x in ["ALL_FILES", "COMPARISONS", "COMBINED", "TS_vs_CT", "ENHANCED", "KS_", "R2_", "TEST_"])]
    combined_reports = [f for f in html_files if "ALL_FILES" in f or "COMPARISONS" in f]
    enhanced_reports = [f for f in html_files if "COMBINED" in f or "TS_vs_CT" in f or "ENHANCED" in f or "R2_ANALYSIS" in f]
    ks_reports = [f for f in html_files if "KS_" in f and "VISUALIZATION" in f]
    test_reports = [f for f in html_files if f.startswith("TEST_")]
    individual_ks = [f for f in html_files if "KS_" in f and "VISUALIZATION" not in f]
    
    # Discover Barabási reports in subdirectory
    barabasi_reports = []
    barabasi_dir = os.path.join(out_dir, "barabasi_analysis")
    if os.path.exists(barabasi_dir):
        barabasi_files = glob.glob(os.path.join(barabasi_dir, "*.html"))
        barabasi_reports = [f"barabasi_analysis/{os.path.basename(f)}" for f in barabasi_files]
    
    other_reports = [f for f in html_files if f not in per_file_reports + combined_reports + enhanced_reports + ks_reports + test_reports + individual_ks]
    
    current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    
    index_html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Power-Law GoF Analysis Reports - Navigation Index</title>
    <style>
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
            margin: 0;
            padding: 20px;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            min-height: 100vh;
        }}
        .container {{
            max-width: 1200px;
            margin: 0 auto;
            background: white;
            border-radius: 15px;
            box-shadow: 0 10px 30px rgba(0,0,0,0.2);
            overflow: hidden;
        }}
        .header {{
            background: linear-gradient(90deg, #2c3e50, #34495e);
            color: white;
            padding: 30px;
            text-align: center;
        }}
        .header h1 {{
            margin: 0 0 10px 0;
            font-size: 2.5em;
            font-weight: 300;
        }}
        .header p {{
            margin: 0;
            opacity: 0.8;
            font-size: 1.1em;
        }}
        .content {{
            padding: 30px;
        }}
        .stats {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px;
            margin-bottom: 30px;
        }}
        .stat-card {{
            background: #f8f9fa;
            padding: 20px;
            border-radius: 10px;
            text-align: center;
            border-left: 4px solid #3498db;
        }}
        .stat-number {{
            font-size: 2em;
            font-weight: bold;
            color: #2c3e50;
            margin-bottom: 5px;
        }}
        .stat-label {{
            color: #7f8c8d;
            font-size: 0.9em;
            text-transform: uppercase;
            letter-spacing: 1px;
        }}
        .section {{
            margin-bottom: 40px;
        }}
        .section h2 {{
            color: #2c3e50;
            border-bottom: 2px solid #ecf0f1;
            padding-bottom: 10px;
            margin-bottom: 20px;
            display: flex;
            align-items: center;
        }}
        .section h2::before {{
            content: "📊";
            margin-right: 10px;
            font-size: 1.2em;
        }}
        .reports-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fill, minmax(300px, 1fr));
            gap: 20px;
        }}
        .report-card {{
            background: white;
            border: 1px solid #ecf0f1;
            border-radius: 8px;
            padding: 20px;
            transition: all 0.3s ease;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }}
        .report-card:hover {{
            transform: translateY(-5px);
            box-shadow: 0 5px 20px rgba(0,0,0,0.15);
            border-color: #3498db;
        }}
        .report-title {{
            font-weight: 600;
            color: #2c3e50;
            margin-bottom: 8px;
            font-size: 1.1em;
        }}
        .report-description {{
            color: #7f8c8d;
            font-size: 0.9em;
            margin-bottom: 15px;
            line-height: 1.4;
        }}
        .report-link {{
            display: inline-block;
            background: #3498db;
            color: white;
            padding: 8px 16px;
            border-radius: 5px;
            text-decoration: none;
            font-size: 0.9em;
            transition: background 0.3s ease;
        }}
        .report-link:hover {{
            background: #2980b9;
        }}
        .footer {{
            margin-top: 40px;
            padding-top: 20px;
            border-top: 1px solid #ecf0f1;
            text-align: center;
            color: #7f8c8d;
            font-size: 0.9em;
        }}
        .enhanced {{
            border-left: 4px solid #e74c3c;
        }}
        .combined {{
            border-left: 4px solid #f39c12;
        }}
        .individual {{
            border-left: 4px solid #27ae60;
        }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>🚀 Power-Law GoF Analysis</h1>
            <p>Navigation Index for All Generated Reports</p>
        </div>
        
        <div class="content">
            <div class="stats">
                <div class="stat-card">
                    <div class="stat-number">{len(analyses)}</div>
                    <div class="stat-label">Files Analyzed</div>
                </div>
                <div class="stat-card">
                    <div class="stat-number">{sum(len(content["blocks"]) for content in analyses.values())}</div>
                    <div class="stat-label">Total Blocks</div>
                </div>
                <div class="stat-card">
                    <div class="stat-number">{report_count}</div>
                    <div class="stat-label">Reports Generated</div>
                </div>
                <div class="stat-card">
                    <div class="stat-number">{len(html_files)}</div>
                    <div class="stat-label">HTML Files</div>
                </div>
            </div>"""
    
    # Individual File Reports Section
    if per_file_reports:
        index_html += f"""
            <div class="section">
                <h2>Individual File Reports</h2>
                <div class="reports-grid">"""
        for file in sorted(per_file_reports):
            # Extract clean name from filename
            clean_name = file.replace("GoF_FAST_", "").replace(".html", "").replace("_", " ")
            index_html += f"""
                    <div class="report-card individual">
                        <div class="report-title">{clean_name}</div>
                        <div class="report-description">Detailed power-law goodness-of-fit analysis for individual file with block-level statistics, aggregated summaries, and integrated visualizations.</div>
                        <a href="{file}" class="report-link">View Report →</a>
                    </div>"""
        index_html += """
                </div>
            </div>"""
    
    # Combined Reports Section
    if combined_reports:
        index_html += f"""
            <div class="section">
                <h2>Combined Analysis Reports</h2>
                <div class="reports-grid">"""
        for file in sorted(combined_reports):
            if "ALL_FILES" in file:
                title = "All Files Combined"
                desc = "Comprehensive overview of all analyzed files with comparative statistics and cross-file analysis."
            elif "COMPARISONS" in file:
                title = "Statistical Comparisons"
                desc = "Statistical comparisons between files, tissues, and region pairs with p-values and significance tests."
            else:
                title = file.replace("GoF_FAST_", "").replace(".html", "").replace("_", " ")
                desc = "Combined analysis report"
            
            index_html += f"""
                    <div class="report-card combined">
                        <div class="report-title">{title}</div>
                        <div class="report-description">{desc}</div>
                        <a href="{file}" class="report-link">View Report →</a>
                    </div>"""
        index_html += """
                </div>
            </div>"""
    
    # Enhanced Visualizations Section
    if enhanced_reports:
        index_html += f"""
            <div class="section">
                <h2>Enhanced Interactive Visualizations</h2>
                <div class="reports-grid">"""
        for file in sorted(enhanced_reports):
            if "COMBINED_TISSUES" in file:
                title = "🎨 Combined Tissues Visualization"
                desc = "Interactive multi-tissue comparison with integrated power-law fitting curves, distribution plots, and statistical overlays."
            elif "TS_vs_CT" in file:
                title = "⚡ TS vs CT Comparison"
                desc = "Direct comparison between TS (tissue-specific) and CT (cross-tissue) analyses with side-by-side visualizations."
            else:
                title = file.replace("GoF_FAST_", "").replace(".html", "").replace("_", " ")
                desc = "Enhanced interactive visualization"
            
            index_html += f"""
                    <div class="report-card enhanced">
                        <div class="report-title">{title}</div>
                        <div class="report-description">{desc}</div>
                        <a href="{file}" class="report-link">View Interactive Plot →</a>
                    </div>"""
        index_html += """
                </div>
            </div>"""
    
    # Statistical Test Reports Section
    if test_reports:
        index_html += f"""
            <div class="section">
                <h2>🧪 Individual Statistical Tests</h2>
                <div class="reports-grid">"""
        for file in sorted(test_reports):
            if "MannWhitney" in file:
                title = "📊 Mann-Whitney U Test"
                desc = "Pairwise comparison between tissues with effect sizes, distributions, and detailed statistical analysis."
                parts = file.replace("TEST_MannWhitney_", "").replace(".html", "").split("_vs_")
                if len(parts) == 2:
                    title += f": {parts[0]} vs {parts[1]}"
            elif "KruskalWallis" in file:
                title = "🔬 Kruskal-Wallis Test"
                desc = "Multi-group comparison across all tissues with rank analysis and post-hoc insights."
            elif "Correlations" in file:
                title = "📈 Correlation Analysis"
                desc = "Comprehensive correlation analysis between alpha, KS-D, R², and p-values with regression summaries."
                tissue_name = file.replace("TEST_Correlations_", "").replace(".html", "")
                title += f": {tissue_name}"
            else:
                title = file.replace("TEST_", "").replace(".html", "").replace("_", " ").title()
                desc = "Individual statistical test analysis with detailed visualizations."
            
            index_html += f"""
                    <div class="report-card enhanced">
                        <div class="report-title">{title}</div>
                        <div class="report-description">{desc}</div>
                        <a href="{file}" class="report-link">View Test Results →</a>
                    </div>"""
        index_html += """
                </div>
            </div>"""
    
    # KS Visualization Reports Section
    if ks_reports or individual_ks:
        index_html += f"""
            <div class="section">
                <h2>🎯 KS Statistic Visualizations</h2>
                <div class="reports-grid">"""
        
        # Main KS visualization report
        for file in sorted(ks_reports):
            index_html += f"""
                    <div class="report-card enhanced">
                        <div class="report-title">🎯 KS Test Overview</div>
                        <div class="report-description">Comprehensive KS statistic visualization across all tissues with empirical vs theoretical CDF comparisons.</div>
                        <a href="{file}" class="report-link">View KS Analysis →</a>
                    </div>"""
        
        # Individual KS reports (first 6 for display)
        for file in sorted(individual_ks)[:6]:
            tissue_pair = file.replace("GoF_FAST_KS_", "").replace(".html", "")
            parts = tissue_pair.split("_")
            if len(parts) >= 4:
                tissue = parts[0]
                region_pair = "_".join(parts[1:-1])  # Everything except last part (index)
                index_html += f"""
                    <div class="report-card individual">
                        <div class="report-title">📊 KS Detail: {tissue}</div>
                        <div class="report-description">Individual KS analysis for {region_pair} with detailed CDF comparison and power-law fit visualization.</div>
                        <a href="{file}" class="report-link">View Details →</a>
                    </div>"""
        
        # Show count if there are more individual reports
        if len(individual_ks) > 6:
            index_html += f"""
                    <div class="report-card">
                        <div class="report-title">➕ Additional KS Reports</div>
                        <div class="report-description">{len(individual_ks) - 6} more individual KS visualizations available in the output directory.</div>
                    </div>"""
        
        index_html += """
                </div>
            </div>"""
    
    # Barabási Network Science Analysis Section
    if barabasi_reports:
        index_html += f"""
            <div class="section">
                <h2>🔬 Barabási Network Science Analysis</h2>
                <p style="margin-bottom: 20px; color: #555; font-style: italic;">
                    Comprehensive scale-free analysis following Barabási's Network Science Chapter 4 recommendations.
                    Implements three-step goodness-of-fit testing with MLE parameter estimation, 
                    KS/Anderson-Darling tests, and comparison against alternative heavy-tailed distributions.
                </p>
                <div class="reports-grid">"""
        
        for file in sorted(barabasi_reports):
            filename = os.path.basename(file)
            if "comparative" in filename.lower():
                title = "🧪 Comprehensive Comparative Analysis"
                desc = "Multi-tissue scale-free property comparison with Barabási criteria evaluation, parameter distributions, and evidence assessment."
            elif "summary" in filename.lower():
                title = "📊 Statistical Summary"
                desc = "Comprehensive statistical summary of scale-free analysis results across all tissues with evidence distribution."
            else:
                # Extract tissue name from filename
                tissue_name = filename.replace("barabasi_", "").replace(".html", "").replace("_", " ").title()
                title = f"🔬 {tissue_name} Scale-Free Analysis"
                desc = f"Detailed goodness-of-fit analysis for {tissue_name} with power-law fitting, alternative distribution comparison, and Barabási criteria evaluation."
            
            index_html += f"""
                    <div class="report-card enhanced">
                        <div class="report-title">{title}</div>
                        <div class="report-description">{desc}</div>
                        <a href="{file}" class="report-link">View Analysis →</a>
                    </div>"""
        
        index_html += """
                </div>
            </div>"""

    # Other Reports Section
    if other_reports:
        index_html += f"""
            <div class="section">
                <h2>Additional Reports</h2>
                <div class="reports-grid">"""
        for file in sorted(other_reports):
            clean_name = file.replace(".html", "").replace("_", " ").title()
            index_html += f"""
                    <div class="report-card">
                        <div class="report-title">{clean_name}</div>
                        <div class="report-description">Additional analysis report or visualization.</div>
                        <a href="{file}" class="report-link">View Report →</a>
                    </div>"""
        index_html += """
                </div>
            </div>"""
    
    index_html += f"""
            <div class="footer">
                <p>Generated on {current_time} | Total processing completed successfully</p>
            </div>
        </div>
    </div>
</body>
</html>"""
    
    # Write index file
    index_path = os.path.join(out_dir, "index.html")
    with open(index_path, 'w', encoding='utf-8') as f:
        f.write(index_html)
    
    print(f"      📄 Navigation index saved: {os.path.basename(index_path)}")

# -------------------- Enhanced Comparison Functions --------------------

def _create_barabasi_goodness_of_fit_report(
    combined_results: Dict[str, List[Tuple]],
    output_dir: str,
    bin_scheme: str = "quantile"
) -> str:
    """
    Create comprehensive Barabási Network Science-inspired goodness-of-fit report.
    
    This implements the three-step process recommended in Chapter 4:
    1. Estimate power law parameters using MLE
    2. Calculate goodness-of-fit using KS test and alternatives
    3. Compare against alternative heavy-tailed distributions
    """
    if not BARABASI_AVAILABLE:
        print("⚠ Barabási analysis not available - skipping comprehensive GoF report")
        return ""
    
    print("\n🔬 Creating Barabási-inspired comprehensive goodness-of-fit analysis...")
    
    # Create output directory for Barabási analysis
    barabasi_dir = Path(output_dir) / "barabasi_analysis"
    barabasi_dir.mkdir(exist_ok=True)
    
    # Collect all degree sequences
    all_analyses = []
    tissue_data = {}
    
    for tissue_type, results_list in combined_results.items():
        print(f"  📊 Processing {tissue_type} data ({len(results_list)} files)...")
        
        degree_sequences = []
        for file_result in results_list:
            if len(file_result) >= 4:
                # Extract connectivity values from each block
                block_data = file_result[3] if len(file_result) > 3 else []
                for block_result in block_data:
                    if hasattr(block_result, 'combined_values') and len(block_result.combined_values) > 0:
                        degree_sequences.extend(block_result.combined_values)
        
        if degree_sequences:
            # Convert to numpy array and filter
            degrees = np.array(degree_sequences)
            degrees = degrees[degrees > 0]  # Remove zeros
            degrees = degrees[np.isfinite(degrees)]  # Remove inf/nan
            
            if len(degrees) > 100:  # Minimum for reliable analysis
                print(f"    🎯 Analyzing {len(degrees):,} degree values for {tissue_type}")
                
                # Perform comprehensive Barabási analysis
                analysis = analyze_scale_free_properties(degrees, xmin_method='clauset')
                analysis.tissue_type = tissue_type
                analysis.bin_scheme = bin_scheme
                
                all_analyses.append(analysis)
                tissue_data[tissue_type] = degrees
                
                # Save individual tissue analysis
                tissue_output = str(barabasi_dir / f"barabasi_{tissue_type}")
                html_path, summary_path = save_scale_free_analysis_report(analysis, tissue_output)
                
                print(f"    ✅ {tissue_type}: Evidence={analysis.scale_free_evidence}, α={analysis.alpha:.3f}, p={_format_pvalue(analysis.power_law_fit.ks_pvalue)}")
            else:
                print(f"    ⚠ Insufficient data for {tissue_type}: {len(degrees)} values")
    
    if not all_analyses:
        print("❌ No sufficient data for Barabási analysis")
        return ""
    
    # Create comparative analysis
    print("\n🔬 Creating comparative Barabási analysis...")
    comparative_fig = _create_barabasi_comparative_analysis(all_analyses, tissue_data)
    
    # Save comparative report
    comparative_path = str(barabasi_dir / "barabasi_comparative_analysis.html")
    comparative_fig.write_html(comparative_path)
    
    # Create summary statistics
    summary_stats = _generate_barabasi_summary_statistics(all_analyses)
    summary_path = str(barabasi_dir / "barabasi_summary_statistics.txt")
    
    with open(summary_path, 'w') as f:
        f.write("COMPREHENSIVE BARABÁSI NETWORK SCIENCE ANALYSIS SUMMARY\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"Analysis performed following Barabási's Network Science Chapter 4 recommendations\n")
        f.write(f"Bin scheme: {bin_scheme}\n")
        f.write(f"Total tissues analyzed: {len(all_analyses)}\n\n")
        
        f.write("SCALE-FREE EVIDENCE SUMMARY:\n")
        f.write("-" * 30 + "\n")
        evidence_counts = {}
        for analysis in all_analyses:
            evidence = analysis.scale_free_evidence
            evidence_counts[evidence] = evidence_counts.get(evidence, 0) + 1
            
        for evidence, count in evidence_counts.items():
            f.write(f"{evidence}: {count} tissues\n")
        
        f.write("\nPER-TISSUE ANALYSIS RESULTS:\n")
        f.write("-" * 35 + "\n")
        f.write("Tissue\t\tα\tKS p-value\tEvidence\tBest Dist\n")
        f.write("-" * 60 + "\n")
        
        for analysis in all_analyses:
            f.write(f"{analysis.tissue_type[:12]:<12}\t{analysis.alpha:.3f}\t{_format_pvalue(analysis.power_law_fit.ks_pvalue)}\t\t{analysis.scale_free_evidence:<8}\t{analysis.best_distribution}\n")
        
        f.write("\nBARABÁSI CRITERIA SATISFACTION:\n")
        f.write("-" * 35 + "\n")
        criteria_stats = {}
        for analysis in all_analyses:
            for criterion, satisfied in analysis.barabasi_criteria.items():
                if criterion not in criteria_stats:
                    criteria_stats[criterion] = {'satisfied': 0, 'total': 0}
                criteria_stats[criterion]['total'] += 1
                if satisfied:
                    criteria_stats[criterion]['satisfied'] += 1
        
        for criterion, stats in criteria_stats.items():
            satisfaction_rate = stats['satisfied'] / stats['total'] * 100
            f.write(f"{criterion}: {stats['satisfied']}/{stats['total']} ({satisfaction_rate:.1f}%)\n")
        
        f.write(f"\n{summary_stats}")
    
    print(f"✅ Barabási comprehensive analysis complete!")
    print(f"📁 Results saved to: {barabasi_dir}")
    print(f"📊 Comparative analysis: {comparative_path}")
    print(f"📋 Summary statistics: {summary_path}")
    
    return str(barabasi_dir)


def _create_barabasi_comparative_analysis(analyses: List[ScaleFreeAnalysis], tissue_data: Dict) -> go.Figure:
    """Create comparative Barabási analysis across tissues."""
    
    fig = make_subplots(
        rows=3, cols=3,
        subplot_titles=(
            "Power Law Exponents (α)", "KS Test P-values", "Scale-Free Evidence",
            "AIC Comparison", "Alternative Distribution Performance", "Barabási Criteria Heatmap",
            "Parameter Distributions", "Goodness-of-Fit Comparison", "Network Science Assessment"
        ),
        specs=[[{"secondary_y": False}, {"secondary_y": False}, {"secondary_y": False}],
               [{"secondary_y": False}, {"secondary_y": False}, {"secondary_y": False}],
               [{"secondary_y": False}, {"secondary_y": False}, {"secondary_y": False}]]
    )
    
    tissues = [a.tissue_type for a in analyses]
    alphas = [a.alpha for a in analyses]
    ks_pvals = [a.power_law_fit.ks_pvalue for a in analyses]
    evidence = [a.scale_free_evidence for a in analyses]
    
    # 1. Power law exponents
    colors = ['red' if 2.0 <= alpha <= 3.5 else 'blue' for alpha in alphas]
    fig.add_trace(go.Bar(x=tissues, y=alphas, name='α values', marker=dict(color=colors)),
                  row=1, col=1)
    fig.add_hline(y=2.0, line_dash="dash", line_color="gray", row=1, col=1)
    fig.add_hline(y=3.5, line_dash="dash", line_color="gray", row=1, col=1)
    
    # 2. KS test p-values
    colors = ['green' if p > 0.1 else 'orange' if p > 0.05 else 'red' for p in ks_pvals]
    fig.add_trace(go.Bar(x=tissues, y=ks_pvals, name='KS p-values', marker=dict(color=colors)),
                  row=1, col=2)
    fig.add_hline(y=0.1, line_dash="dash", line_color="gray", row=1, col=2)
    fig.add_hline(y=0.05, line_dash="dash", line_color="gray", row=1, col=2)
    
    # 3. Scale-free evidence
    evidence_map = {'Strong': 4, 'Moderate': 3, 'Weak': 2, 'None': 1}
    evidence_values = [evidence_map.get(e, 0) for e in evidence]
    colors = ['darkgreen' if e == 'Strong' else 'green' if e == 'Moderate' 
              else 'orange' if e == 'Weak' else 'red' for e in evidence]
    
    fig.add_trace(go.Bar(x=tissues, y=evidence_values, name='Evidence Strength',
                        marker=dict(color=colors), 
                        customdata=evidence,
                        hovertemplate='<b>%{x}</b><br>Evidence: %{customdata}<extra></extra>'),
                  row=1, col=3)
    
    # 4. AIC comparison
    power_law_aics = [a.power_law_fit.aic for a in analyses if not np.isnan(a.power_law_fit.aic)]
    best_alt_aics = []
    
    for a in analyses:
        if a.alternative_fits:
            alt_aics = [alt.aic for alt in a.alternative_fits if not np.isnan(alt.aic)]
            best_alt_aics.append(min(alt_aics) if alt_aics else np.nan)
        else:
            best_alt_aics.append(np.nan)
    
    valid_indices = [i for i, (pl, alt) in enumerate(zip(power_law_aics, best_alt_aics)) 
                     if not np.isnan(pl) and not np.isnan(alt)]
    
    if valid_indices:
        fig.add_trace(go.Scatter(x=[power_law_aics[i] for i in valid_indices],
                                y=[best_alt_aics[i] for i in valid_indices],
                                mode='markers+text',
                                text=[tissues[i] for i in valid_indices],
                                textposition="top center",
                                name='AIC: Power Law vs Best Alternative',
                                marker=dict(size=10)), row=2, col=1)
        
        # Add diagonal line
        min_aic = min(min(power_law_aics), min([x for x in best_alt_aics if not np.isnan(x)]))
        max_aic = max(max(power_law_aics), max([x for x in best_alt_aics if not np.isnan(x)]))
        fig.add_trace(go.Scatter(x=[min_aic, max_aic], y=[min_aic, max_aic],
                                mode='lines', line=dict(dash='dash', color='gray'),
                                name='Equal AIC', showlegend=False), row=2, col=1)
    
    # 5. Alternative distribution performance
    alt_names = set()
    for a in analyses:
        for alt in a.alternative_fits:
            alt_names.add(alt.name)
    
    alt_names = list(alt_names)
    alt_performance = []
    
    for alt_name in alt_names:
        performance = []
        for a in analyses:
            alt_fit = next((alt for alt in a.alternative_fits if alt.name == alt_name), None)
            if alt_fit and not np.isnan(alt_fit.ks_pvalue):
                performance.append(alt_fit.ks_pvalue)
        
        if performance:
            fig.add_trace(go.Box(y=performance, name=alt_name), row=2, col=2)
    
    # 6. Barabási criteria heatmap
    criteria_names = list(analyses[0].barabasi_criteria.keys()) if analyses else []
    criteria_matrix = []
    
    for tissue in tissues:
        analysis = next(a for a in analyses if a.tissue_type == tissue)
        criteria_row = [1 if analysis.barabasi_criteria.get(criterion, False) else 0 
                       for criterion in criteria_names]
        criteria_matrix.append(criteria_row)
    
    if criteria_matrix:
        fig.add_trace(go.Heatmap(z=criteria_matrix, x=criteria_names, y=tissues,
                                colorscale='RdYlGn', showscale=False,
                                hovertemplate='<b>%{y}</b><br>%{x}: %{z}<extra></extra>'),
                      row=2, col=3)
    
    # 7-9: Additional analyses (simplified for now)
    for row, col, title in [(3, 1, "Parameter Distributions"), (3, 2, "Goodness-of-Fit Comparison"), 
                           (3, 3, "Network Science Assessment")]:
        fig.add_annotation(text=f"{title}<br>(Enhanced analysis)", 
                          xref="x domain", yref="y domain", x=0.5, y=0.5, 
                          font=dict(size=14), row=row, col=col)
    
    # Update layout
    fig.update_layout(
        title="Comprehensive Barabási Network Science Analysis Comparison",
        width=1800, height=1400,
        showlegend=True
    )
    
    # Update y-axis labels
    fig.update_yaxes(title_text="Power Law Exponent (α)", row=1, col=1)
    fig.update_yaxes(title_text="KS Test P-value", row=1, col=2)
    fig.update_yaxes(title_text="Evidence Level", row=1, col=3)
    fig.update_xaxes(title_text="Power Law AIC", row=2, col=1)
    fig.update_yaxes(title_text="Best Alternative AIC", row=2, col=1)
    fig.update_yaxes(title_text="P-value", row=2, col=2)
    
    return fig


def _generate_barabasi_summary_statistics(analyses: List[ScaleFreeAnalysis]) -> str:
    """Generate comprehensive summary statistics for Barabási analysis."""
    
    if not analyses:
        return "No analyses available for summary statistics."
    
    # Extract key metrics
    alphas = [a.alpha for a in analyses if not np.isnan(a.alpha)]
    ks_pvals = [a.power_law_fit.ks_pvalue for a in analyses if not np.isnan(a.power_law_fit.ks_pvalue)]
    goodness_scores = [a.power_law_fit.goodness_score for a in analyses if not np.isnan(a.power_law_fit.goodness_score)]
    
    summary = "\nSTATISTICAL SUMMARY:\n"
    summary += "-" * 25 + "\n"
    
    if alphas:
        summary += f"Power Law Exponent (α):\n"
        summary += f"  Mean: {np.mean(alphas):.3f} ± {np.std(alphas):.3f}\n"
        summary += f"  Median: {np.median(alphas):.3f}\n"
        summary += f"  Range: [{np.min(alphas):.3f}, {np.max(alphas):.3f}]\n"
        summary += f"  Scale-free range (2.0-3.5): {np.sum((2.0 <= np.array(alphas)) & (np.array(alphas) <= 3.5))}/{len(alphas)} ({100*np.sum((2.0 <= np.array(alphas)) & (np.array(alphas) <= 3.5))/len(alphas):.1f}%)\n\n"
    
    if ks_pvals:
        summary += f"KS Test P-values:\n"
        summary += f"  Mean: {_format_pvalue(np.mean(ks_pvals))} ± {_format_pvalue(np.std(ks_pvals))}\n"
        summary += f"  Median: {_format_pvalue(np.median(ks_pvals))}\n"
        summary += f"  P > 0.1 (good fit): {np.sum(np.array(ks_pvals) > 0.1)}/{len(ks_pvals)} ({100*np.sum(np.array(ks_pvals) > 0.1)/len(ks_pvals):.1f}%)\n"
        summary += f"  P > 0.05 (acceptable): {np.sum(np.array(ks_pvals) > 0.05)}/{len(ks_pvals)} ({100*np.sum(np.array(ks_pvals) > 0.05)/len(ks_pvals):.1f}%)\n\n"
    
    if goodness_scores:
        summary += f"Combined Goodness Scores:\n"
        summary += f"  Mean: {np.mean(goodness_scores):.3f} ± {np.std(goodness_scores):.3f}\n"
        summary += f"  Median: {np.median(goodness_scores):.3f}\n"
        summary += f"  Range: [{np.min(goodness_scores):.3f}, {np.max(goodness_scores):.3f}]\n\n"
    
    # Evidence distribution
    evidence_counts = {}
    for analysis in analyses:
        evidence = analysis.scale_free_evidence
        evidence_counts[evidence] = evidence_counts.get(evidence, 0) + 1
    
    summary += f"Scale-Free Evidence Distribution:\n"
    total = len(analyses)
    for evidence in ['Strong', 'Moderate', 'Weak', 'None']:
        count = evidence_counts.get(evidence, 0)
        percentage = 100 * count / total if total > 0 else 0
        summary += f"  {evidence}: {count}/{total} ({percentage:.1f}%)\n"
    
    return summary


def _create_tissue_alpha_comparison(analyses: Dict[str, Dict]) -> go.Figure:
    """Create focused comparison of alpha distributions across tissues."""
    
    tissues = list(analyses.keys())[:6]  # Limit for clarity
    if len(tissues) < 2:
        fig = go.Figure()
        fig.add_annotation(text="Need at least 2 tissues for comparison", 
                          xref="paper", yref="paper", x=0.5, y=0.5, showarrow=False)
        return fig

    fig = make_subplots(
        rows=2, cols=2,
        subplot_titles=(
            "📊 Alpha Distribution by Tissue<br><sub>Violin plots showing power-law exponent distributions</sub>", 
            "📈 Alpha vs R² Scatter<br><sub>Relationship between exponent and goodness-of-fit</sub>", 
            "🎯 Statistical Significance<br><sub>P-value distributions for power-law fits</sub>", 
            "📋 Summary Statistics<br><sub>Mean, median, and variability comparison</sub>"
        ),
        vertical_spacing=0.15, horizontal_spacing=0.12
    )
    
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b']
    
    # Collect data
    tissue_data = {}
    for i, tissue in enumerate(tissues):
        blocks = analyses[tissue]["blocks"]
        alphas = np.array([b.alpha for b in blocks if np.isfinite(b.alpha)])
        r2_vals = np.array([b.R2_k_loglog for b in blocks if np.isfinite(b.R2_k_loglog)])
        p_vals = np.array([b.ks_pvalue for b in blocks if np.isfinite(b.ks_pvalue)])
        
        tissue_data[tissue] = {
            'alphas': alphas, 'r2_vals': r2_vals, 'p_vals': p_vals,
            'color': colors[i % len(colors)]
        }
    
    # Plot 1: Alpha distributions
    for tissue, data in tissue_data.items():
        if len(data['alphas']) > 0:
            fig.add_trace(go.Violin(
                y=data['alphas'], name=tissue, 
                line_color=data['color'], fillcolor=data['color'], opacity=0.6,
                box_visible=True, meanline_visible=True,
                hovertemplate=f"<b>{tissue}</b><br>Alpha: %{{y:.3f}}<extra></extra>"
            ), row=1, col=1)
    
    # Plot 2: Alpha vs R² scatter
    for tissue, data in tissue_data.items():
        if len(data['alphas']) > 0 and len(data['r2_vals']) > 0:
            # Match arrays by length
            min_len = min(len(data['alphas']), len(data['r2_vals']))
            fig.add_trace(go.Scatter(
                x=data['alphas'][:min_len], y=data['r2_vals'][:min_len],
                mode='markers', name=tissue, 
                marker=dict(color=data['color'], size=6, opacity=0.7),
                hovertemplate=f"<b>{tissue}</b><br>Alpha: %{{x:.3f}}<br>R²: %{{y:.3f}}<extra></extra>"
            ), row=1, col=2)
    
    # Plot 3: P-value distributions
    for tissue, data in tissue_data.items():
        if len(data['p_vals']) > 0:
            fig.add_trace(go.Histogram(
                x=data['p_vals'], name=tissue, 
                marker_color=data['color'], opacity=0.7, nbinsx=20,
                hovertemplate=f"<b>{tissue}</b><br>P-value: %{{x:.2e}}<extra></extra>"
            ), row=2, col=1)
    
    # Plot 4: Summary statistics
    summary_data = []
    for tissue, data in tissue_data.items():
        if len(data['alphas']) > 0:
            summary_data.append({
                'tissue': tissue,
                'mean': np.mean(data['alphas']),
                'median': np.median(data['alphas']),
                'std': np.std(data['alphas']),
                'color': data['color']
            })
    
    if summary_data:
        tissues_list = [d['tissue'] for d in summary_data]
        means = [d['mean'] for d in summary_data]
        medians = [d['median'] for d in summary_data]
        stds = [d['std'] for d in summary_data]
        colors_list = [d['color'] for d in summary_data]
        
        fig.add_trace(go.Bar(
            x=tissues_list, y=means, name='Mean', 
            marker_color=colors_list, opacity=0.8,
            hovertemplate="<b>%{x}</b><br>Mean: %{y:.3f}<extra></extra>"
        ), row=2, col=2)
    
    fig.update_layout(
        title="📊 Tissue Alpha Comparison - Focused View",
        height=800, showlegend=True,
        plot_bgcolor='white', paper_bgcolor='#fafbfc'
    )
    
    return fig


def _create_tissue_statistical_tests(analyses: Dict[str, Dict]) -> go.Figure:
    """Create focused statistical comparison tests between tissues."""
    
    tissues = list(analyses.keys())[:6]
    if len(tissues) < 2:
        fig = go.Figure()
        fig.add_annotation(text="Need at least 2 tissues for statistical tests", 
                          xref="paper", yref="paper", x=0.5, y=0.5, showarrow=False)
        return fig

    fig = make_subplots(
        rows=2, cols=2,
        subplot_titles=(
            "🧪 Mann-Whitney U Tests<br><sub>Pairwise comparisons between tissues (lower p = more different)</sub>", 
            "📊 Effect Sizes (Cohen's d)<br><sub>Magnitude of differences between tissues</sub>", 
            "🌡️ Distance Matrix<br><sub>Earth Mover's Distance between tissue distributions</sub>", 
            "📈 Power Analysis<br><sub>Statistical power and sample sizes</sub>"
        ),
        specs=[[{"secondary_y": False}, {"secondary_y": False}],
               [{"type": "heatmap"}, {"secondary_y": False}]],
        vertical_spacing=0.15, horizontal_spacing=0.12
    )
    
    # Collect alpha data
    tissue_alphas = {}
    for tissue in tissues:
        blocks = analyses[tissue]["blocks"]
        alphas = np.array([b.alpha for b in blocks if np.isfinite(b.alpha)])
        tissue_alphas[tissue] = alphas
    
    # Statistical tests
    comparisons = []
    p_values = []
    effect_sizes = []
    
    for i in range(len(tissues)):
        for j in range(i+1, len(tissues)):
            if len(comparisons) >= 12:  # Limit for readability
                break
            
            tissue1, tissue2 = tissues[i], tissues[j]
            alphas1 = tissue_alphas[tissue1] 
            alphas2 = tissue_alphas[tissue2]
            
            if len(alphas1) > 0 and len(alphas2) > 0:
                try:
                    from scipy.stats import mannwhitneyu
                    _, p_val = mannwhitneyu(alphas1, alphas2, alternative='two-sided')
                    
                    # Cohen's d
                    pooled_std = np.sqrt(((len(alphas1) - 1) * np.var(alphas1, ddof=1) + 
                                         (len(alphas2) - 1) * np.var(alphas2, ddof=1)) / 
                                        (len(alphas1) + len(alphas2) - 2))
                    cohens_d = (np.mean(alphas1) - np.mean(alphas2)) / pooled_std if pooled_std > 0 else 0
                    
                    comparisons.append(f"{tissue1} vs {tissue2}")
                    p_values.append(p_val)
                    effect_sizes.append(abs(cohens_d))
                except:
                    continue
    
    # Plot 1: Mann-Whitney U tests
    if comparisons:
        fig.add_trace(go.Bar(
            x=comparisons, y=p_values, name="Mann-Whitney p-values",
            marker_color='lightblue', opacity=0.8,
            hovertemplate="<b>%{x}</b><br>p-value: %{y:.2e}<br>Significant: %{text}<extra></extra>",
            text=[f"{'Yes' if p < 0.05 else 'No'}" for p in p_values]
        ), row=1, col=1)
        
        fig.add_hline(y=0.05, line_dash="dash", line_color="red", 
                     annotation_text="p = 0.05", row=1, col=1)
    
    # Plot 2: Effect sizes
    if effect_sizes:
        fig.add_trace(go.Bar(
            x=comparisons, y=effect_sizes, name="Cohen's d",
            marker_color='mediumpurple', opacity=0.8,
            hovertemplate="<b>%{x}</b><br>Effect size: %{y:.3f}<extra></extra>"
        ), row=1, col=2)
        
        # Effect size interpretation lines
        fig.add_hline(y=0.2, line_dash="dash", line_color="green", 
                     annotation_text="Small (0.2)", row=1, col=2)
        fig.add_hline(y=0.5, line_dash="dash", line_color="orange", 
                     annotation_text="Medium (0.5)", row=1, col=2)
        fig.add_hline(y=0.8, line_dash="dash", line_color="red", 
                     annotation_text="Large (0.8)", row=1, col=2)
    
    # Plot 3: Distance matrix
    if len(tissues) >= 2:
        distance_matrix = np.zeros((len(tissues), len(tissues)))
        for i, tissue1 in enumerate(tissues):
            for j, tissue2 in enumerate(tissues):
                if i != j:
                    alphas1 = tissue_alphas[tissue1]
                    alphas2 = tissue_alphas[tissue2]
                    if len(alphas1) > 0 and len(alphas2) > 0:
                        try:
                            from scipy.stats import wasserstein_distance
                            dist = wasserstein_distance(alphas1, alphas2)
                            distance_matrix[i, j] = dist
                        except:
                            distance_matrix[i, j] = 0
        
        fig.add_trace(go.Heatmap(
            z=distance_matrix, x=tissues, y=tissues,
            colorscale='Viridis', name="Distance Matrix",
            hovertemplate="<b>%{x} vs %{y}</b><br>Distance: %{z:.4f}<extra></extra>",
            colorbar=dict(title="EMD Distance")
        ), row=2, col=1)
    
    # Plot 4: Sample sizes and power
    sample_sizes = [len(tissue_alphas[tissue]) for tissue in tissues]
    fig.add_trace(go.Bar(
        x=tissues, y=sample_sizes, name="Sample Sizes",
        marker_color='lightgreen', opacity=0.8,
        hovertemplate="<b>%{x}</b><br>Sample size: %{y}<extra></extra>"
    ), row=2, col=2)
    
    fig.update_layout(
        title="🧪 Statistical Tests Between Tissues",
        height=800, showlegend=True,
        plot_bgcolor='white', paper_bgcolor='#fafbfc'
    )
    
    return fig


def _create_enhanced_tissue_comparison(analyses: Dict[str, Dict]) -> go.Figure:
    """Create clean, organized tissue comparison with clear sections."""
    
    tissues = list(analyses.keys())[:6]  # Limit for clarity
    
    if len(tissues) < 2:
        fig = go.Figure()
        fig.add_annotation(text="Need at least 2 tissues for comparison", 
                          xref="paper", yref="paper", x=0.5, y=0.5,
                          showarrow=False, font=dict(size=16))
        fig.update_layout(title="Enhanced Tissue Comparison", height=400)
        return fig

    # Create main overview with 2x2 layout - much cleaner!
    fig = make_subplots(
        rows=2, cols=2,
        subplot_titles=(
            "📊 Power-Law Exponents (α)<br><sub>Lower values = more scale-free behavior</sub>", 
            "🎯 Goodness-of-Fit (R²)<br><sub>Higher values = better power-law fit</sub>", 
            "📈 Statistical Significance<br><sub>P-values > 0.1 indicate good fits</sub>", 
            "🔬 Tissue Comparisons<br><sub>Statistical differences between tissues</sub>"
        ),
        vertical_spacing=0.15, horizontal_spacing=0.12
    )
    
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b']
    
    # Collect data for each tissue
    tissue_data = {}
    for i, tissue in enumerate(tissues):
        blocks = analyses[tissue]["blocks"]
        alphas = np.array([b.alpha for b in blocks if np.isfinite(b.alpha)])
        r2_vals = np.array([b.R2_k_loglog for b in blocks if np.isfinite(b.R2_k_loglog)])
        p_vals = np.array([b.ks_pvalue for b in blocks if np.isfinite(b.ks_pvalue)])
        
        tissue_data[tissue] = {
            'alphas': alphas, 'r2_vals': r2_vals, 'p_vals': p_vals,
            'color': colors[i % len(colors)]
        }
    
    # Plot 1: Alpha distributions - Clean violin plots
    for tissue, data in tissue_data.items():
        if len(data['alphas']) > 0:
            fig.add_trace(go.Violin(
                y=data['alphas'], name=tissue, 
                line_color=data['color'], fillcolor=data['color'], opacity=0.6,
                box_visible=True, meanline_visible=True,
                hovertemplate=f"<b>{tissue}</b><br>" +
                             f"Alpha: %{{y:.3f}}<br>" +
                             f"Mean: {np.mean(data['alphas']):.3f}<br>" +
                             f"N: {len(data['alphas'])}<extra></extra>"
            ), row=1, col=1)
    
    # Plot 2: R² distributions - Clean box plots  
    for tissue, data in tissue_data.items():
        if len(data['r2_vals']) > 0:
            fig.add_trace(go.Box(
                y=data['r2_vals'], name=tissue, 
                marker_color=data['color'], opacity=0.7, boxmean=True,
                hovertemplate=f"<b>{tissue}</b><br>" +
                             f"R²: %{{y:.3f}}<br>" +
                             f"Median: {np.median(data['r2_vals']):.3f}<extra></extra>"
            ), row=1, col=2)
    
    # Plot 3: P-value significance - Clean scatter
    for i, (tissue, data) in enumerate(tissue_data.items()):
        if len(data['p_vals']) > 0:
            # Add jitter for better visualization
            x_pos = np.full(len(data['p_vals']), i) + np.random.normal(0, 0.1, len(data['p_vals']))
            fig.add_trace(go.Scatter(
                x=x_pos, y=data['p_vals'], mode='markers', name=tissue,
                marker=dict(color=data['color'], size=6, opacity=0.7),
                hovertemplate=f"<b>{tissue}</b><br>P-value: %{{y:.2e}}<extra></extra>"
            ), row=2, col=1)
    
    # Add significance lines
    fig.add_hline(y=0.1, line_dash="dash", line_color="green", 
                 annotation_text="Good fit (p > 0.1)", row=2, col=1)
    fig.add_hline(y=0.05, line_dash="dash", line_color="orange", 
                 annotation_text="Marginal (p > 0.05)", row=2, col=1)
    
    # Plot 4: Statistical comparisons between tissues
    if len(tissues) >= 2:
        comparisons = []
        p_values = []
        
        for i in range(len(tissues)):
            for j in range(i+1, len(tissues)):
                if len(comparisons) >= 8:  # Limit for readability
                    break
                    
                tissue1, tissue2 = tissues[i], tissues[j]
                alphas1 = tissue_data[tissue1]['alphas']
                alphas2 = tissue_data[tissue2]['alphas']
                
                if len(alphas1) > 0 and len(alphas2) > 0:
                    try:
                        from scipy.stats import mannwhitneyu
                        _, p_val = mannwhitneyu(alphas1, alphas2, alternative='two-sided')
                        comparisons.append(f"{tissue1}\nvs\n{tissue2}")
                        p_values.append(p_val)
                    except:
                        continue
        
        if comparisons:
            fig.add_trace(go.Bar(
                x=comparisons, y=p_values, 
                marker_color=['red' if p < 0.05 else 'lightblue' for p in p_values], 
                opacity=0.8, name="Mann-Whitney U",
                hovertemplate="<b>%{x}</b><br>p-value: %{y:.2e}<br>" +
                             "Significant: %{text}<extra></extra>",
                text=['Yes' if p < 0.05 else 'No' for p in p_values]
            ), row=2, col=2)
            
            fig.add_hline(y=0.05, line_dash="dash", line_color="red", 
                         annotation_text="Significance (p = 0.05)", row=2, col=2)
    
    # Clean layout
    fig.update_layout(
        title={
            'text': "📊 Tissue Comparison - Clean & Organized",
            'font': {'size': 20, 'color': '#2c3e50'}
        },
        height=800,
        showlegend=True,
        legend=dict(
            orientation="h", yanchor="bottom", y=1.02, xanchor="center", x=0.5,
            bgcolor="rgba(255,255,255,0.9)", bordercolor="#e1e8ed", borderwidth=1
        ),
        plot_bgcolor='white',
        paper_bgcolor='#fafbfc'
    )
    
    # Clear axis labels
    fig.update_xaxes(title_text="Tissues", row=1, col=1)
    fig.update_yaxes(title_text="Alpha (Power-Law Exponent)", row=1, col=1)
    
    fig.update_xaxes(title_text="Tissues", row=1, col=2) 
    fig.update_yaxes(title_text="R² (Coefficient of Determination)", row=1, col=2)
    
    fig.update_xaxes(title_text="Tissue Index", row=2, col=1)
    fig.update_yaxes(title_text="P-values (log scale)", type="log", row=2, col=1)
    
    fig.update_xaxes(title_text="Tissue Comparisons", row=2, col=2)
    fig.update_yaxes(title_text="Statistical Test P-values", row=2, col=2)
    
    return fig


def _create_ks_visualization_report(analyses: Dict[str, Dict]) -> go.Figure:
    """Create detailed KS statistic visualization for all tissues."""
    tissues = list(analyses.keys())
    n_tissues = len(tissues)
    
    if n_tissues < 1:
        fig = go.Figure()
        fig.add_annotation(text="No tissues available for KS visualization", 
                          xref="paper", yref="paper", x=0.5, y=0.5,
                          showarrow=False, font=dict(size=16))
        fig.update_layout(title="KS Statistic Visualization", height=400)
        return fig
    
    # Define limited tissues early for performance limits
    MAX_TISSUES = 6  # Limit to prevent exponential complexity
    MAX_PAIRS_PER_TISSUE = 5  # Limit pairs per tissue
    
    limited_tissues = tissues[:MAX_TISSUES] if len(tissues) > MAX_TISSUES else tissues
    if len(tissues) > MAX_TISSUES:
        print(f"   ⚠️  Limiting tissue comparison to first {MAX_TISSUES} tissues for performance")
    
    # Collect and calculate metrics for each tissue
    for i, tissue in enumerate(tissues):
        blocks = analyses[tissue]["blocks"]
        alphas = np.array([b.alpha for b in blocks if np.isfinite(b.alpha)])
        ks_vals = np.array([b.ks_D_tail for b in blocks if np.isfinite(b.ks_D_tail)])
        r2_vals = np.array([b.R2_k_loglog for b in blocks if np.isfinite(b.R2_k_loglog)])
        p_vals = np.array([b.ks_pvalue for b in blocks if np.isfinite(b.ks_pvalue)])
        
        all_alphas_by_tissue[tissue] = alphas
        all_ks_by_tissue[tissue] = ks_vals
        
        # Calculate enhanced metrics
        alpha_metrics = _calculate_enhanced_metrics(alphas)
        tissue_metrics[tissue] = alpha_metrics
        
        # Enhanced color and styling
        color = colors[i % len(colors)]
        marker_symbol = marker_styles[i % len(marker_styles)]
        
        # Plot 1: Alpha distributions with better styling
        fig.add_trace(go.Histogram(
            x=alphas, name=f"{tissue} (n={len(alphas)})", 
            opacity=0.7, marker_color=color, nbinsx=25,
            hovertemplate="<b>%{fullData.name}</b><br>α: %{x:.3f}<br>Count: %{y}<extra></extra>"
        ), row=1, col=1)
        
        # Plot 2: KS-D distributions with enhanced box plots
        fig.add_trace(go.Box(
            y=ks_vals, name=f"{tissue}", marker_color=color, 
            boxpoints='outliers', jitter=0.3, pointpos=-1.8,
            hovertemplate="<b>%{fullData.name}</b><br>KS-D: %{y:.4f}<extra></extra>"
        ), row=1, col=2)
        
        # Plot 3: R² distributions
        fig.add_trace(go.Box(
            y=r2_vals, name=f"{tissue}", marker_color=color,
            boxpoints='outliers', jitter=0.3, pointpos=-1.8,
            hovertemplate="<b>%{fullData.name}</b><br>R²: %{y:.4f}<extra></extra>"
        ), row=1, col=3)
        
        # Plot 4: P-value distributions
        fig.add_trace(go.Histogram(
            x=p_vals, name=f"{tissue}", opacity=0.7,
            marker_color=color, nbinsx=25,
            hovertemplate="<b>%{fullData.name}</b><br>p-value: " + 
                         "%{customdata}<br>Count: %{y}<extra></extra>",
            customdata=[_format_pvalue_hover(p) for p in p_vals]
        ), row=1, col=4)
        
        # Plot 5: Enhanced metrics comparison with better formatting
        metrics_names = ['Mean', 'Median', 'P5', 'P95', 'F>0.2', 'F>0.4', 'F>0.6', 'F>0.8']
        metrics_keys = ['mean', 'median', 'p5', 'p95', 'frac_02', 'frac_04', 'frac_06', 'frac_08']
        metrics_vals = [alpha_metrics[m] for m in metrics_keys]
        
        fig.add_trace(go.Scatter(
            x=metrics_names, y=metrics_vals, mode='lines+markers',
            name=f"{tissue}", line=dict(color=color, width=3),
            marker=dict(color=color, size=8, symbol=marker_symbol, line=dict(width=2, color='white')),
            hovertemplate="<b>%{fullData.name}</b><br>Metric: %{x}<br>Value: %{y:.4f}<extra></extra>"
        ), row=2, col=1)
        
        # Plot 9: Mean vs Median with enhanced markers
        fig.add_trace(go.Scatter(
            x=[alpha_metrics['mean']], y=[alpha_metrics['median']], 
            mode='markers+text', name=f"{tissue}", 
            marker=dict(color=color, size=15, symbol=marker_symbol, 
                       line=dict(width=2, color='white')),
            text=[tissue], textposition="top center", textfont=dict(size=10),
            hovertemplate="<b>%{text}</b><br>Mean: %{x:.3f}<br>Median: %{y:.3f}<extra></extra>"
        ), row=3, col=1)
        
        # Plot 10: P5 vs P95 with enhanced visualization
        fig.add_trace(go.Scatter(
            x=[alpha_metrics['p5']], y=[alpha_metrics['p95']], 
            mode='markers+text', name=f"{tissue}", 
            marker=dict(color=color, size=15, symbol=marker_symbol,
                       line=dict(width=2, color='white')),
            text=[tissue], textposition="top center", textfont=dict(size=10),
            hovertemplate="<b>%{text}</b><br>P5: %{x:.3f}<br>P95: %{y:.3f}<extra></extra>"
        ), row=3, col=2)
        
        # Plot 11: Fraction analysis with better bars
        fractions = [alpha_metrics['frac_02'], alpha_metrics['frac_04'], 
                    alpha_metrics['frac_06'], alpha_metrics['frac_08']]
        fig.add_trace(go.Bar(
            x=['> 0.2', '> 0.4', '> 0.6', '> 0.8'], y=fractions,
            name=f"{tissue}", marker_color=color, opacity=0.8,
            hovertemplate="<b>%{fullData.name}</b><br>Threshold: %{x}<br>Fraction: %{y:.3f}<extra></extra>"
        ), row=3, col=3)
        
        # Plot 12: CV vs Skewness with enhanced styling
        fig.add_trace(go.Scatter(
            x=[alpha_metrics['cv']], y=[alpha_metrics['skewness']], 
            mode='markers+text', name=f"{tissue}", 
            marker=dict(color=color, size=15, symbol=marker_symbol,
                       line=dict(width=2, color='white')),
            text=[tissue], textposition="top center", textfont=dict(size=10),
            hovertemplate="<b>%{text}</b><br>CV: %{x:.3f}<br>Skewness: %{y:.3f}<extra></extra>"
        ), row=3, col=4)
    
    # Calculate EMD matrix for heatmap (Plot 6) with enhanced visualization
    emd_matrix = np.zeros((n_tissues, n_tissues))
    for i, tissue1 in enumerate(tissues):
        for j, tissue2 in enumerate(tissues):
            if i != j:
                alphas1 = all_alphas_by_tissue[tissue1]
                alphas2 = all_alphas_by_tissue[tissue2]
                emd_matrix[i, j] = _calculate_emd(alphas1, alphas2)
    
    fig.add_trace(go.Heatmap(
        z=emd_matrix, x=tissues, y=tissues, 
        colorscale='RdYlBu_r', name="EMD", showlegend=False,
        hovertemplate="<b>EMD Distance</b><br>%{y} → %{x}<br>Distance: %{z:.4f}<extra></extra>",
        colorbar=dict(title="EMD Distance")  # Removed titleside - not supported
    ), row=2, col=2)
    
    # Enhanced statistical tests (Plot 7) - limit to prevent performance issues  
    p_values_mw = []
    p_values_kw = []
    tissue_pairs = []
    
    # Limit tissue combinations for performance
    max_comparisons = 15  # Limit total comparisons
    comparison_count = 0
    
    for i in range(len(limited_tissues)):
        for j in range(i+1, len(limited_tissues)):
            if comparison_count >= max_comparisons:
                break
                
            tissue1, tissue2 = limited_tissues[i], limited_tissues[j]
            alphas1 = all_alphas_by_tissue[tissue1]
            alphas2 = all_alphas_by_tissue[tissue2]
            tissue_pairs.append(f"{tissue1} vs {tissue2}")
            comparison_count += 1
            
            if len(alphas1) > 0 and len(alphas2) > 0:
                try:
                    _, p_val_mw = mannwhitneyu(alphas1, alphas2, alternative='two-sided')
                    p_values_mw.append(p_val_mw)
                except:
                    p_values_mw.append(np.nan)
                
                try:
                    # Kruskal-Wallis for comparison
                    _, p_val_kw = sstats.kruskal(alphas1, alphas2)
                    p_values_kw.append(p_val_kw)
                except:
                    p_values_kw.append(np.nan)
            else:
                p_values_mw.append(np.nan)
                p_values_kw.append(np.nan)
        if comparison_count >= max_comparisons:
            break
    
    if p_values_mw:
        # Create bar plot for statistical tests
        x_pos = list(range(len(tissue_pairs)))
        fig.add_trace(go.Bar(
            x=tissue_pairs, y=p_values_mw, name="Mann-Whitney U", 
            marker_color='lightblue', opacity=0.7,
            hovertemplate="<b>%{x}</b><br>Mann-Whitney p: " + 
                         "%{customdata}<extra></extra>",
            customdata=[_format_pvalue_hover(p) for p in p_values_mw]
        ), row=2, col=3)
        
        fig.add_trace(go.Bar(
            x=tissue_pairs, y=p_values_kw, name="Kruskal-Wallis", 
            marker_color='lightcoral', opacity=0.7,
            hovertemplate="<b>%{x}</b><br>Kruskal-Wallis p: " + 
                         "%{customdata}<extra></extra>",
            customdata=[_format_pvalue_hover(p) for p in p_values_kw]
        ), row=2, col=3)
        
        # Add significance lines
        fig.add_hline(y=0.05, line_dash="dash", line_color="red", 
                     annotation_text="p = 0.05", row=2, col=3)
        fig.add_hline(y=0.01, line_dash="dash", line_color="darkred", 
                     annotation_text="p = 0.01", row=2, col=3)
    
    # Enhanced effect sizes (Plot 8) - limit for performance
    effect_sizes = []
    effect_labels = []
    effect_count = 0
    
    for i in range(len(limited_tissues)):
        for j in range(i+1, len(limited_tissues)):
            if effect_count >= max_comparisons:
                break
                
            tissue1, tissue2 = limited_tissues[i], limited_tissues[j]
            alphas1 = all_alphas_by_tissue[tissue1]
            alphas2 = all_alphas_by_tissue[tissue2]
            effect_count += 1
            
            if len(alphas1) > 0 and len(alphas2) > 0:
                # Cohen's d
                pooled_std = np.sqrt(((len(alphas1) - 1) * np.var(alphas1, ddof=1) + 
                                     (len(alphas2) - 1) * np.var(alphas2, ddof=1)) / 
                                    (len(alphas1) + len(alphas2) - 2))
                if pooled_std > 0:
                    cohens_d = (np.mean(alphas1) - np.mean(alphas2)) / pooled_std
                    effect_sizes.append(abs(cohens_d))
                    effect_labels.append(f"{tissue1} vs {tissue2}")
        if effect_count >= max_comparisons:
            break
    
    if effect_sizes:
        fig.add_trace(go.Bar(
            x=effect_labels, y=effect_sizes, name="Effect Size (|Cohen's d|)", 
            marker_color='mediumpurple', opacity=0.8,
            hovertemplate="<b>%{x}</b><br>Effect Size: %{y:.3f}<extra></extra>"
        ), row=2, col=4)
        
        # Add effect size interpretation lines
        fig.add_hline(y=0.2, line_dash="dash", line_color="green", 
                     annotation_text="Small (0.2)", row=2, col=4)
        fig.add_hline(y=0.5, line_dash="dash", line_color="orange", 
                     annotation_text="Medium (0.5)", row=2, col=4)
        fig.add_hline(y=0.8, line_dash="dash", line_color="red", 
                     annotation_text="Large (0.8)", row=2, col=4)
    
    # Enhanced pair-wise comparisons (Plots 13-15) with better organization
    tissue_pair_data = {}
    
    for tissue in limited_tissues:
        blocks = analyses[tissue]["blocks"]
        pairs = [(b.region_a, b.region_b) for b in blocks]
        unique_pairs = list(set(pairs))
        
        # Sort pairs alphabetically for consistency and limit
        unique_pairs.sort()
        
        for pair in unique_pairs[:MAX_PAIRS_PER_TISSUE]:  # Limit pairs per tissue
            pair_blocks = [b for b in blocks if (b.region_a, b.region_b) == pair]
            if pair_blocks:
                pair_label = f"{pair[0]}→{pair[1]}"
                
                alphas = [b.alpha for b in pair_blocks if np.isfinite(b.alpha)]
                ks_vals = [b.ks_D_tail for b in pair_blocks if np.isfinite(b.ks_D_tail)]
                r2_vals = [b.R2_k_loglog for b in pair_blocks if np.isfinite(b.R2_k_loglog)]
                
                if pair_label not in tissue_pair_data:
                    tissue_pair_data[pair_label] = {'tissues': [], 'alphas': [], 'ks': [], 'r2': []}
                
                tissue_pair_data[pair_label]['tissues'].append(tissue)
                tissue_pair_data[pair_label]['alphas'].extend(alphas)
                tissue_pair_data[pair_label]['ks'].extend(ks_vals)
                tissue_pair_data[pair_label]['r2'].extend(r2_vals)
    
    # Plot organized pair-wise data
    pair_labels = sorted(tissue_pair_data.keys())[:10]  # Top 10 pairs
    
    for i, pair_label in enumerate(pair_labels):
        data = tissue_pair_data[pair_label]
        color = colors[i % len(colors)]
        
        if data['alphas']:
            fig.add_trace(go.Box(
                x=[pair_label] * len(data['alphas']), y=data['alphas'], 
                name=f"{pair_label}", marker_color=color, opacity=0.7,
                hovertemplate="<b>%{x}</b><br>α: %{y:.3f}<extra></extra>"
            ), row=4, col=1)
            
        if data['ks']:
            fig.add_trace(go.Box(
                x=[pair_label] * len(data['ks']), y=data['ks'], 
                name=f"{pair_label}", marker_color=color, opacity=0.7,
                hovertemplate="<b>%{x}</b><br>KS-D: %{y:.4f}<extra></extra>"
            ), row=4, col=2)
            
        if data['r2']:
            fig.add_trace(go.Box(
                x=[pair_label] * len(data['r2']), y=data['r2'], 
                name=f"{pair_label}", marker_color=color, opacity=0.7,
                hovertemplate="<b>%{x}</b><br>R²: %{y:.4f}<extra></extra>"
            ), row=4, col=3)
    
    # Combined metrics summary (Plot 16)
    summary_text = f"""
    <b>📊 Analysis Summary</b><br>
    • Tissues analyzed: {len(tissues)}<br>
    • Total comparisons: {len(tissue_pairs)}<br>
    • Significant differences (p<0.05): {sum(1 for p in p_values_mw if not np.isnan(p) and p < 0.05)}<br>
    • Large effect sizes (d>0.8): {sum(1 for es in effect_sizes if es > 0.8)}<br>
    • Mean EMD distance: {np.mean(emd_matrix[emd_matrix > 0]):.4f}<br>
    """
    
    fig.add_annotation(
        text=summary_text, xref="x domain", yref="y domain",
        x=0.1, y=0.9, font=dict(family="monospace", size=12),
        bgcolor="lightgray", bordercolor="black", borderwidth=1,
        row=4, col=4
    )
    
    # Add reference lines to scatter plots
    # Mean vs Median diagonal
    min_val = min([min(all_alphas_by_tissue[t]) for t in tissues if len(all_alphas_by_tissue[t]) > 0])
    max_val = max([max(all_alphas_by_tissue[t]) for t in tissues if len(all_alphas_by_tissue[t]) > 0])
    fig.add_shape(type="line", x0=min_val, y0=min_val, x1=max_val, y1=max_val, 
                  line=dict(color="gray", dash="dash"), row=3, col=1)
    
    # Update layout with enhanced styling
    fig.update_layout(
        title={
            'text': "🔬 Enhanced Tissue Comparison Analysis<br><sub>Comprehensive Statistical Comparison of Power-Law Parameters Across Tissues</sub>",
            'x': 0.5,
            'xanchor': 'center',
            'font': {'size': 20}
        },
        width=1900, height=1700,
        showlegend=True,
        font=dict(family="Arial, sans-serif", size=10),
        plot_bgcolor='white',
        paper_bgcolor='white'
    )
    
    # Update axes labels
    fig.update_xaxes(title_text="Power-law Exponent (α)", row=1, col=1)
    fig.update_xaxes(title_text="Tissue", row=1, col=2)
    fig.update_xaxes(title_text="Tissue", row=1, col=3)
    fig.update_xaxes(title_text="P-value", row=1, col=4)
    fig.update_xaxes(title_text="Statistical Measure", row=2, col=1)
    fig.update_xaxes(title_text="Tissue Comparisons", row=2, col=3)
    fig.update_xaxes(title_text="Tissue Comparisons", row=2, col=4)
    fig.update_xaxes(title_text="Mean α", row=3, col=1)
    fig.update_xaxes(title_text="5th Percentile", row=3, col=2)
    fig.update_xaxes(title_text="Threshold", row=3, col=3)
    fig.update_xaxes(title_text="Coefficient of Variation", row=3, col=4)
    fig.update_xaxes(title_text="Region Pairs", row=4, col=1)
    fig.update_xaxes(title_text="Region Pairs", row=4, col=2)
    fig.update_xaxes(title_text="Region Pairs", row=4, col=3)
    
    fig.update_yaxes(title_text="Frequency", row=1, col=1)
    fig.update_yaxes(title_text="KS-D Statistic", row=1, col=2)
    fig.update_yaxes(title_text="R² Value", row=1, col=3)
    fig.update_yaxes(title_text="Frequency", row=1, col=4)
    fig.update_yaxes(title_text="Metric Value", row=2, col=1)
    fig.update_yaxes(title_text="P-value", row=2, col=3)
    fig.update_yaxes(title_text="Effect Size", row=2, col=4)
    fig.update_yaxes(title_text="Median α", row=3, col=1)
    fig.update_yaxes(title_text="95th Percentile", row=3, col=2)
    fig.update_yaxes(title_text="Fraction", row=3, col=3)
    fig.update_yaxes(title_text="Skewness", row=3, col=4)
    fig.update_yaxes(title_text="α Value", row=4, col=1)
    fig.update_yaxes(title_text="KS-D Value", row=4, col=2)
    fig.update_yaxes(title_text="R² Value", row=4, col=3)
    
    return fig

def _create_ks_visualization_report(analyses: Dict[str, Dict]) -> go.Figure:
    """Create detailed KS statistic visualization for all tissues."""
    tissues = list(analyses.keys())
    n_tissues = len(tissues)
    
    # Calculate grid size
    cols = min(3, n_tissues)
    rows = math.ceil(n_tissues / cols)
    
    subplot_titles = [f"KS Analysis: {tissue}" for tissue in tissues]
    
    fig = make_subplots(
        rows=rows, cols=cols,
        subplot_titles=subplot_titles,
        specs=[[{"secondary_y": False} for _ in range(cols)] for _ in range(rows)],
        vertical_spacing=0.1, horizontal_spacing=0.1
    )
    
    for idx, tissue in enumerate(tissues):
        row = idx // cols + 1
        col = idx % cols + 1
        
        if row > rows:
            break
            
        blocks = analyses[tissue]["blocks"]
        # Take representative block for KS visualization
        if blocks:
            representative_block = blocks[0]  # Use first block
            if hasattr(representative_block, 'v_data') and len(representative_block.v_data) > 0:
                x_data = np.array(representative_block.v_data)
                xmin = representative_block.xmin
                alpha = representative_block.alpha
                
                # Create mini KS visualization
                x_sorted = np.sort(x_data[x_data >= xmin])
                if len(x_sorted) > 0:
                    n = len(x_sorted)
                    empirical_cdf = np.arange(1, n + 1) / n
                    theoretical_cdf = 1 - (x_sorted / xmin) ** (1 - alpha)
                    ks_differences = np.abs(empirical_cdf - theoretical_cdf)
                    
                    # Plot empirical vs theoretical CDF
                    fig.add_trace(go.Scatter(x=x_sorted, y=empirical_cdf, name=f"{tissue} Empirical", 
                                            line=dict(color='blue', width=2), showlegend=(idx==0)), row=row, col=col)
                    fig.add_trace(go.Scatter(x=x_sorted, y=theoretical_cdf, name=f"{tissue} Theoretical", 
                                            line=dict(color='red', width=2), showlegend=(idx==0)), row=row, col=col)
                    
                    # Highlight maximum difference
                    ks_stat = np.max(ks_differences)
                    ks_idx = np.argmax(ks_differences)
                    fig.add_trace(go.Scatter(x=[x_sorted[ks_idx]], y=[empirical_cdf[ks_idx]], 
                                            name=f"KS={ks_stat:.3f}", mode='markers',
                                            marker=dict(color='orange', size=8), showlegend=(idx==0)), row=row, col=col)
        
        # Update axes
        fig.update_xaxes(title_text="Value" if row == rows else "", row=row, col=col)
        fig.update_yaxes(title_text="CDF" if col == 1 else "", row=row, col=col)
    
    fig.update_layout(title="KS Statistic Visualization Across Tissues", 
                      width=1400, height=400*rows)
    
    return fig

def _create_r2_log_scale_analysis(analyses: Dict[str, Dict]) -> go.Figure:
    """Create R² analysis with log10 scale, linear and logarithmic binning."""
    tissues = list(analyses.keys())[:6]  # Limit for clarity
    
    if len(tissues) < 1:
        fig = go.Figure()
        fig.add_annotation(text="No tissues available for R² analysis", 
                          xref="paper", yref="paper", x=0.5, y=0.5,
                          showarrow=False, font=dict(size=16))
        fig.update_layout(title="R² Log Scale Analysis", height=400)
        return fig

    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b']
    
    # Create 1x2 subplot layout - only two plots as requested
    fig = make_subplots(
        rows=1, cols=2,
        subplot_titles=(
            "📊 R² Distribution - Linear Binning<br><sub>Log10 scale with uniform bin spacing</sub>",
            "📈 R² Distribution - Logarithmic Binning<br><sub>Log10 scale with exponential bin jumps (10^x)</sub>"
        ),
        horizontal_spacing=0.15
    )
    
    # Collect all R² data
    all_r2_data = {}
    for i, tissue in enumerate(tissues):
        blocks = analyses[tissue]["blocks"]
        r2_vals = np.array([b.R2_k_loglog for b in blocks if np.isfinite(b.R2_k_loglog) and b.R2_k_loglog > 0])
        
        if len(r2_vals) > 0:
            all_r2_data[tissue] = {
                'r2_vals': r2_vals,
                'color': colors[i % len(colors)]
            }
    
    if not all_r2_data:
        fig.add_annotation(text="No valid R² data available", 
                          xref="paper", yref="paper", x=0.5, y=0.5,
                          showarrow=False, font=dict(size=16))
        return fig
    
    # Plot 1: Linear binning on log scale
    for tissue, data in all_r2_data.items():
        # Convert to log10 for plotting
        log_r2 = np.log10(data['r2_vals'])
        
        fig.add_trace(go.Histogram(
            x=log_r2,
            name=f"{tissue}",
            marker_color=data['color'],
            opacity=0.7,
            nbinsx=50,  # Linear binning - uniform spacing
            hovertemplate=f"<b>{tissue}</b><br>" +
                         "Log10(R²): %{x:.5f}<br>" +
                         "R²: %{customdata:.2e}<br>" +
                         "Count: %{y}<extra></extra>",
            customdata=data['r2_vals']  # Original R² values for hover
        ), row=1, col=1)
    
    # Plot 2: Logarithmic binning on log scale  
    # Create logarithmic bins (powers of 10)
    log_bin_edges = np.arange(-10, 1, 0.5)  # From 10^-10 to 10^0 in 0.5 increments
    
    for tissue, data in all_r2_data.items():
        log_r2 = np.log10(data['r2_vals'])
        
        # Create histogram with logarithmic bins
        counts, bin_edges = np.histogram(log_r2, bins=log_bin_edges)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        
        # Only show non-zero counts
        valid_indices = counts > 0
        if np.any(valid_indices):
            fig.add_trace(go.Bar(
                x=bin_centers[valid_indices],
                y=counts[valid_indices],
                name=f"{tissue}",
                marker_color=data['color'],
                opacity=0.7,
                width=0.4,  # Bar width for logarithmic bins
                hovertemplate=f"<b>{tissue}</b><br>" +
                             "Log10(R²): %{x:.1f}<br>" +
                             "R² Range: %{customdata}<br>" +
                             "Count: %{y}<extra></extra>",
                customdata=[f"10^{edge:.1f} - 10^{bin_edges[i+1]:.1f}" 
                           for i, edge in enumerate(bin_edges[:-1]) if counts[i] > 0]
            ), row=1, col=2)
    
    # Update layout with log10 scale and scientific notation
    fig.update_layout(
        title={
            'text': "📊 R² Analysis - Log10 Scale with Linear & Logarithmic Binning",
            'font': {'size': 18, 'color': '#2c3e50'},
            'y': 0.92  # Move title down slightly to avoid overlap
        },
        height=600,  # Increased height for better spacing
        showlegend=True,
        legend=dict(
            orientation="h", 
            yanchor="bottom", 
            y=-0.18, # Move legend below the plots
            xanchor="center", 
            x=0.5,
            bgcolor="rgba(255,255,255,0.9)", 
            bordercolor="#e1e8ed", 
            borderwidth=1
        ),
        plot_bgcolor='white',
        paper_bgcolor='#fafbfc',
        margin=dict(t=100, b=100, l=60, r=60)  # Increased margins
    )
    
    # Configure axes with proper scientific notation
    fig.update_xaxes(
        title_text="Log10(R²)",
        row=1, col=1,
        tickformat=".1f",  # Show as decimal
        tickvals=np.arange(-8, 1, 1),  # Major ticks at integers
        ticktext=[f"10^{i}" for i in range(-8, 1)]  # Scientific notation labels
    )
    fig.update_yaxes(title_text="Count (Linear Bins)", row=1, col=1)
    
    fig.update_xaxes(
        title_text="Log10(R²)", 
        row=1, col=2,
        tickformat=".1f",
        tickvals=np.arange(-8, 1, 1),
        ticktext=[f"10^{i}" for i in range(-8, 1)]
    )
    fig.update_yaxes(title_text="Count (Log Bins)", row=1, col=2)
    
    # Add reference lines for quality thresholds
    quality_thresholds = [0.8, 0.9, 0.95]
    for threshold in quality_thresholds:
        log_threshold = np.log10(threshold)
        fig.add_vline(
            x=log_threshold, 
            line_dash="dash", 
            line_color="red" if threshold == 0.8 else "orange" if threshold == 0.9 else "green",
            annotation_text=f"R² = {threshold}",
            row=1, col=1
        )
        fig.add_vline(
            x=log_threshold, 
            line_dash="dash", 
            line_color="red" if threshold == 0.8 else "orange" if threshold == 0.9 else "green",
            annotation_text=f"R² = {threshold}",
            row=1, col=2
        )
    
    return fig


def _create_pvalue_scientific_notation_plot(analyses: Dict[str, Dict]) -> go.Figure:
    """Create p-value visualization with proper scientific notation for small numbers."""
    tissues = list(analyses.keys())[:6]  # Limit for clarity
    
    if len(tissues) < 1:
        fig = go.Figure()
        fig.add_annotation(text="No tissues available for p-value analysis", 
                          xref="paper", yref="paper", x=0.5, y=0.5,
                          showarrow=False, font=dict(size=16))
        fig.update_layout(title="P-value Scientific Notation Analysis", height=400)
        return fig

    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b']
    
    # Create 2x2 subplot layout for comprehensive p-value analysis
    fig = make_subplots(
        rows=2, cols=2,
        subplot_titles=(
            "📊 P-values by Tissue<br><sub>Distribution with scientific notation</sub>",
            "📈 P-value Significance Assessment<br><sub>Scatter plot with significance lines</sub>", 
            "📉 Cumulative P-value Distribution<br><sub>ECDF showing significance patterns</sub>",
            "🔬 Statistical Summary<br><sub>Significance counts by tissue</sub>"
        ),
        vertical_spacing=0.25, horizontal_spacing=0.15  # Increased spacing for titles
    )
    
    # Collect p-value data
    tissue_pvalue_data = {}
    for i, tissue in enumerate(tissues):
        blocks = analyses[tissue]["blocks"]
        p_vals = np.array([b.ks_pvalue for b in blocks if np.isfinite(b.ks_pvalue) and b.ks_pvalue > 0])
        
        if len(p_vals) > 0:
            tissue_pvalue_data[tissue] = {
                'p_vals': p_vals,
                'color': colors[i % len(colors)]
            }
    
    if not tissue_pvalue_data:
        fig.add_annotation(text="No valid p-value data available", 
                          xref="paper", yref="paper", x=0.5, y=0.5,
                          showarrow=False, font=dict(size=16))
        return fig
    
    # Plot 1: Box plots of p-values by tissue
    for tissue, data in tissue_pvalue_data.items():
        fig.add_trace(go.Box(
            y=data['p_vals'],
            name=tissue,
            marker_color=data['color'],
            opacity=0.7,
            boxmean=True,
            hovertemplate=f"<b>{tissue}</b><br>" +
                         "P-value: %{y:.2e}<br>" +
                         f"Median: {np.median(data['p_vals']):.2e}<extra></extra>"
        ), row=1, col=1)
    
    # Plot 2: Scatter plot with jitter
    for i, (tissue, data) in enumerate(tissue_pvalue_data.items()):
        # Add jitter for better visualization
        x_pos = np.full(len(data['p_vals']), i) + np.random.normal(0, 0.1, len(data['p_vals']))
        
        fig.add_trace(go.Scatter(
            x=x_pos,
            y=data['p_vals'],
            mode='markers',
            name=tissue,
            marker=dict(color=data['color'], size=4, opacity=0.6),
            hovertemplate=f"<b>{tissue}</b><br>P-value: %{{y:.2e}}<extra></extra>"
        ), row=1, col=2)
    
    # Plot 3: Cumulative distribution
    for tissue, data in tissue_pvalue_data.items():
        sorted_pvals = np.sort(data['p_vals'])
        y_vals = np.arange(1, len(sorted_pvals) + 1) / len(sorted_pvals)
        
        fig.add_trace(go.Scatter(
            x=sorted_pvals,
            y=y_vals,
            mode='lines',
            name=tissue,
            line=dict(color=data['color'], width=2),
            hovertemplate=f"<b>{tissue}</b><br>" +
                         "P-value: %{x:.2e}<br>" +
                         "Cumulative fraction: %{y:.3f}<extra></extra>"
        ), row=2, col=1)
    
    # Plot 4: Significance summary bars
    significance_levels = [0.001, 0.01, 0.05, 0.1]
    significance_data = {level: [] for level in significance_levels}
    tissue_names = list(tissue_pvalue_data.keys())
    
    for tissue in tissue_names:
        p_vals = tissue_pvalue_data[tissue]['p_vals']
        total = len(p_vals)
        
        for level in significance_levels:
            significant_count = np.sum(p_vals < level)
            significance_data[level].append(significant_count / total if total > 0 else 0)
    
    for level in significance_levels:
        fig.add_trace(go.Bar(
            x=tissue_names,
            y=significance_data[level],
            name=f"p < {level:.3f}",
            opacity=0.8,
            hovertemplate=f"<b>%{{x}}</b><br>Fraction p < {level:.3f}: %{{y:.3f}}<extra></extra>"
        ), row=2, col=2)
    
    # Add significance reference lines
    significance_thresholds = [0.001, 0.01, 0.05, 0.1]
    colors_ref = ['darkred', 'red', 'orange', 'green']
    
    for threshold, color in zip(significance_thresholds, colors_ref):
        # Add to box plot
        fig.add_hline(
            y=threshold, 
            line_dash="dash", 
            line_color=color,
            annotation_text=f"p = {threshold:.3f}",
            row=1, col=1
        )
        # Add to scatter plot  
        fig.add_hline(
            y=threshold, 
            line_dash="dash", 
            line_color=color,
            annotation_text=f"p = {threshold:.3f}",
            row=1, col=2
        )
        # Add to cumulative plot
        fig.add_vline(
            x=threshold, 
            line_dash="dash", 
            line_color=color,
            annotation_text=f"p = {threshold:.3f}",
            row=2, col=1
        )
    
    # Update layout with scientific notation
    fig.update_layout(
        title={
            'text': "📊 P-value Analysis with Scientific Notation",
            'font': {'size': 18, 'color': '#2c3e50'},
            'y': 0.95  # Move title down to avoid overlap
        },
        height=900,  # Increased height for better spacing
        showlegend=True,
        legend=dict(
            orientation="h", 
            yanchor="bottom", 
            y=-0.15,  # Move legend below the plots
            xanchor="center", 
            x=0.5,
            bgcolor="rgba(255,255,255,0.9)", 
            bordercolor="#e1e8ed", 
            borderwidth=1
        ),
        plot_bgcolor='white',
        paper_bgcolor='#fafbfc',
        margin=dict(t=120, b=120, l=60, r=60)  # Increased margins for titles and legend
    )
    
    # Configure axes with scientific notation for p-values
    fig.update_xaxes(title_text="Tissues", row=1, col=1)
    fig.update_yaxes(title_text="P-values", type="log", exponentformat="e", row=1, col=1)
    
    fig.update_xaxes(title_text="Tissue Index", row=1, col=2)
    fig.update_yaxes(title_text="P-values", type="log", exponentformat="e", row=1, col=2)
    
    fig.update_xaxes(title_text="P-values", type="log", exponentformat="e", row=2, col=1)
    fig.update_yaxes(title_text="Cumulative Probability", row=2, col=1)
    
    fig.update_xaxes(title_text="Tissues", row=2, col=2)
    fig.update_yaxes(title_text="Fraction Significant", row=2, col=2)
    
    return fig
    """Create clean R² analysis with focused 2x2 layout."""
    tissues = list(analyses.keys())[:6]  # Limit tissues for clarity
    
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b']
    
    # Clean 2x2 subplot grid
    fig = make_subplots(
        rows=2, cols=2,
        subplot_titles=(
            "📊 R² Linear vs Log-Log Fit<br><sub>Comparing linear and log-log fit quality</sub>",
            "📈 R² Distribution by Tissue<br><sub>Box plots of R² values across tissues</sub>", 
            "📉 R² Quality Assessment<br><sub>High-quality fits (R² > 0.8) by tissue</sub>",
            "📦 R² Correlation Analysis<br><sub>Relationship between fit types</sub>"
        ),
        vertical_spacing=0.15, horizontal_spacing=0.12
    )
    
    # Collect data
    tissue_data = {}
    for i, tissue in enumerate(tissues):
        blocks = analyses[tissue]["blocks"]
        r2_linear = np.array([b.R2_k_linear for b in blocks if np.isfinite(b.R2_k_linear)])
        r2_loglog = np.array([b.R2_k_loglog for b in blocks if np.isfinite(b.R2_k_loglog)])
        
        tissue_data[tissue] = {
            'r2_linear': r2_linear,
            'r2_loglog': r2_loglog,
            'color': colors[i % len(colors)]
        }
    
    # Plot 1: Linear vs Log-Log scatter comparison
    for tissue, data in tissue_data.items():
        if len(data['r2_linear']) > 0 and len(data['r2_loglog']) > 0:
            # Match lengths for comparison
            min_len = min(len(data['r2_linear']), len(data['r2_loglog']))
            fig.add_trace(go.Scatter(
                x=data['r2_linear'][:min_len], 
                y=data['r2_loglog'][:min_len],
                mode='markers', name=tissue,
                marker=dict(color=data['color'], size=6, opacity=0.7),
                hovertemplate=f"<b>{tissue}</b><br>R² Linear: %{{x:.3f}}<br>R² Log-Log: %{{y:.3f}}<extra></extra>"
            ), row=1, col=1)
    
    # Add diagonal reference line
    fig.add_shape(type="line", x0=0, y0=0, x1=1, y1=1, 
                  line=dict(color="gray", dash="dash"), row=1, col=1)
    
    # Plot 2: R² distributions by tissue
    for tissue, data in tissue_data.items():
        if len(data['r2_loglog']) > 0:
            fig.add_trace(go.Box(
                y=data['r2_loglog'], name=tissue, 
                marker_color=data['color'], opacity=0.7, boxmean=True,
                hovertemplate=f"<b>{tissue}</b><br>R²: %{{y:.3f}}<br>Median: {np.median(data['r2_loglog']):.3f}<extra></extra>"
            ), row=1, col=2)
    
    # Plot 3: High-quality fits assessment
    high_quality_counts = []
    tissue_names = []
    for tissue, data in tissue_data.items():
        if len(data['r2_loglog']) > 0:
            high_quality = np.sum(data['r2_loglog'] > 0.8)
            total = len(data['r2_loglog'])
            high_quality_counts.append(high_quality / total if total > 0 else 0)
            tissue_names.append(tissue)
    
    if tissue_names:
        fig.add_trace(go.Bar(
            x=tissue_names, y=high_quality_counts,
            marker=dict(color=[tissue_data[t]['color'] for t in tissue_names]),
            opacity=0.8, name="High Quality Fits",
            hovertemplate="<b>%{x}</b><br>Fraction R² > 0.8: %{y:.3f}<extra></extra>"
        ), row=2, col=1)
    
    # Add quality threshold line
    fig.add_hline(y=0.5, line_dash="dash", line_color="green", 
                 annotation_text="50% threshold", row=2, col=1)
    
    # Plot 4: Correlation analysis - aggregate view
    all_linear = []
    all_loglog = []
    tissue_colors = []
    
    for tissue, data in tissue_data.items():
        if len(data['r2_linear']) > 0 and len(data['r2_loglog']) > 0:
            min_len = min(len(data['r2_linear']), len(data['r2_loglog']))
            all_linear.extend(data['r2_linear'][:min_len])
            all_loglog.extend(data['r2_loglog'][:min_len])
            tissue_colors.extend([data['color']] * min_len)
    
    if all_linear:
        # Create 2D histogram for correlation visualization
        fig.add_trace(go.Histogram2d(
            x=all_linear, y=all_loglog, 
            colorscale='Viridis', opacity=0.8,
            nbinsx=20, nbinsy=20,
            hovertemplate="Linear R²: %{x:.3f}<br>Log-Log R²: %{y:.3f}<br>Count: %{z}<extra></extra>"
        ), row=2, col=2)
        
        # Add correlation line
        if len(all_linear) > 1:
            correlation = np.corrcoef(all_linear, all_loglog)[0, 1]
            fig.add_annotation(
                text=f"Correlation: {correlation:.3f}", 
                xref="x4", yref="y4", x=0.1, y=0.9, 
                showarrow=False, font=dict(size=12, color="white"),
                bgcolor="rgba(0,0,0,0.5)", bordercolor="white", borderwidth=1
            )
    
    # Clean layout
    fig.update_layout(
        title={
            'text': "� R² Analysis - Clean & Focused",
            'font': {'size': 20, 'color': '#2c3e50'}
        },
        height=800,
        showlegend=True,
        legend=dict(
            orientation="h", yanchor="bottom", y=1.02, xanchor="center", x=0.5,
            bgcolor="rgba(255,255,255,0.9)", bordercolor="#e1e8ed", borderwidth=1
        ),
        plot_bgcolor='white',
        paper_bgcolor='#fafbfc'
    )
    
    # Clear axis labels
    fig.update_xaxes(title_text="R² Linear Fit", row=1, col=1)
    fig.update_yaxes(title_text="R² Log-Log Fit", row=1, col=1)
    
    fig.update_xaxes(title_text="Tissues", row=1, col=2)
    fig.update_yaxes(title_text="R² Values", row=1, col=2)
    
    fig.update_xaxes(title_text="Tissues", row=2, col=1)
    fig.update_yaxes(title_text="Fraction High Quality (R² > 0.8)", row=2, col=1)
    
    fig.update_xaxes(title_text="R² Linear Fit", row=2, col=2)
    fig.update_yaxes(title_text="R² Log-Log Fit", row=2, col=2)
    
    return fig


def _create_individual_statistical_test_plots(analyses: Dict[str, Dict], out_dir: str) -> int:
    
    all_r2_data = {
        'linear': [], 'loglog': [], 'semilog': [], 
        'tissues': [], 'alphas': [], 'pairs': []
    }
    tissue_r2_stats = {}
    
    # Collect R² data for all tissues with enhanced organization
    for i, tissue in enumerate(tissues):
        blocks = analyses[tissue]["blocks"]
        color = colors[i % len(colors)]
        marker_symbol = marker_styles[i % len(marker_styles)]
        
        # Collect data
        r2_linear = [b.R2_k_linear for b in blocks if np.isfinite(b.R2_k_linear)]
        r2_loglog = [b.R2_k_loglog for b in blocks if np.isfinite(b.R2_k_loglog)]
        r2_semilog = [b.R2_k_semilog for b in blocks if np.isfinite(b.R2_k_semilog)]
        alphas = [b.alpha for b in blocks if np.isfinite(b.alpha)]
        pairs = [f"{b.region_a}→{b.region_b}" for b in blocks]
        
        # Store data
        all_r2_data['linear'].extend(r2_linear)
        all_r2_data['loglog'].extend(r2_loglog)
        all_r2_data['semilog'].extend(r2_semilog)
        all_r2_data['tissues'].extend([tissue] * len(r2_linear))
        all_r2_data['alphas'].extend(alphas[:len(r2_linear)])
        all_r2_data['pairs'].extend(pairs[:len(r2_linear)])
        
        # Calculate tissue-specific statistics
        tissue_r2_stats[tissue] = {
            'linear': {'mean': np.mean(r2_linear), 'median': np.median(r2_linear), 'std': np.std(r2_linear)} if r2_linear else None,
            'loglog': {'mean': np.mean(r2_loglog), 'median': np.median(r2_loglog), 'std': np.std(r2_loglog)} if r2_loglog else None,
            'semilog': {'mean': np.mean(r2_semilog), 'median': np.median(r2_semilog), 'std': np.std(r2_semilog)} if r2_semilog else None
        }
        
        # Ensure equal lengths for comparisons
        min_len = min(len(r2_linear), len(r2_loglog), len(r2_semilog))
        if min_len > 0:
            r2_linear = r2_linear[:min_len]
            r2_loglog = r2_loglog[:min_len]
            r2_semilog = r2_semilog[:min_len]
            alphas = alphas[:min_len]
        
        # Enhanced scatter plots with better markers and hover info
        if len(r2_linear) > 0 and len(r2_loglog) > 0:
            fig.add_trace(go.Scatter(
                x=r2_linear, y=r2_loglog, mode='markers',
                name=f"{tissue} (n={len(r2_linear)})", 
                marker=dict(color=color, size=8, symbol=marker_symbol, 
                           line=dict(width=1, color='white'), opacity=0.7),
                hovertemplate="<b>%{fullData.name}</b><br>R² Linear: %{x:.4f}<br>R² Log-Log: %{y:.4f}<extra></extra>"
            ), row=1, col=1)
        
        if len(r2_linear) > 0 and len(r2_semilog) > 0:
            fig.add_trace(go.Scatter(
                x=r2_linear, y=r2_semilog, mode='markers',
                name=f"{tissue}", 
                marker=dict(color=color, size=8, symbol=marker_symbol, 
                           line=dict(width=1, color='white'), opacity=0.7),
                hovertemplate="<b>%{fullData.name}</b><br>R² Linear: %{x:.4f}<br>R² Semi-Log: %{y:.4f}<extra></extra>"
            ), row=1, col=2)
        
        if len(r2_loglog) > 0 and len(r2_semilog) > 0:
            fig.add_trace(go.Scatter(
                x=r2_loglog, y=r2_semilog, mode='markers',
                name=f"{tissue}", 
                marker=dict(color=color, size=8, symbol=marker_symbol, 
                           line=dict(width=1, color='white'), opacity=0.7),
                hovertemplate="<b>%{fullData.name}</b><br>R² Log-Log: %{x:.4f}<br>R² Semi-Log: %{y:.4f}<extra></extra>"
            ), row=1, col=3)
        
        # Enhanced box plots with organized grouping
        if r2_linear:
            fig.add_trace(go.Box(
                y=r2_linear, name=f"{tissue} Linear", 
                marker_color=color, opacity=0.6, boxpoints='outliers',
                hovertemplate="<b>%{fullData.name}</b><br>R²: %{y:.4f}<extra></extra>",
                offsetgroup=f"{tissue}_linear"
            ), row=2, col=1)
            
        if r2_loglog:
            fig.add_trace(go.Box(
                y=r2_loglog, name=f"{tissue} LogLog", 
                marker_color=color, opacity=0.8, boxpoints='outliers',
                hovertemplate="<b>%{fullData.name}</b><br>R²: %{y:.4f}<extra></extra>",
                offsetgroup=f"{tissue}_loglog"
            ), row=2, col=1)
            
        if r2_semilog:
            fig.add_trace(go.Box(
                y=r2_semilog, name=f"{tissue} SemiLog", 
                marker_color=color, opacity=0.4, boxpoints='outliers',
                hovertemplate="<b>%{fullData.name}</b><br>R²: %{y:.4f}<extra></extra>",
                offsetgroup=f"{tissue}_semilog"
            ), row=2, col=1)
        
        # Enhanced R² vs Alpha plot
        if len(r2_loglog) > 0 and len(alphas) > 0:
            fig.add_trace(go.Scatter(
                x=alphas, y=r2_loglog, mode='markers',
                name=f"{tissue} R²-α", 
                marker=dict(color=color, size=8, symbol=marker_symbol,
                           line=dict(width=1, color='white'), opacity=0.7),
                hovertemplate="<b>%{fullData.name}</b><br>α: %{x:.3f}<br>R²: %{y:.4f}<extra></extra>"
            ), row=2, col=3)
    
    # Enhanced correlation matrix (Plot 5)
    r2_types = ['Linear', 'Log-Log', 'Semi-Log']
    r2_keys = ['linear', 'loglog', 'semilog']
    r2_corr_matrix = np.eye(3)
    
    for i, key1 in enumerate(r2_keys):
        for j, key2 in enumerate(r2_keys):
            if i != j:
                vals1 = np.array(all_r2_data[key1])
                vals2 = np.array(all_r2_data[key2])
                if len(vals1) == len(vals2) and len(vals1) > 0:
                    r2_corr_matrix[i, j] = np.corrcoef(vals1, vals2)[0, 1]
    
    fig.add_trace(go.Heatmap(
        z=r2_corr_matrix, x=r2_types, y=r2_types,
        colorscale='RdBu', zmid=0, showlegend=False,
        hovertemplate="<b>Correlation</b><br>%{y} vs %{x}<br>r = %{z:.3f}<extra></extra>",
        colorbar=dict(title="Pearson r")  # Removed titleside - not supported
    ), row=2, col=2)
    
    # Add correlation values as text
    for i in range(3):
        for j in range(3):
            fig.add_annotation(
                x=r2_types[j], y=r2_types[i],
                text=f"{r2_corr_matrix[i, j]:.3f}",
                showarrow=False, font=dict(color="white" if abs(r2_corr_matrix[i, j]) > 0.5 else "black"),
                row=2, col=2
            )
    
    # Summary statistics table (Plot 7)
    summary_data = []
    for tissue in tissues:
        stats = tissue_r2_stats[tissue]
        for fit_type in ['linear', 'loglog', 'semilog']:
            if stats[fit_type]:
                summary_data.append({
                    'tissue': tissue,
                    'fit_type': fit_type.title(),
                    'mean': stats[fit_type]['mean'],
                    'median': stats[fit_type]['median'],
                    'std': stats[fit_type]['std']
                })
    
    if summary_data:
        # Create grouped bar chart for summary statistics
        for i, stat in enumerate(['mean', 'median']):
            stat_values = []
            labels = []
            colors_stat = []
            
            for item in summary_data:
                stat_values.append(item[stat])
                labels.append(f"{item['tissue']}\n{item['fit_type']}")
                colors_stat.append(colors[hash(item['tissue']) % len(colors)])
            
            fig.add_trace(go.Bar(
                x=labels, y=stat_values, 
                name=f"{stat.title()}", 
                marker_color=colors_stat if stat == 'mean' else None,
                opacity=0.8 if stat == 'mean' else 0.5,
                hovertemplate=f"<b>%{{x}}</b><br>{stat.title()}: %{{y:.4f}}<extra></extra>"
            ), row=3, col=1)
    
    # Quality assessment (Plot 8)
    thresholds = [0.8, 0.9, 0.95]
    quality_data = {'threshold': [], 'proportion': [], 'fit_type': []}
    
    for fit_type, r2_values in [('Linear', all_r2_data['linear']), 
                                ('Log-Log', all_r2_data['loglog']), 
                                ('Semi-Log', all_r2_data['semilog'])]:
        if r2_values:
            r2_array = np.array(r2_values)
            for threshold in thresholds:
                proportion = np.mean(r2_array >= threshold)
                quality_data['threshold'].append(f"R² ≥ {threshold}")
                quality_data['proportion'].append(proportion)
                quality_data['fit_type'].append(fit_type)
    
    if quality_data['threshold']:
        fit_types = list(set(quality_data['fit_type']))
        for i, fit_type in enumerate(fit_types):
            mask = [ft == fit_type for ft in quality_data['fit_type']]
            thresholds_filtered = [quality_data['threshold'][j] for j in range(len(mask)) if mask[j]]
            proportions_filtered = [quality_data['proportion'][j] for j in range(len(mask)) if mask[j]]
            
            fig.add_trace(go.Bar(
                x=thresholds_filtered, y=proportions_filtered, 
                name=f"{fit_type} Quality", 
                marker_color=colors[i % len(colors)], opacity=0.7,
                hovertemplate=f"<b>%{{x}}</b><br>{fit_type}: %{{y:.3f}}<extra></extra>"
            ), row=3, col=2)
    
    # Improvement analysis (Plot 9)
    if len(all_r2_data['linear']) > 0 and len(all_r2_data['loglog']) > 0:
        linear_arr = np.array(all_r2_data['linear'])
        loglog_arr = np.array(all_r2_data['loglog'])
        
        if len(linear_arr) == len(loglog_arr):
            # Calculate improvement ratio
            improvement_ratio = loglog_arr / (linear_arr + 1e-10)  # Avoid division by zero
            
            fig.add_trace(go.Histogram(
                x=improvement_ratio, nbinsx=30, 
                name="Log-Log/Linear Ratio",
                marker_color='steelblue', opacity=0.7,
                hovertemplate="<b>Improvement Ratio</b><br>Ratio: %{x:.5f}<br>Count: %{y}<extra></extra>"
            ), row=3, col=3)
            
            # Add reference line at ratio = 1 (no improvement)
            fig.add_vline(x=1, line_dash="dash", line_color="red", 
                         annotation_text="No improvement", row=3, col=3)
    
    # Add diagonal reference lines to scatter plots
    for row, col in [(1, 1), (1, 2), (1, 3)]:
        fig.add_shape(type="line", x0=0, y0=0, x1=1, y1=1, 
                      line=dict(color="gray", dash="dash", width=2), 
                      row=row, col=col)
        
        # Add perfect fit reference lines
        fig.add_hline(y=0.95, line_dash="dot", line_color="green", 
                     annotation_text="Excellent (0.95)", row=row, col=col)
        fig.add_hline(y=0.8, line_dash="dot", line_color="orange", 
                     annotation_text="Good (0.8)", row=row, col=col)
        fig.add_vline(x=0.95, line_dash="dot", line_color="green", row=row, col=col)
        fig.add_vline(x=0.8, line_dash="dot", line_color="orange", row=row, col=col)
    
    # Enhanced layout
    fig.update_layout(
        title={
            'text': "📊 Comprehensive R² Goodness-of-Fit Analysis<br><sub>Detailed Comparison of Fitting Quality Across Linear, Log-Log, and Semi-Log Approaches</sub>",
            'x': 0.5,
            'xanchor': 'center',
            'font': {'size': 20}
        },
        width=1900, height=1500,
        showlegend=True,
        font=dict(family="Arial, sans-serif", size=10),
        plot_bgcolor='white',
        paper_bgcolor='white'
    )
    
    # Update axis labels with enhanced formatting
    axis_updates = [
        (1, 1, "R² Linear Scale", "R² Log-Log Scale"),
        (1, 2, "R² Linear Scale", "R² Semi-Log Scale"),
        (1, 3, "R² Log-Log Scale", "R² Semi-Log Scale"),
        (2, 1, "Fit Method", "R² Value"),
        (2, 3, "Power-Law Exponent (α)", "R² Value"),
        (3, 1, "Tissue & Fit Type", "R² Value"),
        (3, 2, "Quality Threshold", "Proportion"),
        (3, 3, "Log-Log / Linear Ratio", "Frequency")
    ]
    
    for row, col, xlabel, ylabel in axis_updates:
        fig.update_xaxes(title_text=xlabel, row=row, col=col)
        fig.update_yaxes(title_text=ylabel, row=row, col=col)
    
    return fig
    
    return fig

def _create_individual_statistical_test_plots(analyses: Dict[str, Dict], out_dir: str) -> int:
    """Create individual plots for each statistical test and comparison."""
    plot_count = 0
    tissues = list(analyses.keys())
    
    print(f"   🧪 Creating individual statistical test plots...")
    
    # 1. Mann-Whitney U tests for each pair of tissues
    print(f"      📊 Mann-Whitney U tests between tissues...")
    for i in range(len(tissues)):
        for j in range(i+1, len(tissues)):
            tissue1, tissue2 = tissues[i], tissues[j]
            blocks1 = analyses[tissue1]["blocks"]
            blocks2 = analyses[tissue2]["blocks"]
            
            alphas1 = np.array([b.alpha for b in blocks1 if np.isfinite(b.alpha)])
            alphas2 = np.array([b.alpha for b in blocks2 if np.isfinite(b.alpha)])
            
            if len(alphas1) > 0 and len(alphas2) > 0:
                # Perform Mann-Whitney U test
                try:
                    statistic, p_value = mannwhitneyu(alphas1, alphas2, alternative='two-sided')
                    
                    # Create detailed comparison plot
                    fig = make_subplots(
                        rows=2, cols=3,
                        subplot_titles=(
                            f"Alpha Distributions", f"Box Plot Comparison", f"Q-Q Plot",
                            f"Cumulative Distributions", f"Effect Size", f"Test Summary"
                        )
                    )
                    
                    # Plot 1: Overlapping histograms
                    fig.add_trace(go.Histogram(x=alphas1, name=tissue1, opacity=0.7, 
                                              marker_color='blue', nbinsx=20), row=1, col=1)
                    fig.add_trace(go.Histogram(x=alphas2, name=tissue2, opacity=0.7, 
                                              marker_color='red', nbinsx=20), row=1, col=1)
                    
                    # Plot 2: Box plots
                    fig.add_trace(go.Box(y=alphas1, name=tissue1, marker_color='blue'), row=1, col=2)
                    fig.add_trace(go.Box(y=alphas2, name=tissue2, marker_color='red'), row=1, col=2)
                    
                    # Plot 3: Q-Q plot
                    from scipy.stats import probplot
                    q1 = np.percentile(alphas1, np.linspace(0, 100, min(len(alphas1), 50)))
                    q2 = np.percentile(alphas2, np.linspace(0, 100, min(len(alphas2), 50)))
                    min_len = min(len(q1), len(q2))
                    fig.add_trace(go.Scatter(x=q1[:min_len], y=q2[:min_len], mode='markers',
                                            name='Q-Q', marker=dict(color='green', size=6)), row=1, col=3)
                    # Add diagonal line
                    min_val, max_val = min(np.min(q1), np.min(q2)), max(np.max(q1), np.max(q2))
                    fig.add_trace(go.Scatter(x=[min_val, max_val], y=[min_val, max_val],
                                            mode='lines', name='Perfect Match', 
                                            line=dict(color='gray', dash='dash')), row=1, col=3)
                    
                    # Plot 4: Cumulative distributions
                    alphas1_sorted = np.sort(alphas1)
                    alphas2_sorted = np.sort(alphas2)
                    cdf1 = np.arange(1, len(alphas1_sorted) + 1) / len(alphas1_sorted)
                    cdf2 = np.arange(1, len(alphas2_sorted) + 1) / len(alphas2_sorted)
                    
                    fig.add_trace(go.Scatter(x=alphas1_sorted, y=cdf1, mode='lines',
                                            name=f'{tissue1} CDF', line=dict(color='blue')), row=2, col=1)
                    fig.add_trace(go.Scatter(x=alphas2_sorted, y=cdf2, mode='lines',
                                            name=f'{tissue2} CDF', line=dict(color='red')), row=2, col=1)
                    
                    # Plot 5: Effect size (Cohen's d approximation)
                    mean1, mean2 = np.mean(alphas1), np.mean(alphas2)
                    std1, std2 = np.std(alphas1), np.std(alphas2)
                    pooled_std = np.sqrt(((len(alphas1)-1)*std1**2 + (len(alphas2)-1)*std2**2) / 
                                       (len(alphas1) + len(alphas2) - 2))
                    cohens_d = (mean1 - mean2) / pooled_std if pooled_std > 0 else 0
                    
                    fig.add_trace(go.Bar(x=['Effect Size'], y=[cohens_d], 
                                        marker_color='orange', name='Cohen\'s d'), row=2, col=2)
                    
                    # Plot 6: Test summary text
                    summary_text = f"""
                    Mann-Whitney U Test Results:
                    
                    Statistic: {statistic:.5f}
                    P-value: {_format_pvalue(p_value)}
                    Effect Size (Cohen's d): {cohens_d:.3f}
                    
                    {tissue1}: n={len(alphas1)}, μ={mean1:.3f}, σ={std1:.3f}
                    {tissue2}: n={len(alphas2)}, μ={mean2:.3f}, σ={std2:.3f}
                    
                    Interpretation:
                    {'Significant' if p_value < 0.05 else 'Not significant'} difference
                    Effect size: {'Small' if abs(cohens_d) < 0.5 else 'Medium' if abs(cohens_d) < 0.8 else 'Large'}
                    """
                    
                    fig.add_annotation(text=summary_text, xref="x domain", yref="y domain",
                                     x=0.5, y=0.5, showarrow=False, 
                                     font=dict(size=10, family="monospace"),
                                     bgcolor="lightgray", row=2, col=3)
                    
                    fig.update_layout(
                        title=f"Mann-Whitney U Test: {tissue1} vs {tissue2}",
                        width=1400, height=900,
                        showlegend=True
                    )
                    
                    # Save plot
                    plot_path = os.path.join(out_dir, f"TEST_MannWhitney_{tissue1}_vs_{tissue2}.html")
                    plot_path = plot_path.replace(" ", "_").replace("/", "_").replace("\\", "_")
                    fig.write_html(plot_path)
                    plot_count += 1
                    
                except Exception as e:
                    print(f"         Warning: Could not create Mann-Whitney test for {tissue1} vs {tissue2}: {e}")
    
    # 2. Kruskal-Wallis test for all tissues
    print(f"      🔬 Kruskal-Wallis test across all tissues...")
    all_alphas_by_tissue = []
    tissue_labels = []
    
    for tissue in tissues:
        blocks = analyses[tissue]["blocks"]
        alphas = [b.alpha for b in blocks if np.isfinite(b.alpha)]
        if len(alphas) > 0:
            all_alphas_by_tissue.append(alphas)
            tissue_labels.extend([tissue] * len(alphas))
    
    if len(all_alphas_by_tissue) > 2:
        try:
            kw_statistic, kw_p_value = kruskal(*all_alphas_by_tissue)
            
            # Create Kruskal-Wallis visualization
            fig = make_subplots(
                rows=2, cols=2,
                subplot_titles=(
                    "Box Plots by Tissue", "Violin Plots", 
                    "Mean Ranks", "Test Summary"
                )
            )
            
            colors = px.colors.qualitative.Set1
            
            # Plot 1: Box plots
            for i, (tissue, alphas) in enumerate(zip(tissues, all_alphas_by_tissue)):
                color = colors[i % len(colors)]
                fig.add_trace(go.Box(y=alphas, name=tissue, marker_color=color), row=1, col=1)
            
            # Plot 2: Violin plots
            for i, (tissue, alphas) in enumerate(zip(tissues, all_alphas_by_tissue)):
                color = colors[i % len(colors)]
                fig.add_trace(go.Violin(y=alphas, name=tissue, 
                                       marker_color=color, showlegend=False), row=1, col=2)
            
            # Plot 3: Mean ranks
            from scipy.stats import rankdata
            all_values = np.concatenate(all_alphas_by_tissue)
            ranks = rankdata(all_values)
            
            start_idx = 0
            mean_ranks = []
            for alphas in all_alphas_by_tissue:
                end_idx = start_idx + len(alphas)
                mean_rank = np.mean(ranks[start_idx:end_idx])
                mean_ranks.append(mean_rank)
                start_idx = end_idx
            
            fig.add_trace(go.Bar(x=tissues, y=mean_ranks, marker_color=colors[:len(tissues)],
                                name='Mean Ranks', showlegend=False), row=2, col=1)
            
            # Plot 4: Test summary
            kw_summary = f"""
            Kruskal-Wallis H Test Results:
            
            H Statistic: {kw_statistic:.3f}
            P-value: {_format_pvalue(kw_p_value)}
            Degrees of Freedom: {len(tissues)-1}
            
            Sample Sizes:
            {chr(10).join([f"{tissue}: n={len(alphas)}" for tissue, alphas in zip(tissues, all_alphas_by_tissue)])}
            
            Interpretation:
            {'Significant' if kw_p_value < 0.05 else 'Not significant'} difference between groups
            
            Effect Size (η²): {(kw_statistic - len(tissues) + 1) / (len(all_values) - len(tissues)):.3f}
            """
            
            fig.add_annotation(text=kw_summary, xref="x domain", yref="y domain",
                             x=0.5, y=0.5, showarrow=False,
                             font=dict(size=10, family="monospace"),
                             bgcolor="lightblue", row=2, col=2)
            
            fig.update_layout(
                title="Kruskal-Wallis Test: Alpha Differences Across All Tissues",
                width=1200, height=800,
                showlegend=True
            )
            
            kw_path = os.path.join(out_dir, "TEST_KruskalWallis_AllTissues.html")
            fig.write_html(kw_path)
            plot_count += 1
            
        except Exception as e:
            print(f"         Warning: Could not create Kruskal-Wallis test: {e}")
    
    # 3. Correlation tests between metrics
    print(f"      📈 Correlation analysis plots...")
    for tissue in tissues:
        blocks = analyses[tissue]["blocks"]
        
        # Extract multiple metrics
        alphas = np.array([b.alpha for b in blocks if np.isfinite(b.alpha)])
        ks_vals = np.array([b.ks_D_tail for b in blocks if np.isfinite(b.ks_D_tail)])
        r2_vals = np.array([b.R2_k_loglog for b in blocks if np.isfinite(b.R2_k_loglog)])
        p_vals = np.array([b.ks_pvalue for b in blocks if np.isfinite(b.ks_pvalue)])
        
        min_len = min(len(alphas), len(ks_vals), len(r2_vals), len(p_vals))
        if min_len > 5:  # Need reasonable sample size
            # Trim all arrays to same length
            alphas = alphas[:min_len]
            ks_vals = ks_vals[:min_len]
            r2_vals = r2_vals[:min_len]
            p_vals = p_vals[:min_len]
            
            # Create correlation matrix plot
            fig = make_subplots(
                rows=2, cols=3,
                subplot_titles=(
                    "Alpha vs KS-D", "Alpha vs R²", "Alpha vs P-value",
                    "KS-D vs R²", "Correlation Matrix", "Regression Summary"
                )
            )
            
            # Scatter plots with trend lines
            from scipy.stats import pearsonr, spearmanr
            
            # Alpha vs KS-D
            corr_alpha_ks, p_alpha_ks = pearsonr(alphas, ks_vals)
            fig.add_trace(go.Scatter(x=alphas, y=ks_vals, mode='markers',
                                    name=f'r={corr_alpha_ks:.3f}', marker=dict(size=6)), row=1, col=1)
            
            # Add trend line
            z = np.polyfit(alphas, ks_vals, 1)
            p_func = np.poly1d(z)
            fig.add_trace(go.Scatter(x=alphas, y=p_func(alphas), mode='lines',
                                    name='Trend', line=dict(color='red', dash='dash')), row=1, col=1)
            
            # Alpha vs R²
            corr_alpha_r2, p_alpha_r2 = pearsonr(alphas, r2_vals)
            fig.add_trace(go.Scatter(x=alphas, y=r2_vals, mode='markers',
                                    name=f'r={corr_alpha_r2:.3f}', marker=dict(size=6)), row=1, col=2)
            
            # Alpha vs P-value
            corr_alpha_p, p_alpha_p = pearsonr(alphas, p_vals)
            fig.add_trace(go.Scatter(x=alphas, y=p_vals, mode='markers',
                                    name=f'r={corr_alpha_p:.3f}', marker=dict(size=6)), row=1, col=3)
            
            # KS-D vs R²
            corr_ks_r2, p_ks_r2 = pearsonr(ks_vals, r2_vals)
            fig.add_trace(go.Scatter(x=ks_vals, y=r2_vals, mode='markers',
                                    name=f'r={corr_ks_r2:.3f}', marker=dict(size=6)), row=2, col=1)
            
            # Correlation matrix heatmap
            metrics = np.column_stack([alphas, ks_vals, r2_vals, p_vals])
            corr_matrix = np.corrcoef(metrics, rowvar=False)
            
            fig.add_trace(go.Heatmap(
                z=corr_matrix,
                x=['Alpha', 'KS-D', 'R²', 'P-value'],
                y=['Alpha', 'KS-D', 'R²', 'P-value'],
                colorscale='RdBu',
                zmid=0,
                showlegend=False
            ), row=2, col=2)
            
            # Summary statistics
            corr_summary = f"""
            Correlation Analysis: {tissue}
            
            Pearson Correlations:
            Alpha vs KS-D: r={corr_alpha_ks:.3f}, p={p_alpha_ks:.4f}
            Alpha vs R²: r={corr_alpha_r2:.3f}, p={p_alpha_r2:.4f}
            Alpha vs P-val: r={corr_alpha_p:.3f}, p={p_alpha_p:.4f}
            KS-D vs R²: r={corr_ks_r2:.3f}, p={p_ks_r2:.4f}
            
            Sample Size: n={min_len}
            
            Significant Correlations:
            {chr(10).join([f"- {pair}" for pair, p_val in [
                ("Alpha-KS", p_alpha_ks), ("Alpha-R²", p_alpha_r2), 
                ("Alpha-P", p_alpha_p), ("KS-R²", p_ks_r2)
            ] if p_val < 0.05])}
            """
            
            fig.add_annotation(text=corr_summary, xref="x domain", yref="y domain",
                             x=0.5, y=0.5, showarrow=False,
                             font=dict(size=9, family="monospace"),
                             bgcolor="lightyellow", row=2, col=3)
            
            fig.update_layout(
                title=f"Correlation Analysis: {tissue}",
                width=1500, height=900,
                showlegend=True
            )
            
            corr_path = os.path.join(out_dir, f"TEST_Correlations_{tissue}.html")
            corr_path = corr_path.replace(" ", "_").replace("/", "_").replace("\\", "_")
            fig.write_html(corr_path)
            plot_count += 1
    
    print(f"      ✅ Created {plot_count} individual statistical test plots")
    return plot_count

def make_powerlaw_gof_multi_reports_fast(
    csv_paths: Sequence[str],
    out_dir: str,
    pairs: Optional[Sequence[Tuple[str, str]]] = None,   # None => discover all shared prefixes (upper-tri incl. diagonal)
    max_rows: int = 400,
    max_cols: int = 400,
    seed: int = 42,
    # FAST knobs
    fast: bool = True,
    html_only_mode: bool = True,          # Skip PDF generation for speed
    n_jobs: int = 1,                      # parallel per-pair analysis
    max_cells_per_block: int = 200_000,   # downsample ceiling per block for stats/KS
    xmin_quantiles: Tuple[float,float,int] = (0.80, 0.96, 9),  # smaller grid is faster
    B_bootstrap: int = 60,                # fewer bootstraps
    min_tail: int = 150,                  # slightly higher skips tiny tails (fewer fits)
    file_display_names: Optional[Dict[str, str]] = None,
    write_per_file: bool = True,
    write_all_files: bool = True,
    write_comparisons_only: bool = True,
):
    """
    Faster build:
      - Per-file PDF+HTML (tables + integrated figures)
      - All-files PDF+HTML (tables only)
      - Comparisons-only PDF+HTML (rank tables)
    """
    print("🚀 Starting Power-Law GoF Multi-Reports Analysis (FAST)")
    print("=" * 60)
    print(f"📁 Output directory: {out_dir}")
    print(f"📊 Input files: {len(csv_paths)}")
    for i, path in enumerate(csv_paths, 1):
        print(f"   {i}. {os.path.basename(path)}")
    print(f"⚡ Configuration: n_jobs={n_jobs}, max_rows={max_rows}, max_cols={max_cols}")
    print(f"🔬 Analysis: B_bootstrap={B_bootstrap}, min_tail={min_tail}")
    print("=" * 60)
    
    os.makedirs(out_dir, exist_ok=True)
    random.seed(seed)
    
    print("🔍 Step 1/6: Discovering pairs and building canonical samples...")
    requested_pairs_global = _pairs_auto(csv_paths, pairs)
    print(f"   Found {len(requested_pairs_global)} unique pairs across files")
    
    canonical = _build_canonical_samples(csv_paths, requested_pairs_global, max_rows, max_cols, seed)
    print(f"   ✅ Canonical samples built successfully")

    analyses: Dict[str, Dict] = {}
    blocks_all: List[BlockGoF] = []

    print(f"\n📊 Step 2/6: Processing {len(csv_paths)} files...")
    for i, p in enumerate(csv_paths):
        file_num = i + 1
        label = _display_name_for_csv(p, file_display_names, baseline=(i==0))
        print(f"\n   📄 Processing file {file_num}/{len(csv_paths)}: {label}")
        print(f"      File path: {p}")
        
        # ------------ Read each pair ONCE and precompute v,k ------------
        print(f"      🔧 Reading and preparing canonical data...")
        pair_data: Dict[Tuple[str,str], Tuple[np.ndarray, np.ndarray]] = {}
        TS_vals = []; CT_vals = []; TS_k = []; CT_k = []
        # Fetch sequentially (IO-bound); analyze in parallel shortly
        pair_count = len(canonical.items())
        print(f"      Processing {pair_count} canonical pairs...")
        for pair_idx, ((a_n, b_n), (rows, cols)) in enumerate(canonical.items()):
            if (pair_idx + 1) % max(1, pair_count // 4) == 0 or pair_idx == 0 or pair_idx == pair_count - 1:
                print(f"         Reading pair {pair_idx + 1}/{pair_count}: {a_n} vs {b_n}")
            
            A = _fetch_subblock(p, rows, cols).to_numpy(float)
            v = A[np.isfinite(A)].ravel()
            v = _downsample_values(v, max_cells_per_block, seed ^ (hash((a_n,b_n)) & 0x7FFFFFFF))
            k = _degrees_from_block(A)
            k = _downsample_values(k, max(50_000, max_cells_per_block//4), seed ^ 0xABC123)
            pair_data[(a_n,b_n)] = (v, k)
            if a_n == b_n:
                TS_vals.append(v); TS_k.append(k)
            else:
                CT_vals.append(v); CT_k.append(k)

        print(f"      ✅ Data preparation complete. Processed {len(pair_data)} pairs")
        
        TSv  = _concat_finite(TS_vals);  TSk  = _concat_finite(TS_k)
        CTv  = _concat_finite(CT_vals);  CTk  = _concat_finite(CT_k)
        ALLv = _concat_finite([TSv, CTv]); ALLk = _concat_finite([TSk, CTk])

        print(f"      📈 Data summary - TS: {len(TSv):,} values, CT: {len(CTv):,} values, ALL: {len(ALLv):,} values")

        # ------------ Analyze blocks in parallel ------------
        print(f"      🚀 Starting statistical analysis...")
        blocks: List[BlockGoF] = []

        tasks = [((a_n,b_n), pair_data[(a_n,b_n)], xmin_quantiles, B_bootstrap, min_tail, fast, seed + 1000*i + idx)
                 for idx, (a_n,b_n) in enumerate(pair_data.keys())]
        
        print(f"      Created {len(tasks)} analysis tasks")

        if n_jobs and n_jobs > 1:
            print(f"      🔧 Running parallel analysis with {n_jobs} workers...")
            with ProcessPoolExecutor(max_workers=n_jobs) as ex:
                futures = [ex.submit(_analyze_pair_task, t) for t in tasks]
                completed = 0
                for fut in as_completed(futures):
                    completed += 1
                    a_n, b_n, ana = fut.result()
                    blocks.append(BlockGoF(
                        file_label=label, file_path=p, region_a=a_n, region_b=b_n,
                        n_finite=ana["n_finite"], tail_n=ana["tail_n"],
                        alpha=ana["alpha"], xmin=ana["xmin"], ks_D_tail=ana["ks_D_tail"], p_plaus=ana["p_plaus"],
                        ll_pl=ana["ll_pl"], aic_pl=ana["aic_pl"], ll_exp=ana["ll_exp"], aic_exp=ana["aic_exp"],
                        ll_logn=ana["ll_logn"], aic_logn=ana["aic_logn"], dAIC_exp=ana["dAIC_exp"], dAIC_logn=ana["dAIC_logn"],
                        aicw_pl=ana["aicw_pl"], aicw_exp=ana["aicw_exp"], aicw_logn=ana["aicw_logn"],
                        R2_k_linear=ana["R2_k_linear"], R2_k_loglog=ana["R2_k_loglog"], R2_k_semilog=ana["R2_k_semilog"],
                        ks_pvalue=_ks_pvalue_from_statistic(ana["tail_n"], ana["ks_D_tail"]),
                        v_data=ana["v_data"], k_data=ana["k_data"]
                    ))
                    if completed % max(1, len(tasks) // 10) == 0:  # Progress every 10%
                        print(f"         Progress: {completed}/{len(tasks)} pairs completed ({100*completed/len(tasks):.1f}%)")
        else:
            print(f"      🔧 Running sequential analysis...")
            for idx, t in enumerate(tasks):
                a_n, b_n, ana = _analyze_pair_task(t)
                blocks.append(BlockGoF(
                    file_label=label, file_path=p, region_a=a_n, region_b=b_n,
                    n_finite=ana["n_finite"], tail_n=ana["tail_n"],
                    alpha=ana["alpha"], xmin=ana["xmin"], ks_D_tail=ana["ks_D_tail"], p_plaus=ana["p_plaus"],
                    ll_pl=ana["ll_pl"], aic_pl=ana["aic_pl"], ll_exp=ana["ll_exp"], aic_exp=ana["aic_exp"],
                    ll_logn=ana["ll_logn"], aic_logn=ana["aic_logn"], dAIC_exp=ana["dAIC_exp"], dAIC_logn=ana["dAIC_logn"],
                    aicw_pl=ana["aicw_pl"], aicw_exp=ana["aicw_exp"], aicw_logn=ana["aicw_logn"],
                    R2_k_linear=ana["R2_k_linear"], R2_k_loglog=ana["R2_k_loglog"], R2_k_semilog=ana["R2_k_semilog"],
                    ks_pvalue=_ks_pvalue_from_statistic(ana["tail_n"], ana["ks_D_tail"]),
                    v_data=ana["v_data"], k_data=ana["k_data"]
                ))
                if (idx + 1) % max(1, len(tasks) // 10) == 0:  # Progress every 10%
                    print(f"         Progress: {idx + 1}/{len(tasks)} pairs completed ({100*(idx+1)/len(tasks):.1f}%)")

        print(f"      ✅ Statistical analysis complete. Generated {len(blocks)} blocks")

        # ------------ Aggregated TS/CT/ALL analysis using same fast block fit ------------
        print(f"      📊 Generating aggregated TS/CT/ALL analysis...")
        def _agg(label_kind, vals, kvals):
            ana = _fit_block_fast(vals, kvals, xmin_quantiles, B_bootstrap, min_tail, fast_bootstrap=fast, seed=seed+1234)
            return AggGoF(
                file_label=label, kind=label_kind,
                n_finite=ana["n_finite"], tail_n=ana["tail_n"],
                alpha=ana["alpha"], xmin=ana["xmin"], ks_D_tail=ana["ks_D_tail"], 
                ks_pvalue=_ks_pvalue_from_statistic(ana["tail_n"], ana["ks_D_tail"]),
                p_plaus=ana["p_plaus"],
                ll_pl=ana["ll_pl"], aic_pl=ana["aic_pl"], ll_exp=ana["ll_exp"], aic_exp=ana["aic_exp"],
                ll_logn=ana["ll_logn"], aic_logn=ana["aic_logn"], dAIC_exp=ana["dAIC_exp"], dAIC_logn=ana["dAIC_logn"],
                aicw_pl=ana["aicw_pl"], aicw_exp=ana["aicw_exp"], aicw_logn=ana["aicw_logn"],
                R2_k_linear=ana["R2_k_linear"], R2_k_loglog=ana["R2_k_loglog"], R2_k_semilog=ana["R2_k_semilog"],
            )
        aggs = [_agg("TS", TSv, TSk), _agg("CT", CTv, CTk), _agg("ALL", ALLv, ALLk)]
        print(f"      ✅ Aggregated analysis complete")

        # integrated figures
        print(f"      🎨 Generating visualization figures...")
        figs = {
            "TS":  _subplots_block(TSv,  TSk,  f"{label} — Integrated TS"),
            "CT":  _subplots_block(CTv,  CTk,  f"{label} — Integrated CT"),
            "ALL": _subplots_block(ALLv, ALLk, f"{label} — Integrated TS+CT"),
        }
        print(f"      ✅ Generated {len(figs)} visualization figures")

        analyses[label] = {"blocks": blocks, "aggs": aggs, "figs": figs, "path": p}
        blocks_all.extend(blocks)
        print(f"   ✅ File {file_num}/{len(csv_paths)} processing complete: {label}")

    # ---------- Write reports ----------
    print(f"\n📝 Step 3/6: Generating reports...")
    report_count = 0
    
    if write_per_file:
        print(f"   📄 Generating per-file reports...")
        for idx, (label, content) in enumerate(analyses.items()):
            print(f"      Writing report {idx + 1}/{len(analyses)}: {label}")
            print(f"         🔄 Blocks: {len(content['blocks'])}, Figures: {len(content['figs'])}")
            
            # Generate paths
            pdf_path  = os.path.join(out_dir, f"GoF_FAST_{label}.pdf").replace("/","_").replace("\\","_")
            html_path = os.path.join(out_dir, f"GoF_FAST_{label}.html").replace("/","_").replace("\\","_")
            
            # Skip PDF generation in HTML-only mode
            if html_only_mode:
                pdf_path = None
                print(f"         ⚡ HTML-only mode: skipping PDF generation")
            
            # Add timing for debugging
            import time
            start_time = time.time()
            _render_per_file_report(label, content["blocks"], content["aggs"], content["figs"], pdf_path, html_path)
            elapsed = time.time() - start_time
            print(f"         ✅ Report generated in {elapsed:.5f}s")
            report_count += 1 if html_only_mode else 2
        print(f"   ✅ Per-file reports complete")

    if write_all_files:
        print(f"   📋 Generating combined files report...")
        pdf_all  = os.path.join(out_dir, "GoF_FAST_ALL_FILES.pdf")
        html_all = os.path.join(out_dir, "GoF_FAST_ALL_FILES.html")
        _render_all_files_report(analyses, pdf_all, html_all)
        report_count += 2
        print(f"   ✅ Combined files report complete")

    if write_comparisons_only:
        print(f"   📊 Generating comparisons-only report...")
        pdf_cmp  = os.path.join(out_dir, "GoF_FAST_COMPARISONS_ONLY.pdf")
        html_cmp = os.path.join(out_dir, "GoF_FAST_COMPARISONS_ONLY.html")
        _render_comparisons_only(blocks_all, csv_paths, file_display_names, pdf_cmp, html_cmp, analyses)
        report_count += 2
        print(f"   ✅ Comparisons-only report complete")

    # ---------- Generate Combined Reports (New Feature) ----------
    print(f"\n🔄 Step 4/6: Generating enhanced combined visualizations...")
    # Extract tissue types and create combined visualizations
    all_tissues_gofs = []
    tissue_names = []
    ts_gofs = []
    ts_tissues = []
    ct_gofs = []
    ct_tissues = []
    
    print(f"   📊 Processing tissues for combined analysis...")
    # Group data by tissues/regions
    for label, content in analyses.items():
        blocks = content["blocks"]
        # Get first block to represent this tissue/file
        if blocks:
            representative_block = blocks[0]  # Use first block as representative
            all_tissues_gofs.append(representative_block)
            tissue_names.append(label)
            
            # Categorize as TS or CT based on label
            if "TS" in label.upper():
                ts_gofs.append(representative_block)
                ts_tissues.append(label)
            elif "CT" in label.upper():
                ct_gofs.append(representative_block)
                ct_tissues.append(label)
    
    print(f"      Found {len(all_tissues_gofs)} tissues total ({len(ts_gofs)} TS, {len(ct_gofs)} CT)")
    
    # Generate combined tissues report
    if len(all_tissues_gofs) > 1:
        print(f"   🎨 Generating combined tissues visualization...")
        combined_fig = _render_combined_tissues_report(all_tissues_gofs, tissue_names, "All_Tissues")
        combined_html_path = os.path.join(out_dir, "GoF_FAST_COMBINED_TISSUES.html")
        combined_fig.write_html(combined_html_path)
        print(f"      ✅ Combined tissues report saved: {os.path.basename(combined_html_path)}")
        report_count += 1
    
    # Generate TS vs CT comparison if we have both types
    if len(ts_gofs) > 0 and len(ct_gofs) > 0:
        print(f"   ⚡ Generating TS vs CT comparison...")
        ts_ct_fig = _render_ts_vs_ct_comparison(ts_gofs, ct_gofs, ts_tissues, ct_tissues)
        ts_ct_html_path = os.path.join(out_dir, "GoF_FAST_TS_vs_CT_COMPARISON.html")
        ts_ct_fig.write_html(ts_ct_html_path)
        print(f"      ✅ TS vs CT comparison report saved: {os.path.basename(ts_ct_html_path)}")
        report_count += 1

    # ---------- Generate Enhanced Analysis Reports (New Feature) ----------
    print(f"\n🔬 Step 4.5/6: Generating enhanced analysis visualizations...")
    
    # Enhanced tissue comparison with comprehensive metrics
    print(f"   📈 Creating enhanced tissue comparison...")
    enhanced_tissue_fig = _create_enhanced_tissue_comparison(analyses)
    enhanced_tissue_path = os.path.join(out_dir, "GoF_FAST_ENHANCED_TISSUE_COMPARISON.html")
    enhanced_tissue_fig.write_html(enhanced_tissue_path)
    print(f"      ✅ Enhanced tissue comparison saved: {os.path.basename(enhanced_tissue_path)}")
    report_count += 1
    
    # KS statistic visualization
    print(f"   🎯 Creating KS statistic visualizations...")
    ks_viz_fig = _create_ks_visualization_report(analyses)
    ks_viz_path = os.path.join(out_dir, "GoF_FAST_KS_VISUALIZATION.html")
    ks_viz_fig.write_html(ks_viz_path)
    print(f"      ✅ KS visualization report saved: {os.path.basename(ks_viz_path)}")
    report_count += 1
    
    # R² analysis with log10 scale (linear and logarithmic binning)
    print(f"   📊 Creating R² log scale analysis (linear & logarithmic binning)...")
    r2_log_fig = _create_r2_log_scale_analysis(analyses)
    r2_log_path = os.path.join(out_dir, "GoF_FAST_R2_LOG_SCALE_ANALYSIS.html")
    r2_log_fig.write_html(r2_log_path)
    print(f"      ✅ R² log scale analysis saved: {os.path.basename(r2_log_path)}")
    report_count += 1
    
    # P-value analysis with scientific notation
    print(f"   🔬 Creating p-value scientific notation analysis...")
    pvalue_fig = _create_pvalue_scientific_notation_plot(analyses)
    pvalue_path = os.path.join(out_dir, "GoF_FAST_PVALUE_SCIENTIFIC.html")
    pvalue_fig.write_html(pvalue_path)
    print(f"      ✅ P-value scientific notation report saved: {os.path.basename(pvalue_path)}")
    report_count += 1
    
    # Individual KS visualizations for each tissue/pair combination
    print(f"   🔍 Creating individual KS visualizations...")
    ks_individual_count = 0
    for tissue_label, content in analyses.items():
        blocks = content["blocks"]
        for i, block in enumerate(blocks[:3]):  # Limit to first 3 blocks per tissue
            if hasattr(block, 'v_data') and len(block.v_data) > 0:
                pair_label = f"{block.region_a}_vs_{block.region_b}"
                ks_fig = _visualize_ks_statistic(
                    np.array(block.v_data), 
                    block.xmin, 
                    block.alpha, 
                    f"KS Analysis: {tissue_label} - {pair_label}"
                )
                ks_individual_path = os.path.join(out_dir, f"GoF_FAST_KS_{tissue_label}_{pair_label}_{i}.html")
                ks_individual_path = ks_individual_path.replace(" ", "_").replace("/", "_").replace("\\", "_")
                ks_fig.write_html(ks_individual_path)
                ks_individual_count += 1
        
    print(f"      ✅ Created {ks_individual_count} individual KS visualizations")
    report_count += ks_individual_count
    
    # Individual statistical test plots
    print(f"   🧪 Creating comprehensive individual test plots...")
    test_plot_count = _create_individual_statistical_test_plots(analyses, out_dir)
    report_count += test_plot_count

    # Barabási Network Science Comprehensive Goodness-of-Fit Analysis
    print(f"   🔬 Creating Barabási Network Science comprehensive analysis...")
    try:
        # Convert analyses to the format expected by Barabási function
        combined_results_barabasi = {}
        for tissue_label, content in analyses.items():
            tissue_results = []
            blocks = content.get("blocks", [])
            
            # Create file_result format: [filename, tissue_type, pair_info, block_data]
            for block in blocks:
                if hasattr(block, 'combined_values') and len(block.combined_values) > 0:
                    file_result = [
                        f"synthetic_{tissue_label}",  # filename
                        tissue_label,  # tissue_type
                        f"{block.region_a}_vs_{block.region_b}",  # pair_info
                        [block]  # block_data
                    ]
                    tissue_results.append(file_result)
            
            if tissue_results:
                combined_results_barabasi[tissue_label] = tissue_results
        
        if combined_results_barabasi:
            barabasi_dir = _create_barabasi_goodness_of_fit_report(
                combined_results_barabasi, out_dir, "quantile"
            )
            if barabasi_dir:
                barabasi_files = list(Path(barabasi_dir).glob("*.html"))
                report_count += len(barabasi_files)
                print(f"      ✅ Created {len(barabasi_files)} Barabási analysis reports")
            else:
                print(f"      ⚠ Barabási analysis not generated (insufficient data or module unavailable)")
        else:
            print(f"      ⚠ No data available for Barabási analysis")
            
    except Exception as e:
        print(f"      ❌ Error in Barabási analysis: {e}")
        if BARABASI_AVAILABLE:
            import traceback
            print(f"      📋 Traceback: {traceback.format_exc()}")

    # ---------- Generate HTML Navigation Index ----------
    print(f"\n🌐 Step 5/7: Creating HTML navigation index...")
    _create_html_index(out_dir, analyses, report_count)
    print(f"   ✅ HTML navigation index created")

    # ---------- Final Summary ----------
    print(f"\n🎉 Step 6/7: Analysis Complete!")
    print(f"   📊 Files processed: {len(csv_paths)}")
    print(f"   🔬 Total blocks analyzed: {len(blocks_all):,}")
    print(f"   📄 Reports generated: {report_count}")
    print(f"   🌐 Navigation: Open index.html in {out_dir}")
    print(f"   ✅ All outputs saved to: {out_dir}")

    return dict(canonical=canonical, analyses=analyses)

# -------------------- Example usage --------------------
if __name__ == "__main__":
    # Replace with your paths
    CSV_CT2_TS3 = "/media/psylab-6028/DATA/Eden/CoExpression_ReProduction/nbs/xwgcna_rosmap_constBeta_CT2_TS3_adjacency.csv"
    CSV_CT1_TS1 = "/media/psylab-6028/DATA/Eden/CoExpression_ReProduction/nbs/xwgcna_rosmap_constBeta_CT1_TS1_adjacency.csv"
    CSV_Dynamic = "/media/psylab-6028/DATA/Eden/CoExpression_ReProduction/notebooks/WGCNA_new_ROSMAP_AdjFunc_adjacency_no_NEWFUNC_ISONEW_zerDiag2025-09-08.csv"
    out_dir = "/media/psylab-6028/DATA/Eden/CoExpression_ReProduction/outputsValidationPlots/PL_GoF_MULTI_FAST"
    os.makedirs(out_dir, exist_ok=True)

    file_names = {
        os.path.basename(CSV_CT1_TS1): "CT1_TS1 (baseline)",
        os.path.basename(CSV_CT2_TS3): "CT2_TS3",
        os.path.basename(CSV_Dynamic): "Dynamic Beta"
    }

    # Dynamic n_jobs calculation based on current system state
    import psutil
    available_ram_gb = psutil.virtual_memory().available / (1024**3)
    current_load = os.getloadavg()[0] if hasattr(os, 'getloadavg') else 0
    
    # Base recommendation: 30 for this system
    base_n_jobs = 30
    
    # Adjust for current load
    if current_load > 16:  # High load (50% of cores)
        dynamic_n_jobs = max(10, base_n_jobs - 10)
    elif available_ram_gb < 50:  # Low memory
        dynamic_n_jobs = max(10, base_n_jobs - 5)
    else:
        dynamic_n_jobs = base_n_jobs
    
    print(f"Dynamic n_jobs calculation: {dynamic_n_jobs} (RAM: {available_ram_gb:.1f}GB, Load: {current_load:.5f})")

    make_powerlaw_gof_multi_reports_fast(
        csv_paths=[CSV_CT1_TS1, CSV_CT2_TS3, CSV_Dynamic],
        out_dir=out_dir,
        pairs=None,                # discover common prefixes, upper-tri incl. diagonal
        max_rows=500,              # you can lower further
        max_cols=500,
        seed=123,
        fast=True,
        n_jobs=dynamic_n_jobs,     # Dynamic adjustment
        max_cells_per_block=250000,  # you can lower further
        xmin_quantiles=(0.82, 0.95, 7),
        B_bootstrap=100,
        min_tail=100,
        file_display_names=file_names,
        write_per_file=True,
        write_all_files=True,
        write_comparisons_only=True,
    )
