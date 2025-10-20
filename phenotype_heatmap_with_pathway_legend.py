"""
Phenotype × Module heatmap with right‑side pathway annotations & legend
======================================================================

What this script does
---------------------
1) Loads three inputs (your filenames may vary):
   • PHENO_TSV: long table of module–phenotype stats (must include columns for module id, phenotype name, p_adj/qvalue, and optionally tissue)
   • KEGG_CSV: enrichment results per module (must include columns for Cluster (module id), Description (pathway), category, subcategory, and p/q-values)
   • DETAILS_TSV: module details (must include Cluster ID, and any other descriptors you want in the report)

2) Builds a heatmap where:
   • Y‑axis = modules (Cluster IDs)
   • X‑axis = phenotypes
   • Cell value = −log10(p_adj) by default (configurable)
   • Left side: two color strips showing category and subcategory annotations per module (from KEGG enrichment)
   • Right side: a two‑column legend of categories and subcategories (with colors) inferred from KEGG_CSV
   • Output: PDF file

3) Filtering options:
   • Select a subset of tissues via a 'tissue' column (if present)
   • Include/exclude phenotypes by name (substring or regex)
   • Select top N modules by phenotype significance (top_modules parameter)

4) Tabular outputs:
   • module_details_selected.tsv: module rows from DETAILS_TSV for the modules seen in the heatmap (or for a specified subset)
   • top_phenotypes_per_module.tsv: for each module, the top‑K phenotypes (lowest p_adj)
   • top_pathways_per_module.tsv: for each module, the top‑K pathways (lowest qvalue, fallback to p.adjust/pvalue)

Usage inside a Jupyter notebook
-------------------------------
from phenotype_heatmap_with_pathway_legend import run_all

run_all(
    pheno_tsv="/mnt/data/rosmap_ME_pheno_ME_vs_pheno_correlations.tsv",
    kegg_csv="/mnt/data/kegg_rosmap_constBeta_CT2_TS3.csv",
    details_tsv="/mnt/data/xwgcna_rosmap_constBeta_CT2_TS3_Cluster_details.tsv",
    out_prefix="rosmap_constBeta_CT2_TS3",
    tissues=None,                      # e.g., ["AC","MFBA9BA46","PCGBA23"] or None
    include_phenotypes=None,           # e.g., [r"age", r"cognitive|memory"]
    exclude_phenotypes=None,           # e.g., [r"sex|gender"]
    metric="corr",                     # one of {"neglog10_padj","neglog10_q","neglog10_p","corr","beta"}
    cap=1.0,
    top_modules=None,                  # e.g., 100 to select top 100 modules by phenotype significance
    top_k_phenos=10,
    top_k_pathways=5,
    figure_dpi=180,
    figsize_scale=1.0,                 # Increase if labels are crowded
    cluster_by="category_then_subcategory",  # Options: "none", "category", "subcategory", "category_then_subcategory"
)

Notes
-----
• The code tries to auto‑detect column names with flexible mappers, so slight header differences are OK.
• If padj/qvalue are missing, the code falls back to pvalue; if phenotypes lack p_adj entirely, you can use metric="corr" or "beta".
• Colors for category/subcategory are auto‑assigned from matplotlib tab colormaps.
"""

from __future__ import annotations
import re
import json
import math
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.gridspec import GridSpec

# -----------------------------
# Column mappers (robust names)
# -----------------------------

def _lower_map(cols: Sequence[str]) -> Dict[str, str]:
    return {c.lower().strip(): c for c in cols}

# Phenotype table columns
# Expected long format: one row per (module, phenotype[, tissue]) with stats

def _map_pheno_cols(df: pd.DataFrame) -> Dict[str, Optional[str]]:
    L = _lower_map(df.columns)
    def pick(*keys):
        for k in keys:
            if k in L: return L[k]
        return None
    return {
        "module":   pick("cluster id","cluster","module","module_id","module id"),
        "phenotype":pick("phenotype","trait","pheno","pheno_name","name"),
        "tissue":   pick("tissue","tissues"),
        "padj":     pick("p_adj","padj","p.adjust","fdr","adj p","adj_p","adj.p"),
        "qvalue":   pick("qvalue","q","q.value"),
        "pvalue":   pick("pvalue","p","p_val","p-val"),
        "corr":     pick("cor","corr","rho","r","pearson_r","spearman_r"),
        "beta":     pick("beta","effect","estimate"),
    }

# KEGG enrichment columns

def _map_kegg_cols(df: pd.DataFrame) -> Dict[str, Optional[str]]:
    L = _lower_map(df.columns)
    def pick(*keys):
        for k in keys:
            if k in L: return L[k]
        return None
    return {
        "cluster":     pick("cluster","cluster id","module","module id"),
        "desc":        pick("description","term","pathway","name"),
        "category":    pick("category","catagory"),
        "subcategory": pick("subcategory","subcatagory","sub_category","sub-cat"),
        "p":           pick("pvalue","p","p_val"),
        "padj":        pick("p.adjust","padj","fdr"),
        "qvalue":      pick("qvalue","q","q.value"),
    }

# Details table columns

def _map_details_cols(df: pd.DataFrame) -> Dict[str, Optional[str]]:
    L = _lower_map(df.columns)
    def pick(*keys):
        for k in keys:
            if k in L: return L[k]
        return None
    return {
        "cluster_id": pick("cluster id","cluster","module","module id"),
        "cluster_type": pick("cluster type","type"),
        "cluster_size": pick("cluster size","size","gene count","genes","n_genes"),
        "tissue": pick("tissue","dominant tissue","main tissue"),
    }

# -----------------------------
# IO helpers
# -----------------------------

def _read_table(path: str) -> pd.DataFrame:
    if path.lower().endswith(".csv"):
        return pd.read_csv(path, engine="python", on_bad_lines="skip")
    else:
        return pd.read_csv(path, sep="\t", engine="python", on_bad_lines="skip")

# -----------------------------
# Filtering utilities
# -----------------------------

def _regex_any(patterns: Optional[Sequence[str]]) -> Optional[re.Pattern]:
    if not patterns:
        return None
    pat = "|".join(f"({p})" for p in patterns)
    try:
        return re.compile(pat, flags=re.IGNORECASE)
    except re.error:
        # Fallback: escape everything
        pat = "|".join(re.escape(p) for p in patterns)
        return re.compile(pat, flags=re.IGNORECASE)

# -----------------------------
# Build phenotype matrix
# -----------------------------

def build_pheno_matrix(
    pheno_tsv: str,
    tissues: Optional[Sequence[str]] = None,
    include_phenotypes: Optional[Sequence[str]] = None,
    exclude_phenotypes: Optional[Sequence[str]] = None,
    metric: str = "corr",
    cap: float = 1.0,
    top_modules: Optional[int] = None,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Returns (matrix_df, padj_matrix, pheno_long_filtered)
    matrix_df: rows=modules (int), cols=phenotypes (str), values per `metric` (signed correlations)
    padj_matrix: rows=modules, cols=phenotypes, values=adjusted p-values for significance stars
    pheno_long_filtered: the filtered long form table (for downstream tables)
    
    Parameters
    ----------
    metric : str
        Default changed to "corr" for signed correlation values
    top_modules : int, optional
        If provided, select only the top N modules ranked by minimum p-value across all phenotypes
    """
    # Read and validate input
    df_raw = _read_table(pheno_tsv)
    M = _map_pheno_cols(df_raw)
    
    need = ["module","phenotype"]
    for k in need:
        if M[k] is None:
            raise ValueError(f"Missing required column in phenotype file: {k}")

    # Normalize basic columns
    df = df_raw.copy()
    
    # Extract numeric part from module column (handles "M15", "15", etc.)
    module_series = df[M["module"]].astype(str).str.strip()
    # Try to extract numbers, handling formats like "M15" or just "15"
    module_numeric = module_series.str.extract(r'(\d+)', expand=False)
    df[M["module"]] = pd.to_numeric(module_numeric, errors="coerce").astype("Int64")
    
    df = df[df[M["module"]].notna()].copy()
    
    if df.empty:
        raise ValueError(f"No valid module IDs found. Sample module values: {df_raw[M['module']].head(10).tolist()}")
    
    df[M["module"]] = df[M["module"]].astype(int)
    df[M["phenotype"]] = df[M["phenotype"]].astype(str).str.strip()

    # Optional tissue filter
    if tissues and M["tissue"] is not None:
        keep = {str(t).strip() for t in tissues if str(t).strip()}
        if keep:
            filtered = df[df[M["tissue"]].astype(str).str.strip().isin(keep)].copy()
            if not filtered.empty:
                df = filtered

    # Include/exclude phenotypes
    inc_pat = _regex_any(include_phenotypes)
    exc_pat = _regex_any(exclude_phenotypes)
    if inc_pat:
        filtered = df[df[M["phenotype"]].str.contains(inc_pat, na=False, regex=True)].copy()
        if not filtered.empty:
            df = filtered
    if exc_pat:
        df = df[~df[M["phenotype"]].str.contains(exc_pat, na=False, regex=True)].copy()

    # Determine measure (now supporting signed correlations)
    vals: Optional[pd.Series] = None
    if metric == "corr" and M["corr"] is not None:
        vals = pd.to_numeric(df[M["corr"]], errors="coerce")  # Keep sign for diverging colormap
    elif metric == "beta" and M["beta"] is not None:
        vals = pd.to_numeric(df[M["beta"]], errors="coerce")  # Keep sign
    elif metric == "neglog10_padj" and M["padj"] is not None:
        s = pd.to_numeric(df[M["padj"]], errors="coerce")
        vals = -np.log10(s.clip(lower=np.nextafter(0, 1)))
    elif metric == "neglog10_q" and M["qvalue"] is not None:
        s = pd.to_numeric(df[M["qvalue"]], errors="coerce")
        vals = -np.log10(s.clip(lower=np.nextafter(0, 1)))
    elif metric == "neglog10_p" and M["pvalue"] is not None:
        s = pd.to_numeric(df[M["pvalue"]], errors="coerce")
        vals = -np.log10(s.clip(lower=np.nextafter(0, 1)))
    else:
        # Fallback priority
        for key in ("corr","beta","padj","qvalue","pvalue"):
            if M[key] is not None:
                s = pd.to_numeric(df[M[key]], errors="coerce")
                if key in {"padj","qvalue","pvalue"}:
                    vals = -np.log10(s.clip(lower=np.nextafter(0, 1)))
                else:
                    vals = s  # Keep sign for corr/beta
                metric = f"auto_from_{key}"
                break
    if vals is None:
        raise ValueError("Could not find a numeric measure to plot. Please set `metric` to one of the available columns.")
    
    # Get p-adjusted values for significance stars
    padj_vals = pd.to_numeric(df[M["padj"]], errors="coerce") if M["padj"] else None

    # Build main value matrix
    mat_data = pd.DataFrame({
        "module": df[M["module"]],
        "phenotype": df[M["phenotype"]],
        "val": vals
    }).dropna(subset=["val"])  # Drop NaN values to avoid KeyError
    
    # For signed metrics (corr/beta), take the value with max absolute value
    if metric in ("corr", "beta") or "corr" in metric or "beta" in metric:
        def get_max_abs(g):
            if g.empty or g["val"].isna().all():
                return pd.Series({"val": np.nan})
            idx = g["val"].abs().idxmax()
            return g.loc[[idx]].iloc[0]
        
        mat = (mat_data.groupby(["module","phenotype"], as_index=False)
               .apply(get_max_abs, include_groups=False)
               .pivot(index="module", columns="phenotype", values="val")
               .sort_index())
    else:
        mat = (mat_data.groupby(["module","phenotype"], as_index=False)
               .agg(val=("val","max"))
               .pivot(index="module", columns="phenotype", values="val")
               .sort_index())
    
    # Build p-adjusted matrix for significance stars
    padj_mat = None
    if padj_vals is not None:
        padj_data = pd.DataFrame({
            "module": df[M["module"]],
            "phenotype": df[M["phenotype"]],
            "padj": padj_vals
        })
        padj_mat = (padj_data.groupby(["module","phenotype"], as_index=False)
                    .agg(padj=("padj","min"))
                    .pivot(index="module", columns="phenotype", values="padj")
                    .reindex(index=mat.index, columns=mat.columns))
    
    # Select top modules by minimum p-value if requested
    if top_modules is not None and len(mat) > top_modules and padj_mat is not None:
        # Compute minimum p-value per module across all phenotypes
        module_scores = padj_mat.min(axis=1)
        top_mods = module_scores.nsmallest(top_modules).index
        mat = mat.loc[top_mods]
        if padj_mat is not None:
            padj_mat = padj_mat.loc[top_mods]
        # Filter df to match selected modules
        df = df[df[M["module"]].isin(top_mods)].copy()
    
    if mat.size > 0 and np.isfinite(cap):
        mat = mat.clip(lower=-cap, upper=cap)

    return mat, padj_mat, df

# -----------------------------
# Build category/subcategory palettes and legend
# -----------------------------

def _uniq_clean(series: pd.Series) -> List[str]:
    ss = series.astype(str).str.strip()
    uniq = [x for x in pd.unique(ss) if x and x.lower() not in ("nan","none")]
    return uniq


def build_cat_palettes(kegg_csv: str) -> Tuple[Dict[str, Tuple[float,float,float,float]], Dict[str, Tuple[float,float,float,float]]]:
    df = _read_table(kegg_csv)
    M = _map_kegg_cols(df)
    if M["category"] is None and M["subcategory"] is None:
        return {}, {}
    cats = _uniq_clean(df[M["category"]]) if M["category"] else []
    subs = _uniq_clean(df[M["subcategory"]]) if M["subcategory"] else []

    # Use matplotlib tab colormaps (tab20 + tab20b + tab20c + Set3) for variety
    pool = list(plt.cm.tab20.colors) + list(plt.cm.tab20b.colors) + list(plt.cm.tab20c.colors) + list(plt.cm.Set3.colors)
    cat_pal = {lab: pool[i % len(pool)] for i, lab in enumerate(cats)}
    sub_pal = {lab: pool[i % len(pool)] for i, lab in enumerate(subs)}
    return cat_pal, sub_pal


def get_module_annotations(kegg_csv: str, modules: Sequence[int]) -> pd.DataFrame:
    """
    Extract dominant category and subcategory for each module from KEGG enrichment data.
    Returns a DataFrame with index=module_id and columns=[category, subcategory]
    """
    df = _read_table(kegg_csv)
    M = _map_kegg_cols(df)
    
    df = df.copy()
    df[M["cluster"]] = pd.to_numeric(df[M["cluster"]], errors="coerce").astype("Int64")
    df = df[df[M["cluster"]].notna()].copy()
    df[M["cluster"]] = df[M["cluster"]].astype(int)
    
    # Filter to requested modules
    df = df[df[M["cluster"]].isin(modules)].copy()
    
    # For each module, find the most enriched pathway and its category/subcategory
    score_col = None
    for c in (M.get("qvalue"), M.get("padj"), M.get("p")):
        if c and c in df.columns:
            score_col = c; break
    
    if score_col:
        df["__score__"] = pd.to_numeric(df[score_col], errors="coerce")
        # Get the best (lowest score) pathway per module
        idx = df.groupby(M["cluster"])["__score__"].idxmin()
        best = df.loc[idx]
    else:
        # If no score, just take first pathway per module
        best = df.groupby(M["cluster"]).first()
    
    # Build result dataframe
    result = pd.DataFrame(index=pd.Index(modules, name="module"))
    result["category"] = np.nan
    result["subcategory"] = np.nan
    
    if M["category"]:
        cat_map = best.set_index(M["cluster"])[M["category"]].to_dict()
        result["category"] = result.index.map(cat_map)
    if M["subcategory"]:
        sub_map = best.set_index(M["cluster"])[M["subcategory"]].to_dict()
        result["subcategory"] = result.index.map(sub_map)
    
    return result


def draw_two_column_legend(
    ax: plt.Axes,
    cat_pal: Dict[str, Tuple],
    sub_pal: Dict[str, Tuple],
    alpha_levels: Optional[Tuple[Tuple[float, str], ...]] = None,
    title_left: str = "Category",
    title_right: str = "Subcategory",
    title_alpha: str = "Significance",
    fontsize: int = 10,
):
    """
    Draw compact legend using matplotlib's native legend function with small color patches.
    Similar style to enrichmentViz.py
    """
    ax.axis("off")
    
    # Create patch handles for categories and subcategories
    handles_cat = [Patch(facecolor=cat_pal[k], edgecolor='none', label=str(k)) for k in cat_pal]
    handles_sub = [Patch(facecolor=sub_pal[k], edgecolor='none', label=str(k)) for k in sub_pal]
    handles_sig: List[Patch] = []
    if alpha_levels:
        for thresh, symbol in alpha_levels:
            label = f"{symbol}: p ≤ {thresh:g}"
            handles_sig.append(
                Patch(facecolor="none", edgecolor="none", label=label)
            )
    
    def place_leg(ax, handles, title, y_top, fs):
        """Place a legend at the specified vertical position and return the next y position"""
        if not handles:
            return y_top
        # Estimate height needed (each item takes about 0.035 of the axis height)
        est_h = 0.03 * (len(handles) + 1.5)  # Compact spacing
        leg = ax.legend(
            handles=handles, 
            title=title, 
            loc="upper left",
            bbox_to_anchor=(0.0, y_top), 
            fontsize=fs, 
            title_fontsize=fs+1, 
            frameon=True,
            fancybox=False,
            edgecolor='gray',
            framealpha=0.8
        )
        ax.add_artist(leg)
        return max(0.02, y_top - est_h)
    
    # Place legends vertically with controlled gaps
    y = 0.98
    if handles_sig:
        y = place_leg(ax, handles_sig, title_alpha, y, fontsize)
        y = y - 0.04
    y = place_leg(ax, handles_cat, title_left, y, fontsize)
    y = y - 0.05  # Small vertical gap between category and subcategory legends
    place_leg(ax, handles_sub, title_right, y, fontsize)

# -----------------------------
# Clustering utilities
# -----------------------------

def cluster_modules(
    mat: pd.DataFrame,
    module_annotations: pd.DataFrame,
    cluster_by: str = "none",
) -> pd.Index:
    """
    Cluster modules by different strategies.
    
    Parameters
    ----------
    mat : pd.DataFrame
        Matrix with modules as rows, phenotypes as columns
    module_annotations : pd.DataFrame
        DataFrame with index=module_id, columns=[category, subcategory]
    cluster_by : str
        Clustering strategy: "none", "category", "subcategory", "category_then_subcategory"
    
    Returns
    -------
    pd.Index
        Reordered module indices
    """
    import scipy.cluster.hierarchy as sch
    from scipy.spatial.distance import pdist
    
    if cluster_by == "none":
        # Keep original order
        return mat.index
    
    elif cluster_by == "category":
        # Sort by category, then hierarchical clustering within each category
        annot = module_annotations.reindex(mat.index).copy()
        annot["category"] = annot["category"].fillna("Unknown")
        
        # Get unique categories in order of appearance
        cats = []
        for cat in annot["category"]:
            if cat not in cats:
                cats.append(cat)
        
        reordered = []
        for cat in cats:
            cat_mods = annot[annot["category"] == cat].index.tolist()
            if len(cat_mods) > 1:
                # Hierarchical clustering within category
                sub_mat = mat.loc[cat_mods]
                data = np.nan_to_num(sub_mat.values, nan=0.0)
                distances = pdist(data, metric='euclidean')
                linkage = sch.linkage(distances, method='average')
                dendro_idx = sch.dendrogram(linkage, no_plot=True)['leaves']
                reordered.extend([cat_mods[i] for i in dendro_idx])
            else:
                reordered.extend(cat_mods)
        
        return pd.Index(reordered)
    
    elif cluster_by == "subcategory":
        # Sort by subcategory, then hierarchical clustering within each subcategory
        annot = module_annotations.reindex(mat.index).copy()
        annot["subcategory"] = annot["subcategory"].fillna("Unknown")
        
        # Get unique subcategories in order of appearance
        subcats = []
        for subcat in annot["subcategory"]:
            if subcat not in subcats:
                subcats.append(subcat)
        
        reordered = []
        for subcat in subcats:
            sub_mods = annot[annot["subcategory"] == subcat].index.tolist()
            if len(sub_mods) > 1:
                # Hierarchical clustering within subcategory
                sub_mat = mat.loc[sub_mods]
                data = np.nan_to_num(sub_mat.values, nan=0.0)
                distances = pdist(data, metric='euclidean')
                linkage = sch.linkage(distances, method='average')
                dendro_idx = sch.dendrogram(linkage, no_plot=True)['leaves']
                reordered.extend([sub_mods[i] for i in dendro_idx])
            else:
                reordered.extend(sub_mods)
        
        return pd.Index(reordered)
    
    elif cluster_by == "category_then_subcategory":
        # Sort by category, then subcategory, then hierarchical clustering within each subcategory
        annot = module_annotations.reindex(mat.index).copy()
        annot["category"] = annot["category"].fillna("Unknown")
        annot["subcategory"] = annot["subcategory"].fillna("Unknown")
        
        # Get unique categories in order of appearance
        cats = []
        for cat in annot["category"]:
            if cat not in cats:
                cats.append(cat)
        
        reordered = []
        for cat in cats:
            cat_annot = annot[annot["category"] == cat]
            
            # Get subcategories within this category
            subcats = []
            for subcat in cat_annot["subcategory"]:
                if subcat not in subcats:
                    subcats.append(subcat)
            
            for subcat in subcats:
                sub_mods = cat_annot[cat_annot["subcategory"] == subcat].index.tolist()
                if len(sub_mods) > 1:
                    # Hierarchical clustering within subcategory
                    sub_mat = mat.loc[sub_mods]
                    data = np.nan_to_num(sub_mat.values, nan=0.0)
                    distances = pdist(data, metric='euclidean')
                    linkage = sch.linkage(distances, method='average')
                    dendro_idx = sch.dendrogram(linkage, no_plot=True)['leaves']
                    reordered.extend([sub_mods[i] for i in dendro_idx])
                else:
                    reordered.extend(sub_mods)
        
        return pd.Index(reordered)
    
    else:
        raise ValueError(f"Unknown cluster_by option: {cluster_by}. Use 'none', 'category', 'subcategory', or 'category_then_subcategory'")


# -----------------------------
# Plot heatmap (no seaborn, pure matplotlib)
# -----------------------------

def plot_heatmap_with_right_legend(
    mat: pd.DataFrame,
    padj_mat: Optional[pd.DataFrame],
    cat_pal: Dict[str, Tuple],
    sub_pal: Dict[str, Tuple],
    module_annotations: pd.DataFrame,
    out_pdf: str,
    cap: float = 1.0,
    figure_dpi: int = 160,
    figsize_scale: float = 1.0,
    alpha_levels: Tuple[Tuple[float, str], ...] = ((0.001, "***"), (0.01, "**"), (0.05, "*")),
    cluster_by: str = "none",
    transpose: bool = False,
    pad_inches: float = 1.0,
):
    """
    Plot phenotype × module heatmap with category/subcategory color strips on the left
    and a legend on the right. Adds a vertical color scale adjacent to the heatmap and a
    legend explaining significance star thresholds.
    
    Parameters
    ----------
    mat : pd.DataFrame
        Matrix with modules as rows, phenotypes as columns (signed correlation values)
    padj_mat : pd.DataFrame, optional
        Matrix with p-adjusted values for significance stars
    cat_pal : dict
        Category name -> color mapping
    sub_pal : dict
        Subcategory name -> color mapping
    module_annotations : pd.DataFrame
        DataFrame with index=module_id, columns=[category, subcategory]
    out_pdf : str
        Output PDF file path
    alpha_levels : tuple
        Significance thresholds and corresponding star symbols
    cluster_by : str
        Clustering strategy: "none", "category", "subcategory", "category_then_subcategory"
    transpose : bool
        If True, transpose the heatmap (phenotypes on Y-axis, modules on X-axis)
    pad_inches : float
        Extra padding (in inches) added around the saved PDF, increasing page size without
        scaling the plotted elements
    """
    from matplotlib.colors import to_rgb
    import matplotlib.patheffects as pe
    
    # Apply clustering
    reordered_idx = cluster_modules(mat, module_annotations, cluster_by)
    mat = mat.loc[reordered_idx]
    if padj_mat is not None:
        padj_mat = padj_mat.loc[reordered_idx]
    module_annotations = module_annotations.reindex(reordered_idx)
    
    # Transpose if requested
    if transpose:
        mat = mat.T
        if padj_mat is not None:
            padj_mat = padj_mat.T
    
    # Dynamic figure size - adaptive to content
    n_mod, n_pheno = mat.shape
    base_w = max(6.0, 0.20 * n_pheno)
    base_h = max(4.0, 0.30 * n_mod)  # Height per module
    fig_w = (base_w + 5.0) * figsize_scale   # +5 for the legend column
    fig_h = (base_h + 2.0) * figsize_scale  # Increased to 2.0 to accommodate legend

    fig = plt.figure(figsize=(fig_w, fig_h), dpi=figure_dpi)
    
    cbar_width = 0.25

    if transpose:
        # For transposed view: strips on top/bottom of heatmap
        # Strips should have same height ratio as heatmap for proper alignment
        strip_ratio = 4.0 / n_mod  # Proportional height per module
        gap = 0.015
        legend_width = 1.1
        
        gs = GridSpec(nrows=4, ncols=3,
                      height_ratios=[strip_ratio, strip_ratio, 4.0, gap],
                      width_ratios=[4.0, cbar_width, legend_width],
                      hspace=0.01, wspace=0.01, figure=fig)
        
        ax_cat = fig.add_subplot(gs[0,0])  # Category strip (top)
        ax_sub = fig.add_subplot(gs[1,0])  # Subcategory strip (below category)
        ax_hm  = fig.add_subplot(gs[2,0])  # Heatmap
        ax_cbar = fig.add_subplot(gs[2,1])  # Vertical colorbar next to heatmap
        ax_leg = fig.add_subplot(gs[:,2])  # Legend (full right column)
        
    else:
        # Original view: strips on right side of heatmap
        # Strips should have same width ratio as heatmap for proper alignment
        strip_width = 4.0 / n_mod  # Proportional width per module
        gap = 0.2
        legend_width = 1.1
        
        gs = GridSpec(nrows=1, ncols=6, 
                      width_ratios=[4.0, strip_width, strip_width, cbar_width, gap, legend_width], 
                      wspace=0.02, figure=fig)

        ax_hm  = fig.add_subplot(gs[0,0])  # Heatmap
        ax_cat = fig.add_subplot(gs[0,1])  # Category strip (right of heatmap)
        ax_sub = fig.add_subplot(gs[0,2])  # Subcategory strip (right of category strip)
        ax_cbar = fig.add_subplot(gs[0,3])  # Vertical colorbar after pathway strips
        ax_leg = fig.add_subplot(gs[0,5])  # Legend

    # Heatmap with diverging colormap (blue=negative, red=positive)
    data = np.nan_to_num(mat.values, nan=0.0)
    im = ax_hm.imshow(data, aspect="auto", vmin=-cap, vmax=cap, cmap="RdBu_r")

    # Axis labels
    ax_hm.set_yticks(range(n_mod))
    ax_hm.set_yticklabels(mat.index.astype(str), fontsize=8)
    ax_hm.set_xticks(range(n_pheno))
    ax_hm.set_xticklabels(mat.columns.astype(str), rotation=90, fontsize=8)

    # Significance stars (without colorbar)
    if padj_mat is not None:
        def get_star(p):
            if not np.isfinite(p):
                return ""
            for thresh, symbol in alpha_levels:
                if p <= thresh:
                    return symbol
            return ""
        
        for i in range(n_mod):
            for j in range(n_pheno):
                p_val = padj_mat.iloc[i, j]
                star = get_star(p_val)
                if star:
                    ax_hm.text(
                        j, i, star,
                        ha="center", va="center",
                        color="black", fontsize=8, weight="bold",
                        path_effects=[
                            pe.Stroke(linewidth=2, foreground="white"),
                            pe.Normal()
                        ]
                    )

    # Build color arrays for category and subcategory strips
    default_color = (0.92, 0.92, 0.92)  # Light gray for missing values
    
    if transpose:
        # For transposed view, modules are on X-axis, so we need horizontal strips
        def get_color_array(annot_col: str, palette: Dict[str, Tuple]) -> np.ndarray:
            colors = []
            # mat.columns now contains module IDs (after transpose)
            for mod_id in mat.columns:
                if mod_id in module_annotations.index:
                    val = module_annotations.loc[mod_id, annot_col]
                    if pd.notna(val):
                        val_str = str(val).strip()
                        color = palette.get(val_str, default_color)
                    else:
                        color = default_color
                else:
                    color = default_color
                colors.append(to_rgb(color))
            return np.array(colors)
        
        cat_colors = get_color_array("category", cat_pal)
        sub_colors = get_color_array("subcategory", sub_pal)
        
        # Draw horizontal color strips (1 row, n columns)
        # Use same extent as heatmap for perfect alignment
        n_cols = len(mat.columns)
        ax_cat.imshow(cat_colors.reshape(1, -1, 3), aspect="auto", extent=[-0.5, n_cols - 0.5, -0.5, 0.5])
        ax_sub.imshow(sub_colors.reshape(1, -1, 3), aspect="auto", extent=[-0.5, n_cols - 0.5, -0.5, 0.5])
        
        # Sync x-axis limits with heatmap for perfect alignment
        hm_xlim = (-0.5, n_cols - 0.5)
        for ax in (ax_cat, ax_sub):
            ax.set_xlim(hm_xlim)
            ax.set_ylim(-0.5, 0.5)
            ax.set_xticks([])
            ax.set_yticks([])
            for sp in ax.spines.values():
                sp.set_visible(False)
    else:
        # For normal view, modules are on Y-axis, so we need vertical strips
        def get_color_array(annot_col: str, palette: Dict[str, Tuple]) -> np.ndarray:
            colors = []
            for mod_id in mat.index:
                if mod_id in module_annotations.index:
                    val = module_annotations.loc[mod_id, annot_col]
                    if pd.notna(val):
                        val_str = str(val).strip()
                        color = palette.get(val_str, default_color)
                    else:
                        color = default_color
                else:
                    color = default_color
                colors.append(to_rgb(color))
            return np.array(colors)
        
        cat_colors = get_color_array("category", cat_pal)
        sub_colors = get_color_array("subcategory", sub_pal)
        
        # Draw vertical color strips (n rows, 1 column)
        # Use same extent as heatmap for perfect alignment
        n_rows = len(mat.index)
        ax_cat.imshow(cat_colors.reshape(-1, 1, 3), aspect="auto", extent=[-0.5, 0.5, n_rows - 0.5, -0.5])
        ax_sub.imshow(sub_colors.reshape(-1, 1, 3), aspect="auto", extent=[-0.5, 0.5, n_rows - 0.5, -0.5])
        
        # Sync y-axis limits with heatmap for perfect alignment
        hm_ylim = (n_rows - 0.5, -0.5)  # Inverted to match heatmap
        for ax in (ax_cat, ax_sub):
            ax.set_ylim(hm_ylim)
            ax.set_xlim(-0.5, 0.5)
            ax.set_xticks([])
            ax.set_yticks([])
            for sp in ax.spines.values():
                sp.set_visible(False)

    # Right legend for category/subcategory and significance thresholds
    draw_two_column_legend(ax_leg, cat_pal, sub_pal, alpha_levels=alpha_levels)

    # Add horizontal colorbar beneath the heatmap
    cbar = fig.colorbar(im, cax=ax_cbar, orientation="vertical")
    cbar.ax.tick_params(labelsize=8)

    # Save with tight layout to fit all content
    fig.savefig(out_pdf, bbox_inches="tight", pad_inches=pad_inches, format="pdf")
    plt.close(fig)

# -----------------------------
# Tables: module details, top phenos, top pathways
# -----------------------------

def _safe_int(x):
    try:
        return int(x)
    except Exception:
        return np.nan


def build_module_details(details_tsv: str, modules: Sequence[int]) -> pd.DataFrame:
    df = _read_table(details_tsv)
    M = _map_details_cols(df)
    if M["cluster_id"] is None:
        raise ValueError("DETAILS file must have a 'Cluster ID' / module id column")

    df = df.copy()
    df["__cluster_id__"] = pd.to_numeric(df[M["cluster_id"]], errors="coerce").astype("Int64")
    df = df[df["__cluster_id__"].notna()].copy()
    df["__cluster_id__"] = df["__cluster_id__"].astype(int)

    sel = df[df["__cluster_id__"].isin(list(modules))].copy()
    # Present nice canonical columns first if they exist
    wanted = [M.get("cluster_id"), M.get("cluster_type"), M.get("cluster_size"), M.get("tissue")]
    wanted = [w for w in wanted if w and w in sel.columns]
    other_cols = [c for c in sel.columns if c not in wanted + ["__cluster_id__"]]
    out = sel[wanted + other_cols]
    out = out.rename(columns={M.get("cluster_id"):"Cluster ID"} if M.get("cluster_id") else {})
    return out


def top_phenotypes(pheno_long: pd.DataFrame, top_k: int = 10) -> pd.DataFrame:
    M = _map_pheno_cols(pheno_long)
    if M["padj"] is None:
        raise ValueError("Phenotype table lacks an adjusted p column (padj/p.adjust).")
    df = pheno_long[[M["module"], M["phenotype"], M["padj"]]].copy()
    df[M["padj"]] = pd.to_numeric(df[M["padj"]], errors="coerce")
    df = df.dropna(subset=[M["padj"]])

    def take_top(g: pd.DataFrame) -> pd.DataFrame:
        gg = g.nsmallest(top_k, columns=[M["padj"]])
        return gg[[M["phenotype"], M["padj"]]]  # Return only the columns we need

    out = df.groupby(M["module"], group_keys=True).apply(take_top, include_groups=True).reset_index()
    # After groupby with include_groups=True, the module column is preserved
    out = out.rename(columns={M["module"]:"Module", M["phenotype"]:"Phenotype", M["padj"]:"p_adj"})
    out["Module"] = pd.to_numeric(out["Module"], errors="coerce").astype("Int64")
    out = out.sort_values(["Module","p_adj"], ascending=[True, True])
    # Drop any extra index columns
    out = out[[c for c in out.columns if c in ["Module", "Phenotype", "p_adj"]]]
    return out


def top_pathways_per_module(kegg_csv: str, top_k: int = 5) -> pd.DataFrame:
    df = _read_table(kegg_csv)
    M = _map_kegg_cols(df)
    need = ["cluster","desc"]
    for k in need:
        if M[k] is None:
            raise ValueError(f"KEGG file missing required column: {k}")

    df = df.copy()
    df[M["cluster"]] = pd.to_numeric(df[M["cluster"]], errors="coerce").astype("Int64")
    df = df[df[M["cluster"]].notna()].copy()
    df[M["cluster"]] = df[M["cluster"]].astype(int)
    df[M["desc"]] = df[M["desc"]].astype(str).str.strip()

    # Choose score: qvalue, fallback padj, fallback p
    score_col = None
    for c in (M.get("qvalue"), M.get("padj"), M.get("p")):
        if c and c in df.columns:
            score_col = c; break
    if score_col is None:
        raise ValueError("KEGG file lacks qvalue/p.adjust/pvalue columns")

    df["__score__"] = pd.to_numeric(df[score_col], errors="coerce")

    def take_top(g: pd.DataFrame) -> pd.DataFrame:
        gg = g.nsmallest(top_k, columns=["__score__"])
        # Return relevant columns, grouping column will be added by groupby
        keep_cols = [c for c in [M["desc"], M.get("category"), M.get("subcategory"), "__score__"] if c]
        return gg[keep_cols]

    out = (df.groupby(M["cluster"], group_keys=True)
             .apply(take_top, include_groups=True)
             .reset_index())

    rename_dict = {
        M["cluster"]: "Module",
        M["desc"]: "Pathway",
    }
    if M.get("category"):
        rename_dict[M["category"]] = "Category"
    if M.get("subcategory"):
        rename_dict[M["subcategory"]] = "Subcategory"
    # Always rename score_col to "score"
    if "__score__" in out.columns:
        rename_dict["__score__"] = "score"
    
    out = out.rename(columns=rename_dict)
    
    # Ensure these columns exist for downstream usage
    for col in ["Category","Subcategory"]:
        if col not in out.columns:
            out[col] = np.nan

    out["Module"] = pd.to_numeric(out["Module"], errors="coerce").astype("Int64")
    out = out.sort_values(["Module","score"], ascending=[True, True])
    
    # Select final columns
    final_cols = ["Module","Pathway"]
    if "Category" in out.columns:
        final_cols.append("Category")
    if "Subcategory" in out.columns:
        final_cols.append("Subcategory")
    final_cols.append("score")
    
    return out[final_cols]

# -----------------------------
# Pipeline runner
# -----------------------------

def run_all(
    pheno_tsv: str,
    kegg_csv: str,
    details_tsv: str,
    out_prefix: str = "heatmap",
    tissues: Optional[Sequence[str]] = None,
    include_phenotypes: Optional[Sequence[str]] = None,
    exclude_phenotypes: Optional[Sequence[str]] = None,
    metric: str = "corr",
    cap: float = 1.0,
    top_modules: Optional[int] = None,
    top_k_phenos: int = 10,
    top_k_pathways: int = 5,
    figure_dpi: int = 160,
    figsize_scale: float = 1.0,
    cluster_by: str = "none",
    transpose: bool = False,
    pdf_pad_inches: float = 5.5,
) -> Dict[str, str]:
    """
    Main pipeline to generate phenotype × module heatmap with pathway annotations.
    
    Parameters
    ----------
    top_modules : int, optional
        If provided, select only the top N modules ranked by maximum significance across phenotypes
    cluster_by : str
        Clustering strategy: "none" (by module order), "category", "subcategory", or "category_then_subcategory"
    transpose : bool
        If True, transpose the heatmap (phenotypes on Y-axis, modules on X-axis with strips on top/bottom)
    pdf_pad_inches : float
        Extra padding (in inches) around the saved PDF; larger values increase page size while preserving element scale
    """
    # 1) Build phenotype matrix
    mat, padj_mat, pheno_long = build_pheno_matrix(
        pheno_tsv=pheno_tsv,
        tissues=tissues,
        include_phenotypes=include_phenotypes,
        exclude_phenotypes=exclude_phenotypes,
        metric=metric,
        cap=cap,
        top_modules=top_modules,
    )
    if mat.empty:
        raise RuntimeError("After filtering, the phenotype × module matrix is empty.")

    # 2) Legend palettes from KEGG
    cat_pal, sub_pal = build_cat_palettes(kegg_csv)
    
    # 3) Get module annotations (category/subcategory from KEGG)
    modules = list(mat.index.astype(int))
    module_annotations = get_module_annotations(kegg_csv, modules)

    # 4) Plot heatmap with significance stars
    out_pdf = f"{out_prefix}_heatmap.pdf"
    plot_heatmap_with_right_legend(
        mat=mat, padj_mat=padj_mat,
        cat_pal=cat_pal, sub_pal=sub_pal,
        module_annotations=module_annotations,
        out_pdf=out_pdf, cap=cap, figure_dpi=figure_dpi, figsize_scale=figsize_scale,
        cluster_by=cluster_by,
        transpose=transpose,
        pad_inches=pdf_pad_inches,
    )

    # 5) Tables
    details_df = build_module_details(details_tsv, modules)
    details_path = f"{out_prefix}_module_details_selected.tsv"
    details_df.to_csv(details_path, sep="\t", index=False)

    phenos_df = top_phenotypes(pheno_long, top_k=top_k_phenos)
    phenos_path = f"{out_prefix}_top_phenotypes_per_module.tsv"
    phenos_df.to_csv(phenos_path, sep="\t", index=False)

    pathways_df = top_pathways_per_module(kegg_csv, top_k=top_k_pathways)
    pathways_path = f"{out_prefix}_top_pathways_per_module.tsv"
    pathways_df.to_csv(pathways_path, sep="\t", index=False)

    return {
        "heatmap_pdf": out_pdf,
        "module_details_tsv": details_path,
        "top_phenotypes_tsv": phenos_path,
        "top_pathways_tsv": pathways_path,
    }


if __name__ == "__main__":
    # Example: update paths to your local files or keep as uploaded ones
    outputs = run_all(
    pheno_tsv="/Users/edeneldar/CoExpression_ReProduction/notebooks/rosmap_ME_pheno_ME_vs_pheno_correlations.tsv",
    kegg_csv="/Volumes/Transcend/Eden/CoExpression_ReProduction/kegg_rosmap_constBeta_CT2_TS3.csv",
    details_tsv="/Volumes/Transcend/Eden/CoExpression_ReProduction/nbs/xwgcna_rosmap_constBeta_CT2_TS3_Cluster_details.tsv",
    out_prefix="rosmap_constBeta_CT2_TS3",
    tissues=['CROSS'],                 # e.g., ["AC","MFBA9BA46","PCGBA23"]
    include_phenotypes=None,      # e.g., [r"cog|memory", r"age"]
    exclude_phenotypes=[r"msex"],      # e.g., [r"sex|gender"]
    metric="corr",                # Changed to "corr" for signed correlation values
    cap=1.0,                      # Changed to 1.0 for correlation range [-1, 1]
    top_modules=20,             # e.g., 100 to select only top 100 modules by phenotype significance
    top_k_phenos=10,
    top_k_pathways=5,
    figure_dpi=180,
    figsize_scale=1.0,
    cluster_by="category_then_subcategory",  # Options: "none", "category", "subcategory", "category_then_subcategory"
    transpose=False,              # Set to True for phenotypes on Y-axis, modules on X-axis
)
    print(json.dumps(outputs, indent=2))
