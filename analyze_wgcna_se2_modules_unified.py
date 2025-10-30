#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
analyze_wgcna_se2_modules_unified.py
Analysis and visualization of overlaps between WGCNA modules and SpeakEasy2 communities,
treating modules as unified entities with tissue-specific gene annotations.

Inputs:
  --wgcna-file: Path to WGCNA modules table TSV
  --se2-file: Path to SpeakEasy2 clusters table CSV
  --outdir: Output directory for results

Outputs:
  - module_pair_metrics_unified.csv: Detailed metrics for each module-community pair
  - top_intersections_gene_lists_unified.csv: Gene lists for top 30 intersections with tissue annotations
  - bar_top_pairs_unified.png: Bar plot of top N pairs by Jaccard index
  - tissue_composition_per_module.csv: Tissue composition for each module
  - tissue_composition_plots/: Directory with tissue composition plots for each module
"""

import argparse
import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from typing import Dict, Set, List, Tuple, Any
import re
from collections import defaultdict, Counter

# ---------- IO helpers ----------
def read_wgcna_modules(path: str) -> pd.DataFrame:
    """Read WGCNA module table from TSV file."""
    df = pd.read_csv(path, sep='\t')
    return df

def read_se2_clusters(path: str) -> pd.DataFrame:
    """Read SpeakEasy2 cluster table from CSV file."""
    df = pd.read_csv(path)
    return df

def normalize_gene_id(gene_id: str) -> str:
    """
    Normalize gene ID by removing version numbers but preserving tissue information.
    E.g., ENSG00000123456.10 -> ENSG00000123456
    But ENSG00000123456.10_AC -> ENSG00000123456_AC
    """
    gene_id = str(gene_id)
    
    # If the gene ID has a tissue suffix after an underscore
    if '_' in gene_id:
        base_id, tissue = gene_id.split('_', 1)
        if '.' in base_id:  # Remove version number if present
            base_id = base_id.split('.')[0]
        return f"{base_id}_{tissue}"
    
    # If the gene ID has a version number
    if '.' in gene_id:
        return gene_id.split('.')[0]
    
    return gene_id

# ---------- Analysis functions ----------
def prepare_unified_gene_data(wgcna_df: pd.DataFrame, se2_df: pd.DataFrame) -> Tuple[Dict, Dict, Dict, Dict]:
    """
    Prepare gene data for analysis by creating dictionaries mapping modules to gene sets.
    Treats modules as unified entities across tissues, but preserves tissue information for each gene.
    
    Returns:
        wgcna_genes: Dict mapping WGCNA module IDs to sets of tissue-annotated genes
        se2_genes: Dict mapping SE2 community IDs to sets of tissue-annotated genes
        wgcna_module_sizes: Dict with size of each WGCNA module
        se2_community_sizes: Dict with size of each SE2 community
    """
    # For WGCNA, group genes by Cluster ID (ignoring tissue for the module identification)
    wgcna_genes = defaultdict(set)
    
    for _, row in wgcna_df.iterrows():
        module_id = str(row['Cluster ID'])
        tissue = row['Tissue']
        gene = str(row['Gene Symbol'])
        
        # Store the gene with tissue information
        gene_with_tissue = f"{gene}_{tissue}"
        wgcna_genes[module_id].add(gene_with_tissue)
    
    # For SpeakEasy2, group genes by Cluster.ID (ignoring tissue for the module identification)
    se2_genes = defaultdict(set)
    
    for _, row in se2_df.iterrows():
        module_id = str(row['Cluster.ID'])
        tissue = row['Tissue']
        gene = str(row['Gene.Symbol'])
        
        # Store the gene with tissue information
        gene_with_tissue = f"{gene}_{tissue}"
        se2_genes[module_id].add(gene_with_tissue)
    
    # Calculate sizes
    wgcna_module_sizes = {module: len(genes) for module, genes in wgcna_genes.items()}
    se2_community_sizes = {comm: len(genes) for comm, genes in se2_genes.items()}
    
    return wgcna_genes, se2_genes, wgcna_module_sizes, se2_community_sizes

def analyze_tissue_composition(gene_set: Set[str]) -> Dict[str, int]:
    """
    Analyze the tissue composition of a gene set.
    Returns a dictionary mapping tissue to count of genes.
    """
    tissue_counts = Counter()
    
    for gene_with_tissue in gene_set:
        if '_' in gene_with_tissue:
            tissue = gene_with_tissue.split('_')[-1]
            tissue_counts[tissue] += 1
    
    return tissue_counts

def calculate_unified_module_pair_metrics(wgcna_genes: Dict, se2_genes: Dict) -> pd.DataFrame:
    """
    Calculate metrics for all pairs of WGCNA modules and SE2 communities.
    Returns DataFrame with metrics for each pair.
    """
    pairs = []
    
    for wgcna_id, wgcna_set in wgcna_genes.items():
        # Analyze tissue composition
        wgcna_tissue_comp = analyze_tissue_composition(wgcna_set)
        wgcna_dominant_tissue = max(wgcna_tissue_comp.items(), key=lambda x: x[1])[0] if wgcna_tissue_comp else "unknown"
        
        for se2_id, se2_set in se2_genes.items():
            # Analyze tissue composition
            se2_tissue_comp = analyze_tissue_composition(se2_set)
            se2_dominant_tissue = max(se2_tissue_comp.items(), key=lambda x: x[1])[0] if se2_tissue_comp else "unknown"
            
            # Calculate intersection and union
            intersection = wgcna_set.intersection(se2_set)
            union = wgcna_set.union(se2_set)
            
            # Calculate Jaccard index
            jaccard = len(intersection) / len(union) if len(union) > 0 else 0
            
            # Analyze tissue composition of intersection
            intersection_tissue_comp = analyze_tissue_composition(intersection)
            
            # Convert tissue compositions to strings for easier storage
            wgcna_tissue_str = "; ".join([f"{t}:{c}" for t, c in wgcna_tissue_comp.items()])
            se2_tissue_str = "; ".join([f"{t}:{c}" for t, c in se2_tissue_comp.items()])
            intersection_tissue_str = "; ".join([f"{t}:{c}" for t, c in intersection_tissue_comp.items()])
            
            pairs.append({
                "wgcna_module": wgcna_id,
                "wgcna_dominant_tissue": wgcna_dominant_tissue,
                "wgcna_tissue_composition": wgcna_tissue_str,
                "se2_community": se2_id,
                "se2_dominant_tissue": se2_dominant_tissue,
                "se2_tissue_composition": se2_tissue_str,
                "wgcna_size": len(wgcna_set),
                "se2_size": len(se2_set),
                "intersection_size": len(intersection),
                "intersection_tissue_composition": intersection_tissue_str,
                "union_size": len(union),
                "jaccard_index": jaccard,
                "same_dominant_tissue": wgcna_dominant_tissue == se2_dominant_tissue,
                "pair_id": f"W:{wgcna_id} → S:{se2_id}"
            })
    
    return pd.DataFrame(pairs)

def generate_unified_top_intersection_gene_lists(
    metrics_df: pd.DataFrame, 
    wgcna_genes: Dict, 
    se2_genes: Dict, 
    top_n: int = 30
) -> pd.DataFrame:
    """
    Generate a table with gene lists for the top N pairs with largest intersections.
    Preserves tissue information for each gene.
    """
    # Sort by intersection size, descending
    top_pairs = metrics_df.sort_values("intersection_size", ascending=False).head(top_n)
    
    gene_lists = []
    for _, row in top_pairs.iterrows():
        wgcna_id = row['wgcna_module']
        se2_id = row['se2_community']
        
        wgcna_set = wgcna_genes.get(wgcna_id, set())
        se2_set = se2_genes.get(se2_id, set())
        intersection = wgcna_set.intersection(se2_set)
        
        # Group genes by tissue in the intersection
        intersection_by_tissue = defaultdict(list)
        for gene_with_tissue in intersection:
            gene_base, tissue = gene_with_tissue.rsplit('_', 1)
            intersection_by_tissue[tissue].append(gene_base)
        
        # Format tissue-specific gene lists
        tissue_specific_lists = {}
        for tissue, genes in intersection_by_tissue.items():
            tissue_specific_lists[f"intersection_genes_{tissue}"] = ";".join(sorted(genes))
        
        gene_lists.append({
            "wgcna_module": wgcna_id,
            "wgcna_dominant_tissue": row['wgcna_dominant_tissue'],
            "se2_community": se2_id,
            "se2_dominant_tissue": row['se2_dominant_tissue'],
            "intersection_size": len(intersection),
            "intersection_genes": ";".join(sorted(intersection)),
            "wgcna_only_genes": ";".join(sorted(wgcna_set - intersection)),
            "se2_only_genes": ";".join(sorted(se2_set - intersection)),
            **tissue_specific_lists
        })
    
    return pd.DataFrame(gene_lists)

def analyze_tissue_composition_per_module(wgcna_genes: Dict, se2_genes: Dict) -> pd.DataFrame:
    """
    Analyze tissue composition for each module/community.
    Returns a DataFrame with detailed tissue composition.
    """
    rows = []
    
    # Process WGCNA modules
    for module_id, gene_set in wgcna_genes.items():
        tissue_comp = analyze_tissue_composition(gene_set)
        total_genes = len(gene_set)
        
        for tissue, count in tissue_comp.items():
            rows.append({
                "module_type": "WGCNA",
                "module_id": module_id,
                "tissue": tissue,
                "gene_count": count,
                "percentage": (count / total_genes) * 100 if total_genes > 0 else 0,
                "total_genes": total_genes
            })
    
    # Process SE2 communities
    for comm_id, gene_set in se2_genes.items():
        tissue_comp = analyze_tissue_composition(gene_set)
        total_genes = len(gene_set)
        
        for tissue, count in tissue_comp.items():
            rows.append({
                "module_type": "SpeakEasy2",
                "module_id": comm_id,
                "tissue": tissue,
                "gene_count": count,
                "percentage": (count / total_genes) * 100 if total_genes > 0 else 0,
                "total_genes": total_genes
            })
    
    return pd.DataFrame(rows)

# ---------- Visualization functions ----------
def barplot_top_pairs_unified(df: pd.DataFrame, out_png: str, metric: str = "jaccard_index", top_n: int = 150):
    """Generate bar plot showing top N pairs by specified metric."""
    data = df.sort_values(metric, ascending=False).head(top_n)
    n = len(data)
    
    fig, ax = plt.subplots(figsize=(max(8, 0.35*n), 6))
    bars = ax.bar(range(n), data[metric].values)
    
    # Color bars by tissue match
    for i, (_, row) in enumerate(data.iterrows()):
        if row['same_dominant_tissue']:
            bars[i].set_color('forestgreen')
        else:
            bars[i].set_color('lightcoral')
    
    ax.set_xticks(range(n))
    ax.set_xticklabels(data["pair_id"].values, rotation=90)
    ax.set_title(f"Top {top_n} module-community pairs by {metric}\n(green: same dominant tissue, red: different dominant tissues)")
    ax.set_ylabel(metric)
    fig.tight_layout()
    fig.savefig(out_png, dpi=160)
    plt.close(fig)

def barplot_per_module_topk_unified(df: pd.DataFrame, out_png: str, metric: str = "jaccard_index", top_modules: int = 20, k: int = 3):
    """Generate bar plot showing top K communities for each of the top modules."""
    # Select top modules based on their max metric score
    order_mod = (df.groupby(["wgcna_module"])[metric]
                  .max()
                  .sort_values(ascending=False)
                  .head(top_modules)
                  .index.tolist())
    
    data = []
    for wm in order_mod:
        # Get top K communities for this module
        module_data = df[df["wgcna_module"] == wm]
        module_top_k = module_data.sort_values(metric, ascending=False).head(k)
        data.append(module_top_k)
    
    if not data:  # Handle case where there are no matches
        return
        
    sub = pd.concat(data)
    sub = sub.sort_values(["wgcna_module", metric], ascending=[True, False])
    
    # Create bar plot
    n = len(sub)
    fig, ax = plt.subplots(figsize=(max(10, 0.4*n), 7))
    bars = ax.bar(range(n), sub[metric].values)
    
    # Color bars by tissue match
    for i, (_, row) in enumerate(sub.iterrows()):
        if row['same_dominant_tissue']:
            bars[i].set_color('forestgreen')
        else:
            bars[i].set_color('lightcoral')
    
    ax.set_xticks(range(n))
    ax.set_xticklabels(sub["pair_id"].values, rotation=90)
    ax.set_title(f"Top {k} communities per module (showing {len(order_mod)} modules)\n(green: same dominant tissue, red: different dominant tissues)")
    ax.set_ylabel(metric)
    fig.tight_layout()
    fig.savefig(out_png, dpi=160)
    plt.close(fig)

def plot_tissue_composition(df: pd.DataFrame, module_type: str, module_id: str, out_dir: str):
    """Generate a pie chart showing tissue composition for a single module."""
    module_df = df[(df['module_type'] == module_type) & (df['module_id'] == module_id)]
    
    if len(module_df) == 0:
        return
    
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.pie(module_df['gene_count'], labels=module_df['tissue'], autopct='%1.1f%%')
    ax.set_title(f"Tissue Composition for {module_type} Module {module_id}\nTotal: {module_df['total_genes'].iloc[0]} genes")
    
    # Save the plot
    os.makedirs(out_dir, exist_ok=True)
    out_file = os.path.join(out_dir, f"{module_type}_{module_id}_tissue_composition.png")
    fig.savefig(out_file, dpi=160)
    plt.close(fig)

# ---------- Main ----------
def main():
    parser = argparse.ArgumentParser(description="Analyze overlap between WGCNA modules and SpeakEasy2 communities")
    parser.add_argument("--wgcna-file", required=True, help="Path to WGCNA modules table TSV")
    parser.add_argument("--se2-file", required=True, help="Path to SpeakEasy2 clusters table CSV")
    parser.add_argument("--outdir", required=True, help="Output directory for results")
    parser.add_argument("--top-n", type=int, default=150, help="Number of top pairs to show in overall bar plot")
    parser.add_argument("--top-modules", type=int, default=20, help="Number of top modules to show in per-module bar plot")
    parser.add_argument("--k", type=int, default=3, help="Number of top communities to show per module")
    parser.add_argument("--min-module-size", type=int, default=0, help="Minimum WGCNA module size to include")
    parser.add_argument("--min-community-size", type=int, default=0, help="Minimum SE2 community size to include")
    parser.add_argument("--plot-all-tissue-compositions", action="store_true", help="Generate tissue composition plots for all modules")
    
    args = parser.parse_args()
    
    os.makedirs(args.outdir, exist_ok=True)
    tissue_plots_dir = os.path.join(args.outdir, "tissue_composition_plots")
    os.makedirs(tissue_plots_dir, exist_ok=True)
    
    print("Loading input files...")
    wgcna_df = read_wgcna_modules(args.wgcna_file)
    se2_df = read_se2_clusters(args.se2_file)
    
    print(f"WGCNA data: {wgcna_df.shape[0]} genes across {wgcna_df['Cluster ID'].nunique()} modules")
    print(f"SE2 data: {se2_df.shape[0]} genes across {se2_df['Cluster.ID'].nunique()} communities")
    
    # Get unique tissues in each dataset
    wgcna_tissues = wgcna_df['Tissue'].unique()
    se2_tissues = se2_df['Tissue'].unique()
    print(f"\nTissues in WGCNA data: {', '.join(wgcna_tissues)}")
    print(f"Tissues in SE2 data: {', '.join(se2_tissues)}")
    
    # Prepare gene sets and calculate metrics
    print("\nPreparing unified gene data and calculating metrics...")
    wgcna_genes, se2_genes, wgcna_sizes, se2_sizes = prepare_unified_gene_data(wgcna_df, se2_df)
    
    # Apply size filters if specified
    if args.min_module_size > 0:
        original_count = len(wgcna_genes)
        wgcna_genes = {module: genes for module, genes in wgcna_genes.items() 
                       if len(genes) >= args.min_module_size}
        print(f"Filtered WGCNA modules by size: {original_count} -> {len(wgcna_genes)}")
    
    if args.min_community_size > 0:
        original_count = len(se2_genes)
        se2_genes = {comm: genes for comm, genes in se2_genes.items() 
                    if len(genes) >= args.min_community_size}
        print(f"Filtered SE2 communities by size: {original_count} -> {len(se2_genes)}")
    
    # Calculate metrics for all module-community pairs
    print("Calculating metrics for all module-community pairs...")
    metrics_df = calculate_unified_module_pair_metrics(wgcna_genes, se2_genes)
    
    print(f"Generated metrics for {len(metrics_df)} module-community pairs")
    print(f"Same dominant tissue pairs: {metrics_df['same_dominant_tissue'].sum()} out of {len(metrics_df)}")
    
    # Generate gene lists for top intersections
    print("Generating gene lists for top intersections...")
    gene_lists_df = generate_unified_top_intersection_gene_lists(metrics_df, wgcna_genes, se2_genes)
    
    # Analyze tissue composition for each module/community
    print("Analyzing tissue composition for each module/community...")
    tissue_comp_df = analyze_tissue_composition_per_module(wgcna_genes, se2_genes)
    
    # Save results
    print("Saving results...")
    metrics_df.to_csv(os.path.join(args.outdir, "module_pair_metrics_unified.csv"), index=False)
    gene_lists_df.to_csv(os.path.join(args.outdir, "top_intersections_gene_lists_unified.csv"), index=False)
    tissue_comp_df.to_csv(os.path.join(args.outdir, "tissue_composition_per_module.csv"), index=False)
    
    # Generate visualizations
    print("Generating visualizations...")
    barplot_top_pairs_unified(metrics_df, 
                             os.path.join(args.outdir, "bar_top_pairs_unified.png"), 
                             metric="jaccard_index", 
                             top_n=args.top_n)
    
    barplot_per_module_topk_unified(metrics_df, 
                                   os.path.join(args.outdir, "bar_module_topk_unified.png"), 
                                   metric="jaccard_index",
                                   top_modules=args.top_modules, 
                                   k=args.k)
    
    # Generate tissue composition plots for top modules
    print("Generating tissue composition plots...")
    
    # Get top modules based on intersection size
    top_wgcna_modules = (metrics_df.groupby('wgcna_module')['intersection_size']
                         .max().sort_values(ascending=False)
                         .head(10).index.tolist())
    
    top_se2_communities = (metrics_df.groupby('se2_community')['intersection_size']
                          .max().sort_values(ascending=False)
                          .head(10).index.tolist())
    
    print(f"Creating tissue composition plots for top {len(top_wgcna_modules)} WGCNA modules...")
    for module_id in top_wgcna_modules:
        plot_tissue_composition(tissue_comp_df, "WGCNA", module_id, tissue_plots_dir)
    
    print(f"Creating tissue composition plots for top {len(top_se2_communities)} SpeakEasy2 communities...")
    for comm_id in top_se2_communities:
        plot_tissue_composition(tissue_comp_df, "SpeakEasy2", comm_id, tissue_plots_dir)
    
    # Optionally plot all module compositions
    if args.plot_all_tissue_compositions:
        print("Creating tissue composition plots for all modules (this may take a while)...")
        for module_id in wgcna_genes.keys():
            plot_tissue_composition(tissue_comp_df, "WGCNA", module_id, tissue_plots_dir)
        
        for comm_id in se2_genes.keys():
            plot_tissue_composition(tissue_comp_df, "SpeakEasy2", comm_id, tissue_plots_dir)
    
    # Summary statistics
    print("\nSummary Statistics:")
    print(f"Total WGCNA modules: {len(wgcna_genes)}")
    print(f"Total SE2 communities: {len(se2_genes)}")
    print(f"Total module-community pairs: {len(metrics_df)}")
    
    # Statistics on Jaccard indices
    jaccard_mean = metrics_df['jaccard_index'].mean()
    jaccard_median = metrics_df['jaccard_index'].median()
    high_jaccard_pairs = (metrics_df['jaccard_index'] > 0.5).sum()
    
    print(f"Jaccard index - Mean: {jaccard_mean:.4f}, Median: {jaccard_median:.4f}")
    print(f"Pairs with high Jaccard (>0.5): {high_jaccard_pairs}")
    print(f"Pairs with same dominant tissue: {metrics_df['same_dominant_tissue'].sum()} ({metrics_df['same_dominant_tissue'].mean()*100:.1f}%)")
    
    # Tissue statistics
    tissue_counts = tissue_comp_df.groupby(['module_type', 'tissue'])['gene_count'].sum().reset_index()
    print("\nTissue distribution across modules:")
    for module_type in ["WGCNA", "SpeakEasy2"]:
        type_counts = tissue_counts[tissue_counts['module_type'] == module_type]
        print(f"  {module_type}:")
        for _, row in type_counts.iterrows():
            print(f"    {row['tissue']}: {row['gene_count']} genes")
    
    print("\nDone! Results saved to:", args.outdir)

if __name__ == "__main__":
    main()