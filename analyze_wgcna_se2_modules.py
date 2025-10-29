#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
analyze_wgcna_se2_modules.py
Analysis and visualization of overlaps between WGCNA modules and SpeakEasy2 communities.

Inputs:
  --wgcna-file: Path to WGCNA modules table TSV
  --se2-file: Path to SpeakEasy2 clusters table CSV
  --outdir: Output directory for results

Outputs:
  - module_pair_metrics.csv: Detailed metrics for each module-community pair
  - top_intersections_gene_lists.csv: Gene lists for top 30 intersections
  - bar_top_pairs.png: Bar plot of top N pairs by Jaccard index
  - bar_module_topk.png: Bar plot showing top K communities for each module
"""

import argparse
import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from typing import Dict, Set, List, Tuple
import re

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
def prepare_gene_data(wgcna_df: pd.DataFrame, se2_df: pd.DataFrame) -> Tuple[Dict, Dict, Dict, Dict]:
    """
    Prepare gene data for analysis by creating dictionaries mapping modules to gene sets.
    Returns:
        wgcna_genes: Dict mapping WGCNA module IDs to sets of genes
        se2_genes: Dict mapping SE2 community IDs to sets of genes
        wgcna_module_sizes: Dict with size of each WGCNA module
        se2_community_sizes: Dict with size of each SE2 community
    """
    # For WGCNA, group genes by Cluster ID and Tissue
    wgcna_genes = {}
    for tissue in wgcna_df['Tissue'].unique():
        tissue_df = wgcna_df[wgcna_df['Tissue'] == tissue]
        for cluster_id in tissue_df['Cluster ID'].unique():
            # Use cluster_id_tissue as the key to distinguish same IDs across tissues
            key = f"{cluster_id}_{tissue}"
            gene_set = set(tissue_df[tissue_df['Cluster ID'] == cluster_id]['Gene Symbol'])
            wgcna_genes[key] = gene_set
    
    # For SpeakEasy2, group genes by Cluster.ID and Tissue
    se2_genes = {}
    for tissue in se2_df['Tissue'].unique():
        tissue_df = se2_df[se2_df['Tissue'] == tissue]
        for cluster_id in tissue_df['Cluster.ID'].unique():
            # Use cluster_id_tissue as the key
            key = f"{cluster_id}_{tissue}"
            gene_set = set(tissue_df[tissue_df['Cluster.ID'] == cluster_id]['Gene.Symbol'])
            se2_genes[key] = gene_set
    
    # Calculate sizes
    wgcna_module_sizes = {module: len(genes) for module, genes in wgcna_genes.items()}
    se2_community_sizes = {comm: len(genes) for comm, genes in se2_genes.items()}
    
    return wgcna_genes, se2_genes, wgcna_module_sizes, se2_community_sizes

def calculate_module_pair_metrics(wgcna_genes: Dict, se2_genes: Dict) -> pd.DataFrame:
    """
    Calculate metrics for all pairs of WGCNA modules and SE2 communities.
    Returns DataFrame with metrics for each pair.
    """
    pairs = []
    
    for wgcna_id, wgcna_set in wgcna_genes.items():
        wgcna_module_id = wgcna_id.split('_')[0]  # Extract numeric ID
        wgcna_tissue = wgcna_id.split('_')[1]     # Extract tissue
        
        for se2_id, se2_set in se2_genes.items():
            se2_community_id = se2_id.split('_')[0]  # Extract numeric ID
            se2_tissue = se2_id.split('_')[1]        # Extract tissue
            
            # Calculate intersection and union
            intersection = wgcna_set.intersection(se2_set)
            union = wgcna_set.union(se2_set)
            
            # Calculate Jaccard index
            jaccard = len(intersection) / len(union) if len(union) > 0 else 0
            
            pairs.append({
                "wgcna_module": wgcna_module_id,
                "wgcna_tissue": wgcna_tissue,
                "se2_community": se2_community_id,
                "se2_tissue": se2_tissue,
                "wgcna_size": len(wgcna_set),
                "se2_size": len(se2_set),
                "intersection_size": len(intersection),
                "union_size": len(union),
                "jaccard_index": jaccard,
                "same_tissue": wgcna_tissue == se2_tissue,
                "pair_id": f"W:{wgcna_module_id}_{wgcna_tissue} → S:{se2_community_id}_{se2_tissue}"
            })
    
    return pd.DataFrame(pairs)

def generate_top_intersection_gene_lists(
    metrics_df: pd.DataFrame, 
    wgcna_genes: Dict, 
    se2_genes: Dict, 
    top_n: int = 30
) -> pd.DataFrame:
    """
    Generate a table with gene lists for the top N pairs with largest intersections.
    """
    # Sort by intersection size, descending
    top_pairs = metrics_df.sort_values("intersection_size", ascending=False).head(top_n)
    
    gene_lists = []
    for _, row in top_pairs.iterrows():
        wgcna_id = f"{row['wgcna_module']}_{row['wgcna_tissue']}"
        se2_id = f"{row['se2_community']}_{row['se2_tissue']}"
        
        wgcna_set = wgcna_genes.get(wgcna_id, set())
        se2_set = se2_genes.get(se2_id, set())
        intersection = wgcna_set.intersection(se2_set)
        
        gene_lists.append({
            "wgcna_module": row['wgcna_module'],
            "wgcna_tissue": row['wgcna_tissue'],
            "se2_community": row['se2_community'],
            "se2_tissue": row['se2_tissue'],
            "intersection_size": len(intersection),
            "intersection_genes": ";".join(sorted(intersection)),
            "wgcna_only_genes": ";".join(sorted(wgcna_set - intersection)),
            "se2_only_genes": ";".join(sorted(se2_set - intersection))
        })
    
    return pd.DataFrame(gene_lists)

# ---------- Visualization functions ----------
def barplot_top_pairs(df: pd.DataFrame, out_png: str, metric: str = "jaccard_index", top_n: int = 150):
    """Generate bar plot showing top N pairs by specified metric."""
    data = df.sort_values(metric, ascending=False).head(top_n)
    n = len(data)
    
    fig, ax = plt.subplots(figsize=(max(8, 0.35*n), 6))
    ax.bar(range(n), data[metric].values)
    ax.set_xticks(range(n))
    ax.set_xticklabels(data["pair_id"].values, rotation=90)
    ax.set_title(f"Top {top_n} module-community pairs by {metric}")
    ax.set_ylabel(metric)
    fig.tight_layout()
    fig.savefig(out_png, dpi=160)
    plt.close(fig)

def barplot_per_module_topk(df: pd.DataFrame, out_png: str, metric: str = "jaccard_index", top_modules: int = 20, k: int = 3):
    """Generate bar plot showing top K communities for each of the top modules."""
    # Select top modules based on their max metric score
    order_mod = (df.groupby(["wgcna_module", "wgcna_tissue"])[metric]
                  .max()
                  .sort_values(ascending=False)
                  .head(top_modules)
                  .index.tolist())
    
    data = []
    for wm, wt in order_mod:
        # Get top K communities for this module
        module_data = df[(df["wgcna_module"] == wm) & (df["wgcna_tissue"] == wt)]
        module_top_k = module_data.sort_values(metric, ascending=False).head(k)
        data.append(module_top_k)
    
    if not data:  # Handle case where there are no matches
        return
        
    sub = pd.concat(data)
    sub = sub.sort_values(["wgcna_module", "wgcna_tissue", metric], ascending=[True, True, False])
    
    # Create bar plot
    n = len(sub)
    fig, ax = plt.subplots(figsize=(max(10, 0.4*n), 7))
    ax.bar(range(n), sub[metric].values)
    ax.set_xticks(range(n))
    ax.set_xticklabels(sub["pair_id"].values, rotation=90)
    ax.set_title(f"Top {k} communities per module (showing {len(order_mod)} modules)")
    ax.set_ylabel(metric)
    fig.tight_layout()
    fig.savefig(out_png, dpi=160)
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
    
    args = parser.parse_args()
    
    os.makedirs(args.outdir, exist_ok=True)
    
    print("Loading input files...")
    wgcna_df = read_wgcna_modules(args.wgcna_file)
    se2_df = read_se2_clusters(args.se2_file)
    
    print(f"WGCNA data: {wgcna_df.shape[0]} genes across {wgcna_df['Cluster ID'].nunique()} modules")
    print(f"SE2 data: {se2_df.shape[0]} genes across {se2_df['Cluster.ID'].nunique()} communities")
    
    # Inspect the first few rows of each dataset
    print("\nWGCNA data sample:")
    print(wgcna_df.head())
    print("\nSE2 data sample:")
    print(se2_df.head())
    
    # Get unique tissues in each dataset
    wgcna_tissues = wgcna_df['Tissue'].unique()
    se2_tissues = se2_df['Tissue'].unique()
    print(f"\nTissues in WGCNA data: {', '.join(wgcna_tissues)}")
    print(f"Tissues in SE2 data: {', '.join(se2_tissues)}")
    
    # Prepare gene sets and calculate metrics
    print("\nPreparing gene data and calculating metrics...")
    wgcna_genes, se2_genes, wgcna_sizes, se2_sizes = prepare_gene_data(wgcna_df, se2_df)
    
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
    metrics_df = calculate_module_pair_metrics(wgcna_genes, se2_genes)
    
    print(f"Generated metrics for {len(metrics_df)} module-community pairs")
    print(f"Same tissue pairs: {metrics_df['same_tissue'].sum()} out of {len(metrics_df)}")
    
    # Generate gene lists for top intersections
    print("Generating gene lists for top intersections...")
    gene_lists_df = generate_top_intersection_gene_lists(metrics_df, wgcna_genes, se2_genes)
    
    # Save results
    print("Saving results...")
    metrics_df.to_csv(os.path.join(args.outdir, "module_pair_metrics.csv"), index=False)
    gene_lists_df.to_csv(os.path.join(args.outdir, "top_intersections_gene_lists.csv"), index=False)
    
    # Generate visualizations
    print("Generating visualizations...")
    barplot_top_pairs(metrics_df, 
                     os.path.join(args.outdir, "bar_top_pairs.png"), 
                     metric="jaccard_index", 
                     top_n=args.top_n)
    
    barplot_per_module_topk(metrics_df, 
                           os.path.join(args.outdir, "bar_module_topk.png"), 
                           metric="jaccard_index",
                           top_modules=args.top_modules, 
                           k=args.k)
    
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
    print(f"Pairs with same tissue: {metrics_df['same_tissue'].sum()} ({metrics_df['same_tissue'].mean()*100:.1f}%)")
    
    print("\nDone! Results saved to:", args.outdir)

if __name__ == "__main__":
    main()