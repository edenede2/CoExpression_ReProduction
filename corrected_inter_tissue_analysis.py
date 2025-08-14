#!/usr/bin/env python3
"""
Corrected Inter-Tissue Network Analysis
This script demonstrates the correct format for inter-tissue analysis that matches the R script.
"""

import pickle
from pathlib import Path
from rich import print as rprint
from new_script.network_analysis import NetworkAnalysis

def main():
    # Load preprocessed tissue data
    with open("tissue_processed_data.pkl", "rb") as f:
        tissue_processed_data = pickle.load(f)
    
    # Load sample mappings (you'll need to run the preprocessing notebook first)
    target_tissues = ['Brain - Cortex', 'Muscle - Skeletal', 'Adipose - Subcutaneous']
    
    # Initialize network analysis
    network_analysis = NetworkAnalysis()
    network_analysis.set_analysis_type(analysis_type='inter_tissue')

    print("=== Finding Common Genes Across Tissues ===")

    # Get all gene sets from each tissue (without tissue prefixes)
    tissue_gene_sets = {}
    for tissue, tissue_data in tissue_processed_data.items():
        # Remove tissue prefix to get just the gene IDs
        genes_without_prefix = set([gene.split('_', 1)[1] for gene in tissue_data.index])
        tissue_gene_sets[tissue] = genes_without_prefix
        rprint(f"{tissue}: {len(genes_without_prefix)} genes")

    # Find intersection of all tissue gene sets
    common_genes = set.intersection(*tissue_gene_sets.values())
    common_genes = list(common_genes)
    rprint(f"Common genes across all tissues: {len(common_genes)}")

    print("=== Preparing Data for Inter-Tissue Network (R Script Format) ===")

    tissue_expression_data = {}
    tissue_sample_mapping = {}

    # Note: You'll need to load these from your preprocessing
    # For now, creating dummy sample mappings
    tissue_young_samples = {
        'Brain - Cortex': [],
        'Muscle - Skeletal': [], 
        'Adipose - Subcutaneous': []
    }
    tissue_old_samples = {
        'Brain - Cortex': [],
        'Muscle - Skeletal': [], 
        'Adipose - Subcutaneous': []
    }

    for tissue in target_tissues:
        rprint(f"Processing tissue: {tissue}")
        
        tissue_data = tissue_processed_data[tissue]
        # Filter to common genes (keeping tissue prefix in the data)
        common_tissue_genes = [f"{tissue}_{gene}" for gene in common_genes if f"{tissue}_{gene}" in tissue_data.index]
        
        # CRITICAL: Keep data as genes x samples (R script format)
        # The inter_tissue_network_analysis.py will handle transposition internally
        tissue_subset = tissue_data.loc[common_tissue_genes]
        rprint(f"  Tissue data shape: {tissue_subset.shape} (genes × samples)")
        
        # Store in the original format
        tissue_expression_data[tissue] = tissue_subset
        tissue_sample_mapping[tissue] = tissue_young_samples[tissue] + tissue_old_samples[tissue]

    print("=== Loading Data into Network Analysis ===")

    network_analysis.load_tissue_data_for_inter_tissue_analysis(
        tissue_expression_data,
        tissue_sample_mapping
    )

    print("=== Constructing Inter-Tissue Network ===")

    inter_tissue_results = network_analysis.construct_inter_tissue_network(
        ts_power=6,
        ct_power=3,
        correlation_method='pearson',
        min_module_size=30
    )
    
    print("=== Results ===")
    print(f"Adjacency matrix shape: {inter_tissue_results['adjacency_matrix'].shape}")
    print(f"TOM matrix shape: {inter_tissue_results['tom_matrix'].shape}")
    print(f"Number of modules: {len(set(inter_tissue_results['module_assignments']))}")

if __name__ == "__main__":
    main()
