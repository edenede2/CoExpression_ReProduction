# Fix for inter-tissue network data format issue
# The problem: inter-tissue network expects samples×genes but we have genes×samples

import pickle
import pandas as pd
from rich import print as rprint

# Load the saved processed data
with open("tissue_processed_data.pkl", "rb") as f:
    tissue_processed_data = pickle.load(f)

target_tissues = ['Brain - Cortex', 'Muscle - Skeletal', 'Adipose - Subcutaneous']

# Debug: Check current data format
print("=== Current Data Format ===")
for tissue in target_tissues:
    tissue_data = tissue_processed_data[tissue]
    rprint(f"{tissue} shape: {tissue_data.shape}")
    rprint(f"  Index (genes): {tissue_data.index[:3].tolist()}...")
    rprint(f"  Columns (samples): {tissue_data.columns[:3].tolist()}...")
    print()

# The fix: when preparing data for inter-tissue network, transpose it
print("=== Preparing Data for Inter-Tissue Network ===")

# Load common genes (from previous processing) - FIX THE LOGIC
print("=== Finding Common Genes ===")

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

# Prepare correctly formatted data
tissue_expression_data = {}

for tissue in target_tissues:
    rprint(f"Processing {tissue}...")
    
    tissue_data = tissue_processed_data[tissue]
    common_tissue_genes = [f"{tissue}_{gene}" for gene in common_genes if f"{tissue}_{gene}" in tissue_data.index]
    
    # Get the subset of data with common genes
    tissue_subset = tissue_data.loc[common_tissue_genes]
    rprint(f"  Before transpose: {tissue_subset.shape} (genes × samples)")
    
    # CRITICAL FIX: Transpose to samples × genes format
    tissue_expression_data[tissue] = tissue_subset.T  # genes×samples -> samples×genes
    
    # Remove the tissue prefix from gene names for the columns
    tissue_expression_data[tissue].columns = [gene.split('_', 1)[1] for gene in tissue_expression_data[tissue].columns]
    
    rprint(f"  After transpose: {tissue_expression_data[tissue].shape} (samples × genes)")
    rprint(f"  Gene names (columns): {tissue_expression_data[tissue].columns[:3].tolist()}...")
    rprint(f"  Sample names (index): {tissue_expression_data[tissue].index[:3].tolist()}...")
    print()

print("=== Data Format Fixed! ===")
print("Now the data is in the correct format for inter-tissue network analysis:")
print("- Rows = samples")
print("- Columns = genes")
print("- Gene names are the same across tissues (no tissue prefix)")
