# CORRECTED CODE FOR INTER-TISSUE NETWORK ANALYSIS
# Copy this into your notebook to replace the problematic cell

from new_script.network_analysis import NetworkAnalysis

network_analysis = NetworkAnalysis()
network_analysis.set_analysis_type(analysis_type='inter_tissue')

print("=== Finding Common Genes Across Tissues ===")

# Get all gene sets from each tissue (without tissue prefixes) - CORRECTED LOGIC
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

print("=== Preparing Data for Inter-Tissue Network ===")

tissue_expression_data = {}
tissue_sample_mapping = {}

for tissue in target_tissues:
    rprint(f"Processing tissue: {tissue}")
    
    tissue_data = tissue_processed_data[tissue]
    common_tissue_genes = [f"{tissue}_{gene}" for gene in common_genes if f"{tissue}_{gene}" in tissue_data.index]
    
    # Get the subset of data with common genes
    tissue_subset = tissue_data.loc[common_tissue_genes]
    rprint(f"  Before transpose: {tissue_subset.shape} (genes × samples)")
    
    # CRITICAL FIX: Transpose to samples × genes format for inter-tissue analysis
    # The inter-tissue network code expects: rows = samples, columns = genes
    tissue_expression_data[tissue] = tissue_subset.T  # Transpose: genes×samples -> samples×genes
    
    # Remove the tissue prefix from gene names for the columns
    tissue_expression_data[tissue].columns = [gene.split('_', 1)[1] for gene in tissue_expression_data[tissue].columns]

    tissue_sample_mapping[tissue] = tissue_young_samples[tissue] + tissue_old_samples[tissue]
    
    rprint(f"  After transpose: {tissue_expression_data[tissue].shape} (samples × genes)")

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
