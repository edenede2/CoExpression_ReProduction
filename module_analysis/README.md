# WGCNA and SpeakEasy2 Module Comparison Tool

This tool analyzes and compares gene co-expression modules identified by WGCNA and SpeakEasy2 clustering algorithms. It calculates overlap metrics, generates visualizations, and provides detailed gene lists for the top intersections.

## Features

- Processes WGCNA module tables and SpeakEasy2 cluster tables
- Calculates module sizes, intersection sizes, and Jaccard indices
- Handles tissue-specific gene identifiers
- Generates bar plots showing top module pairs
- Creates detailed gene lists for the top intersections
- Supports filtering by module/community size

## Usage

```bash
python analyze_wgcna_se2_modules.py --wgcna-file <path_to_wgcna_file> --se2-file <path_to_se2_file> --outdir <output_directory> [options]
```

### Required Arguments

- `--wgcna-file`: Path to the WGCNA module table TSV file
- `--se2-file`: Path to the SpeakEasy2 clusters table CSV file
- `--outdir`: Directory to save output files

### Optional Arguments

- `--top-n`: Number of top pairs to show in the overall bar plot (default: 150)
- `--top-modules`: Number of top modules to show in the per-module bar plot (default: 20)
- `--k`: Number of top communities to show per module (default: 3)
- `--min-module-size`: Minimum WGCNA module size to include (default: 0)
- `--min-community-size`: Minimum SE2 community size to include (default: 0)

## Output Files

1. `module_pair_metrics.csv`: Detailed metrics for all module-community pairs
   - wgcna_module, wgcna_tissue: WGCNA module ID and tissue
   - se2_community, se2_tissue: SE2 community ID and tissue
   - wgcna_size, se2_size: Number of genes in each module/community
   - intersection_size, union_size: Size of intersection and union gene sets
   - jaccard_index: Jaccard similarity coefficient (intersection/union)
   - same_tissue: Boolean flag indicating if both are from the same tissue
   - pair_id: Combined identifier for the module-community pair

2. `top_intersections_gene_lists.csv`: Gene lists for top 30 intersection pairs
   - intersection_genes: Genes found in both modules
   - wgcna_only_genes: Genes unique to the WGCNA module
   - se2_only_genes: Genes unique to the SE2 community

3. `bar_top_pairs.png`: Bar plot of top N pairs by Jaccard index

4. `bar_module_topk.png`: Bar plot showing top K communities for each module

## Example

```bash
python analyze_wgcna_se2_modules.py \
  --wgcna-file /path/to/wgcna_modules.tsv \
  --se2-file /path/to/se2_clusters.csv \
  --outdir ./results \
  --top-n 100 \
  --min-module-size 10
```

## File Format Requirements

### WGCNA Module Table Format (TSV)
```
Cluster ID  Tissue  Gene Symbol     stringAsFactors
1     AC      ENSG00000001497.16      FALSE
1     AC      ENSG00000002016.17      FALSE
...
```

### SpeakEasy2 Cluster Table Format (CSV)
```
Resolution,ClusterColumn,Cluster.ID,Tissue,Gene.Symbol,Gene.ID
3,cluster3_id,36,AC,ENSG00000000003.14,ENSG00000000003.14
3,cluster3_id,49,AC,ENSG00000000419.12,ENSG00000000419.12
...
```

## Note on Tissue-Specific Analysis

This tool handles tissue-specific gene identifiers and provides a comprehensive analysis of module overlaps, taking into account tissue origin. The `same_tissue` flag in the output helps identify module pairs from the same or different tissues.