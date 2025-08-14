# Inter-Tissue Network Analysis Fix Summary

## Problem Identified

The original error was:
```
ValueError: could not broadcast input array from shape (2444,2444) into shape (242,242)
```

This occurred because the inter-tissue network analysis was incorrectly handling the data format.

## Root Cause

1. **Data Format Confusion**: The notebook was transposing data to `samples × genes`, but the inter-tissue analysis expected `genes × samples`
2. **Tissue Gene Naming**: The R script creates tissue-specific gene names (e.g., `Brain_GENE1`, `Muscle_GENE1`) to treat the same gene from different tissues as separate entities
3. **Adjacency Matrix Structure**: The adjacency matrix should be `(total_genes_all_tissues × total_genes_all_tissues)`, not `(samples × samples)`

## Key Fixes Applied

### 1. Fixed `inter_tissue_network_analysis.py`

**Updated `load_tissue_data()` method:**
- Now properly handles both `genes × samples` and `samples × genes` input formats
- Automatically transposes to the correct internal format (`samples × genes`)
- Creates tissue-specific gene names following R script convention: `{tissue}_{gene}`
- Properly tracks gene positions for adjacency matrix construction

**Updated `construct_inter_tissue_adjacency()` method:**
- Now follows the R script `AdjacencyFromExpr` function exactly
- Correctly implements the nested loop structure from the R script
- Properly handles tissue-specific (TS) and cross-tissue (CT) correlations
- Uses the correct matrix indexing based on `tissue_indices`

### 2. Corrected Data Preparation

**In the notebook/script:**
- Keep preprocessed data in `genes × samples` format (original preprocessing output)
- Don't transpose to `samples × genes` - let `inter_tissue_network_analysis.py` handle format internally
- Filter to common genes while preserving tissue prefixes
- Pass tissue sample mappings correctly

### 3. R Script Alignment

**Key R script patterns now implemented:**
```r
# R script creates tissue-specific gene names
colnames(datExpr) <- paste(tissue_name,"_",colnames(datExpr),sep="")

# R script adjacency matrix structure
adj_mat<-matrix(0,nrow=tissue_index_adj[total_tissues+1], ncol=tissue_index_adj[total_tissues+1])

# R script tissue-specific correlations
adj_mat[(tissue_index_adj[i]+1):tissue_index_adj[i+1], (tissue_index_adj[i]+1):tissue_index_adj[i+1]] <- 
  abs(cor(vector_expr[[i]], method = cor_method))^TS_power

# R script cross-tissue correlations  
adj_mat[(tissue_index_adj[i]+1):tissue_index_adj[i+1],(tissue_index_adj[j]+1):tissue_index_adj[j+1]] <- 
  abs(cor(vector_expr[[i]][common_Samples,],vector_expr[[j]][common_Samples,], method = cor_method))^CT_power
```

## Expected Results

After these fixes:
- **Adjacency Matrix**: `(total_genes_all_tissues × total_genes_all_tissues)` 
  - For 3 tissues with 1687 common genes each: `(5061 × 5061)`
- **Gene Names**: `["Brain - Cortex_GENE1", "Brain - Cortex_GENE2", ..., "Muscle - Skeletal_GENE1", ...]`
- **Matrix Structure**: 
  - Diagonal blocks: tissue-specific correlations
  - Off-diagonal blocks: cross-tissue correlations
- **Powers**: TS power (typically 6) for within-tissue, CT power (typically 3) for between-tissue

## Usage

1. Run preprocessing to get `tissue_processed_data.pkl`
2. Use the corrected format:
   ```python
   # Keep as genes × samples (don't transpose)
   tissue_expression_data[tissue] = tissue_data.loc[common_tissue_genes]
   
   # Load into network analysis
   network_analysis.load_tissue_data_for_inter_tissue_analysis(
       tissue_expression_data, tissue_sample_mapping
   )
   ```

3. The inter-tissue analysis will now correctly:
   - Create tissue-specific gene identifiers
   - Build the proper adjacency matrix structure
   - Apply correct correlation calculations
   - Generate meaningful inter-tissue modules
