# Final Heatmap Updates - Correlation Visualization

## Summary of Latest Changes

The `phenotype_heatmap_with_pathway_legend.py` script has been updated to show **signed correlation values** with a **diverging colormap** and **significance stars**.

### New Features

#### 1. Diverging Colormap (Blue-Red)
- **Red colors** = Positive correlations
- **Blue colors** = Negative correlations  
- **White** = Zero correlation
- Uses matplotlib's `RdBu_r` (Red-Blue reversed) colormap
- Value range: -1.0 to +1.0 for correlations

#### 2. Significance Stars
Three levels of statistical significance are shown with stars overlaid on the heatmap:
- `***` = p_adj ≤ 0.001 (highly significant)
- `**`  = p_adj ≤ 0.01 (very significant)
- `*`   = p_adj ≤ 0.05 (significant)

Stars are displayed with:
- Black color with bold font
- White outline for visibility on colored backgrounds
- Centered on each heatmap cell

#### 3. Signed Correlation Handling
- The `metric` parameter default changed from `"neglog10_padj"` to `"corr"`
- Correlations keep their sign (positive/negative)
- For duplicate module-phenotype pairs, the value with maximum absolute value is selected
- `cap` parameter default changed from `6.0` to `1.0` to match correlation range

#### 4. Module Selection by Significance
- When `top_modules` is specified, modules are now ranked by **minimum p-value** (not correlation)
- Selects the modules with the strongest statistical evidence across phenotypes

## Visual Example

```
Heatmap Layout:
┌────┬────┬─────────────────────────────────┬────┬─────────────┐
│Sub │Cat │   Correlation Heatmap           │    │  Legend:    │
│cat │    │   *** ** *                      │    │  Categories │
│Col │Col │   🔵 = Negative (anticorr)      │    │  &          │
│ors │ors │   ⚪ = Zero                     │    │  Subcats    │
│    │    │   🔴 = Positive (correlated)    │    │             │
└────┴────┴─────────────────────────────────┴────┴─────────────┘
```

## Usage

### Default (Signed Correlations with Stars):
```python
from phenotype_heatmap_with_pathway_legend import run_all

run_all(
    pheno_tsv="rosmap_ME_pheno_correlations.tsv",
    kegg_csv="kegg_results.csv",
    details_tsv="cluster_details.tsv",
    out_prefix="my_analysis",
    tissues=['CROSS'],
    exclude_phenotypes=[r"msex"],
    metric="corr",      # Signed correlations
    cap=1.0,            # Correlation range [-1, 1]
    figure_dpi=180,
)
```

### Select Top 100 Most Significant Modules:
```python
run_all(
    pheno_tsv="rosmap_ME_pheno_correlations.tsv",
    kegg_csv="kegg_results.csv",
    details_tsv="cluster_details.tsv",
    out_prefix="my_analysis_top100",
    metric="corr",
    cap=1.0,
    top_modules=100,    # Select 100 modules with lowest p-values
    figure_dpi=180,
)
```

### Use Absolute Correlations (Old Style):
```python
run_all(
    pheno_tsv="rosmap_ME_pheno_correlations.tsv",
    kegg_csv="kegg_results.csv",
    details_tsv="cluster_details.tsv",
    out_prefix="my_analysis_abs",
    metric="neglog10_padj",  # Use -log10(p_adj) instead
    cap=6.0,
    figure_dpi=180,
)
```

## Technical Details

### Colormap
- **Colormap**: `RdBu_r` (reversed Red-Blue)
- **Range**: -cap to +cap (default: -1.0 to 1.0)
- **Center**: 0.0 (white)
- **Interpretation**: 
  - Values near -1.0 → Dark blue (strong negative correlation)
  - Values near 0.0 → White (no correlation)
  - Values near +1.0 → Dark red (strong positive correlation)

### Significance Testing
- **Source**: Uses `p_adj` column from phenotype TSV file
- **Thresholds**: Can be customized in the `alpha_levels` parameter
- **Display**: Stars with white stroke outline for visibility
- **Font**: Size 8, bold weight

### Data Processing
- NaN values are filtered out before processing
- For multiple observations per module-phenotype pair:
  - Correlation matrix: Takes value with maximum absolute correlation
  - P-value matrix: Takes minimum p-value
- Module ranking uses minimum p-value across all phenotypes

## Output Files

All output files remain the same:
1. **`{prefix}_heatmap.pdf`** - Main visualization (now with diverging colors + stars)
2. **`{prefix}_module_details_selected.tsv`** - Module details
3. **`{prefix}_top_phenotypes_per_module.tsv`** - Top phenotypes per module
4. **`{prefix}_top_pathways_per_module.tsv`** - Top pathways per module

## Parameters Updated

| Parameter | Old Default | New Default | Description |
|-----------|-------------|-------------|-------------|
| `metric` | `"neglog10_padj"` | `"corr"` | Now shows signed correlations |
| `cap` | `6.0` | `1.0` | Adjusted for correlation range |

## File Size

The new PDF with stars is larger (~401KB vs ~280KB) due to the additional text elements (stars) rendered on the heatmap.
