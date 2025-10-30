# Phenotype Heatmap Updates

## Summary of Changes

The `phenotype_heatmap_with_pathway_legend.py` script has been updated with the following new features:

### 1. PDF Output (instead of PNG)
- The heatmap is now saved as a **PDF file** for publication-quality vector graphics
- Output filename: `{out_prefix}_heatmap.pdf`

### 2. Category and Subcategory Annotations
- **Two color strips** are now displayed on the **left side** of the heatmap matrix
- These strips show the dominant KEGG pathway category and subcategory for each module
- Color strips are synchronized with the heatmap rows (one color per module)
- Colors match the legend on the right side for easy interpretation
- Style matches the `enrichmentViz.py` implementation

### 3. Module Selection by Phenotype Significance
- New parameter: **`top_modules`** (optional integer)
- When set, selects only the top N modules ranked by their maximum significance across all phenotypes
- For metrics like `neglog10_padj`, higher values = more significant
- Example: `top_modules=100` will show only the 100 most phenotype-associated modules

## Usage Examples

### Basic usage (all modules):
```python
from phenotype_heatmap_with_pathway_legend import run_all

run_all(
    pheno_tsv="rosmap_ME_pheno_correlations.tsv",
    kegg_csv="kegg_results.csv",
    details_tsv="cluster_details.tsv",
    out_prefix="my_analysis",
    tissues=['CROSS'],
    exclude_phenotypes=[r"msex"],
    metric="neglog10_padj",
    cap=6.0,
    figure_dpi=180,
)
```

### Select top 100 modules by phenotype significance:
```python
run_all(
    pheno_tsv="rosmap_ME_pheno_correlations.tsv",
    kegg_csv="kegg_results.csv",
    details_tsv="cluster_details.tsv",
    out_prefix="my_analysis_top100",
    tissues=['CROSS'],
    exclude_phenotypes=[r"msex"],
    metric="neglog10_padj",
    cap=6.0,
    top_modules=100,  # <-- Select top 100 modules
    figure_dpi=180,
)
```

## Output Files

1. **`{prefix}_heatmap.pdf`** - Main heatmap visualization with:
   - Subcategory color strip (leftmost)
   - Category color strip (second from left)
   - Heatmap matrix (modules × phenotypes)
   - Legend (right side) showing category/subcategory colors

2. **`{prefix}_module_details_selected.tsv`** - Details for modules in heatmap
3. **`{prefix}_top_phenotypes_per_module.tsv`** - Top K phenotypes per module
4. **`{prefix}_top_pathways_per_module.tsv`** - Top K pathways per module

## Visual Layout

```
┌────┬────┬─────────────────────┬────┬─────────────┐
│Sub │Cat │                     │    │             │
│cat │    │                     │    │  Legend:    │
│    │    │     Heatmap         │    │  Categories │
│Col │Col │     Matrix          │    │  &          │
│ors │ors │                     │    │  Subcats    │
│    │    │                     │    │             │
└────┴────┴─────────────────────┴────┴─────────────┘
```

## Parameters Reference

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `pheno_tsv` | str | required | Path to phenotype correlation file |
| `kegg_csv` | str | required | Path to KEGG enrichment results |
| `details_tsv` | str | required | Path to module details file |
| `out_prefix` | str | "heatmap" | Output filename prefix |
| `tissues` | list | None | Filter to specific tissues |
| `include_phenotypes` | list | None | Regex patterns to include phenotypes |
| `exclude_phenotypes` | list | None | Regex patterns to exclude phenotypes |
| `metric` | str | "neglog10_padj" | Metric to plot (neglog10_padj, neglog10_q, neglog10_p, corr, beta) |
| `cap` | float | 6.0 | Maximum value for color scale |
| **`top_modules`** | int | None | **NEW: Select top N modules by phenotype significance** |
| `top_k_phenos` | int | 10 | Number of top phenotypes per module in TSV output |
| `top_k_pathways` | int | 5 | Number of top pathways per module in TSV output |
| `figure_dpi` | int | 160 | Figure resolution |
| `figsize_scale` | float | 1.0 | Scale factor for figure size |

## Technical Details

- Module annotations (category/subcategory) are extracted from the KEGG enrichment file
- For each module, the category/subcategory of the **most significantly enriched pathway** is used
- Missing annotations are shown in light gray
- Color palettes use matplotlib's tab20, tab20b, tab20c, and Set3 colormaps for variety
