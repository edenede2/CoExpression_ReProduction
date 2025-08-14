# CoExpression Network Analysis Workflow Guide

This guide provides detailed steps for running the new Python-based CoExpression Network Analysis pipeline, based on the methodology described in the research article and implemented in the original R scripts.

## 🎯 Overview

The pipeline implements two main analysis approaches:
1. **Traditional WGCNA**: Single-condition coexpression network analysis  
2. **Inter-tissue X-WGCNA**: Cross-tissue network analysis for studying tissue coordination

Both approaches can be followed by **Modular Differential Connectivity (MDC)** analysis to compare connectivity patterns between different conditions (e.g., young vs old).

---

## 📂 Prerequisites and Setup

### 1. Data Requirements
Ensure the following files are present in the `data/` directory:
- `GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct` - GTEx expression data
- `GTEx_Analysis_v8_Annotations_SampleAttributesDS.tsv` - Sample metadata
- `GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.tsv` - Subject phenotypes
- `hgnc_complete_set.tsv` - Protein coding gene annotations

### 2. Environment Setup
```bash
cd new_script
pip install -r requirements.txt
```

---

## 🔄 WORKFLOW 1: Traditional WGCNA Analysis (Single Tissue)

This workflow analyzes coexpression within a single tissue type, comparing different conditions (e.g., young vs old).

### Step 1: Data Loading and Initial Filtering

**Function Order:**
1. `GTExDataReader()` - Initialize data reader
2. `read_gct_file()` - Load expression data
3. `read_sample_attributes()` - Load sample metadata  
4. `read_subject_phenotypes()` - Load subject data
5. `filter_samples_by_metadata()` - Filter samples by tissue and condition

**Arguments and Usage:**
```python
from data_reader import GTExDataReader

# Initialize reader
reader = GTExDataReader()

# Load expression data (optionally limit samples for testing)
expression_df, gene_metadata = reader.read_gct_file(
    "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct",
    sample_limit=None  # Use None for full dataset, or number for subset
)

# Load metadata
sample_attrs = reader.read_sample_attributes()
subject_phenos = reader.read_subject_phenotypes()

# Filter samples for specific tissue and age groups
young_samples = reader.filter_samples_by_metadata(
    sample_attrs, 
    subject_phenos,
    tissue='Brain - Cortex',           # Target tissue
    age_group='young',                 # 'young' (20-59 years) or 'old' (60-79 years)
    sex=None,                         # Optional: 'Male', 'Female'  
    min_samples=30                    # Minimum samples required
)

old_samples = reader.filter_samples_by_metadata(
    sample_attrs,
    subject_phenos, 
    tissue='Brain - Cortex',
    age_group='old',
    sex=None,
    min_samples=30
)
```

**Details:** This step loads the GTEx data and filters samples based on tissue type, age group, and optionally sex. The age groups are defined as young (20-59 years) and old (60-79 years) to match the original study methodology.

### Step 2: Data Preprocessing and Quality Control

**Function Order:**
1. `DataPreprocessor()` - Initialize preprocessor
2. `filter_protein_coding_genes()` - Keep only protein coding genes
3. `preprocess_expression_data()` - Comprehensive preprocessing including:
   - Log transformation
   - Variance filtering
   - RIN score filtering
   - PCA-based outlier detection
   - Sample quality control

**Arguments and Usage:**
```python
from preprocessing import DataPreprocessor

# Initialize preprocessor
preprocessor = DataPreprocessor()

# Combine samples for preprocessing
all_samples = young_samples + old_samples
expression_subset = expression_df[all_samples]

# Load HGNC data for protein coding gene filtering
hgnc_data = reader.read_hgnc_data()

# Filter to protein coding genes only
expression_subset = preprocessor.filter_protein_coding_genes(
    expression_subset, 
    hgnc_data
)

# Comprehensive preprocessing
processed_df, preprocessing_info = preprocessor.preprocess_expression_data(
    expression_subset,
    log_transform=True,              # Apply log2(TPM + 1) transformation
    variance_filter=True,            # Filter low-variance genes
    variance_threshold=0.1,          # Minimum variance threshold
    rin_filter=True,                # Filter samples by RNA Integrity Number
    min_rin=6.0,                    # Minimum RIN score
    outlier_detection=True,          # PCA-based outlier detection
    outlier_threshold=2.5,          # Standard deviations for outlier cutoff
    min_samples_per_group=15        # Minimum samples required per group
)
```

**Details:** This step performs comprehensive quality control similar to the original R pipeline. It filters genes and samples based on variance, RNA quality, and statistical outlier detection. The preprocessing ensures data quality for downstream network analysis.

### Step 3: Network Construction and Module Detection

**Function Order:**
1. `NetworkAnalysis()` - Initialize network analyzer
2. `set_analysis_type('traditional')` - Set to traditional WGCNA mode
3. `select_soft_threshold()` - Determine optimal power parameter
4. `construct_network()` - Build coexpression network and detect modules

**Arguments and Usage:**
```python
from network_analysis import NetworkAnalysis

# Initialize network analyzer  
network_analyzer = NetworkAnalysis()
network_analyzer.set_analysis_type('traditional')

# Select soft threshold power (equivalent to R's pickSoftThreshold)
power_analysis = network_analyzer.select_soft_threshold(
    processed_df,
    powers=range(1, 21),           # Test powers from 1 to 20
    network_type='unsigned',       # 'unsigned' or 'signed'
    correlation_method='pearson'   # 'pearson' or 'spearman'
)

# Use recommended power (typically 6-12 for unsigned networks)
optimal_power = power_analysis['recommended_power']

# Construct network and detect modules
network_results = network_analyzer.construct_network(
    processed_df,
    power=optimal_power,           # Soft threshold power
    network_type='unsigned',       # Network type
    tom_type='unsigned',          # TOM calculation type
    min_module_size=30,           # Minimum genes per module
    deepSplit=2,                  # Module detection sensitivity (0-3)
    merge_threshold=0.25,         # Module merging threshold
    correlation_method='pearson'   # Correlation method
)
```

**Details:** This step constructs the weighted gene coexpression network using the WGCNA methodology. It calculates adjacency matrices, topological overlap matrices (TOM), and performs hierarchical clustering to identify gene modules. The soft threshold power is selected to achieve scale-free topology.

### Step 4: Modular Differential Connectivity (MDC) Analysis

**Function Order:**
1. `MDCAnalyzer()` - Initialize MDC analyzer
2. `analyze_differential_connectivity()` - Compare connectivity between conditions
3. `perform_permutation_test()` - Statistical significance testing
4. `apply_fdr_correction()` - Multiple testing correction

**Arguments and Usage:**
```python
from mdc_analysis import MDCAnalyzer

# Initialize MDC analyzer
mdc_analyzer = MDCAnalyzer()

# Split processed data by condition
young_data = processed_df[young_samples]
old_data = processed_df[old_samples]

# Perform MDC analysis
mdc_results = mdc_analyzer.analyze_differential_connectivity(
    young_data,                    # Condition 1 data
    old_data,                     # Condition 2 data  
    network_results['modules'],    # Module assignments
    optimal_power,                # Network power parameter
    "Young",                      # Condition 1 name
    "Old",                       # Condition 2 name
    permutations=1000,           # Number of permutation tests
    random_seed=42,              # For reproducibility
    min_module_size=10,          # Minimum module size for analysis
    connectivity_type='kWithin'   # 'kWithin', 'kTotal', or 'kOut'
)

# Apply FDR correction
mdc_results_corrected = mdc_analyzer.apply_fdr_correction(
    mdc_results,
    method='fdr_bh'              # Benjamini-Hochberg FDR correction
)
```

**Details:** MDC analysis compares the intra-module connectivity between young and old conditions. It calculates the ratio of connectivity and uses permutation testing to assess statistical significance. FDR correction accounts for multiple testing across modules.

---

## 🔄 WORKFLOW 2: Inter-tissue X-WGCNA Analysis

This workflow studies coordination patterns across multiple tissues, implementing the X-WGCNA methodology for cross-tissue analysis.

### Step 1: Multi-tissue Data Loading

**Function Order:**
1. `GTExDataReader()` - Initialize data reader
2. Load data for multiple tissues
3. `filter_samples_by_metadata()` - Filter samples for each tissue
4. Ensure sample overlap across tissues

**Arguments and Usage:**
```python
from data_reader import GTExDataReader

reader = GTExDataReader()

# Define tissues to analyze
target_tissues = [
    'Brain - Cortex',
    'Adipose - Subcutaneous',
    'Muscle - Skeletal'
]

# Load expression data
expression_df, gene_metadata = reader.read_gct_file(
    "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct"
)

# Load metadata
sample_attrs = reader.read_sample_attributes()
subject_phenos = reader.read_subject_phenotypes()

# Filter samples for each tissue and condition
tissue_young_samples = {}
tissue_old_samples = {}

for tissue in target_tissues:
    tissue_young_samples[tissue] = reader.filter_samples_by_metadata(
        sample_attrs, subject_phenos,
        tissue=tissue,
        age_group='young',
        min_samples=20
    )
    
    tissue_old_samples[tissue] = reader.filter_samples_by_metadata(
        sample_attrs, subject_phenos,
        tissue=tissue, 
        age_group='old',
        min_samples=20
    )
```

### Step 2: Cross-tissue Data Preprocessing  

**Function Order:**
1. `DataPreprocessor()` - Initialize preprocessor
2. Preprocess each tissue separately
3. Identify common genes across tissues
4. Create tissue-prefixed gene identifiers

**Arguments and Usage:**
```python
from preprocessing import DataPreprocessor

preprocessor = DataPreprocessor()

# Preprocess each tissue separately
tissue_processed_data = {}

for tissue in target_tissues:
    # Combine young and old samples for this tissue
    tissue_samples = tissue_young_samples[tissue] + tissue_old_samples[tissue]
    tissue_expression_df = expression_df[tissue_samples]

    # Preprocess tissue expression data
    processed_tissue, info = preprocessor.preprocess_expression_data(
        expression_df=tissue_expression_df,
        hgnc_data=hgnc_data,
        sample_attrs=sample_attrs,
        subject_phenos=subject_phenos,
        apply_log=True,
        filter_low_expr=True,
        filter_low_var=True,
        detect_outliers=True,
        quantile_norm= True,
        regress_confounders= True,
        min_expression=0.1,
        min_variance=0.2,
        outlier_contamination=0.01
    )
    
    # Add tissue prefix to gene names
    processed_tissue.index = [f"{tissue}_{gene}" for gene in processed_tissue.index]
    tissue_processed_data[tissue] = processed_tissue

# Identify common genes across tissues (before prefixing)
common_genes = set(expression_df.index)
for tissue_data in tissue_processed_data.values():
    tissue_genes = set([gene.split('_', 1)[1] for gene in tissue_data.index])
    common_genes = common_genes.intersection(tissue_genes)

print(f"Common genes across tissues: {len(common_genes)}")
```

### Step 3: Inter-tissue Network Construction

**Function Order:**
1. `NetworkAnalysis()` - Initialize network analyzer
2. `set_analysis_type('inter_tissue')` - Set to inter-tissue mode
3. `load_tissue_data_for_inter_tissue_analysis()` - Load multi-tissue data
4. `construct_inter_tissue_network()` - Build cross-tissue network

**Arguments and Usage:**
```python
from network_analysis import NetworkAnalysis

# Initialize network analyzer for inter-tissue analysis
network_analyzer = NetworkAnalysis()
network_analyzer.set_analysis_type('inter_tissue')

# Prepare tissue data for analysis
tissue_expression_data = {}
tissue_sample_mapping = {}

for tissue in target_tissues:
    # Filter to common genes only
    tissue_data = tissue_processed_data[tissue]
    common_tissue_genes = [f"{tissue}_{gene}" for gene in common_genes 
                          if f"{tissue}_{gene}" in tissue_data.index]
    tissue_expression_data[tissue] = tissue_data.loc[common_tissue_genes]
    
    # Map samples to subjects for finding overlaps
    tissue_sample_mapping[tissue] = tissue_young_samples[tissue] + tissue_old_samples[tissue]

# Load data into analyzer
network_analyzer.load_tissue_data_for_inter_tissue_analysis(
    tissue_expression_data,
    tissue_sample_mapping
)

# Construct inter-tissue network
inter_tissue_results = network_analyzer.construct_inter_tissue_network(
    ts_power=6,                   # Within-tissue power
    ct_power=3,                   # Cross-tissue power  
    correlation_method='pearson', # Correlation method
    min_module_size=30,          # Minimum module size
    deepSplit=2,                 # Module detection sensitivity
    merge_threshold=0.25         # Module merging threshold
)
```

**Details:** This step implements X-WGCNA methodology, constructing adjacency matrices with different powers for within-tissue (ts_power) and cross-tissue (ct_power) connections. The resulting network captures both intra-tissue coexpression and inter-tissue coordination patterns.

### Step 4: Inter-tissue MDC Analysis

**Function Order:**
1. `InterTissueMDCAnalysis()` - Initialize inter-tissue MDC analyzer
2. `load_modules()` - Load cross-tissue module assignments
3. `load_condition_networks()` - Load networks for different conditions
4. `compute_inter_tissue_mdc()` - Calculate cross-tissue MDC scores

**Arguments and Usage:**
```python
from inter_tissue_mdc_analysis import InterTissueMDCAnalysis

# Initialize inter-tissue MDC analyzer
inter_mdc_analyzer = InterTissueMDCAnalysis()

# Load module assignments from inter-tissue analysis
inter_mdc_analyzer.load_modules(inter_tissue_results['modules'])

# Build condition-specific networks
young_networks = {}
old_networks = {}

for condition, samples_dict in [('young', tissue_young_samples), ('old', tissue_old_samples)]:
    # Prepare condition-specific data
    condition_tissue_data = {}
    
    for tissue in target_tissues:
        tissue_samples = samples_dict[tissue]
        tissue_expr = tissue_processed_data[tissue][tissue_samples]
        condition_tissue_data[tissue] = tissue_expr
    
    # Build network for this condition
    condition_network = network_analyzer.construct_inter_tissue_network(
        condition_tissue_data,
        ts_power=6,
        ct_power=3,
        correlation_method='pearson'
    )
    
    if condition == 'young':
        young_networks = condition_network
    else:
        old_networks = condition_network

# Load networks into MDC analyzer  
inter_mdc_analyzer.load_condition_networks({
    'Young': young_networks['adjacency_matrix'],
    'Old': old_networks['adjacency_matrix']
})

# Compute inter-tissue MDC
inter_mdc_results = inter_mdc_analyzer.compute_inter_tissue_mdc(
    condition1='Young',
    condition2='Old', 
    min_module_size=10,
    permutations=1000,
    connectivity_type='intra_module'
)
```

**Details:** Inter-tissue MDC analysis compares cross-tissue module connectivity between conditions. This reveals whether coordination patterns between tissues change with aging or other conditions.

---

## 📊 Results Interpretation and Visualization

### Visualization Functions

**Available Visualization Functions:**
```python
from visualization import NetworkVisualizer

visualizer = NetworkVisualizer()

# Module connectivity heatmaps
visualizer.plot_module_connectivity_heatmap(
    mdc_results,
    title="Module Connectivity: Young vs Old"
)

# MDC bar plots
visualizer.plot_mdc_barplot(
    mdc_results_corrected,
    significance_threshold=0.05,
    title="Differential Connectivity Analysis"
)

# Network module visualization
visualizer.plot_network_modules(
    network_results['adjacency_matrix'],
    network_results['modules'],
    layout='spring'
)

# Expression heatmaps for significant modules
visualizer.plot_module_expression_heatmap(
    processed_df,
    network_results['modules'],
    module_id=1,  # Specific module to visualize
    conditions=['Young', 'Old'],
    sample_labels=young_samples + old_samples
)
```

---

## 🚀 Running the Complete Pipeline

### Option 1: Streamlit Web Interface (Recommended)

```bash
cd new_script
streamlit run streamlit_app.py
```

**Workflow in Streamlit:**
1. **Data Loading Page**: Load and explore GTEx data
2. **Data Preprocessing Page**: Configure and run preprocessing
3. **Network Analysis Page**: Choose analysis type and build networks
4. **MDC Analysis Page**: Compare conditions and run statistical tests
5. **Results Visualization Page**: Generate plots and export results

### Option 2: Command Line Execution

Create a custom script following the workflows above:

```python
# example_pipeline.py
from data_reader import GTExDataReader
from preprocessing import DataPreprocessor  
from network_analysis import NetworkAnalysis
from mdc_analysis import MDCAnalyzer

def run_traditional_pipeline():
    # Follow Step 1-4 from Traditional WGCNA workflow
    # (See detailed code above)
    pass

def run_inter_tissue_pipeline():
    # Follow Step 1-4 from Inter-tissue workflow  
    # (See detailed code above)
    pass

if __name__ == "__main__":
    # Choose analysis type
    analysis_type = "traditional"  # or "inter_tissue"
    
    if analysis_type == "traditional":
        run_traditional_pipeline()
    else:
        run_inter_tissue_pipeline()
```

### Option 3: Testing and Validation

```bash
# Quick functionality test
python test_pipeline.py

# Specific component tests
python test_age_filtering.py
python test_inter_tissue.py
```

---

## 🔍 Key Differences from Original R Scripts

### What's New in Python Implementation:

1. **Modern Data Handling**: 
   - Polars for efficient large-data processing (noted in original but implementation uses pandas)
   - Streamlined GTEx GCT file parsing
   - Better memory management for large datasets

2. **Enhanced Quality Control**:
   - Integrated RIN score filtering
   - PCA-based outlier detection  
   - Comprehensive preprocessing pipeline

3. **Interactive Web Interface**:
   - Streamlit app for user-friendly analysis
   - Real-time parameter adjustment
   - Interactive visualizations

4. **Improved Statistical Methods**:
   - Robust permutation testing
   - Multiple FDR correction methods
   - Better handling of small sample sizes

### Missing Features (To Be Added):

1. **Advanced Module Analysis**:
   - Module eigengene calculation (present in R version)
   - Module-trait relationship analysis
   - Kruskal-Wallis testing for module significance

2. **Extended Cross-tissue Methods**:
   - Tissue-specific power selection (R version uses different powers per tissue)
   - Cross-tissue module preservation statistics
   - Tissue coordination index calculations

3. **Enhanced Visualizations**:
   - Cytoscape-style network plots (R version exports to Cytoscape)
   - Module preservation plots across conditions
   - Tissue-specific connectivity plots

4. **Pathway Analysis Integration**:
   - REACTOME pathway enrichment (referenced in R scripts)
   - GO term analysis for significant modules
   - Disease association analysis

---

## 📝 Notes for Implementation

### When Using Traditional WGCNA:
- Ensure adequate sample sizes (minimum 15-20 per group)
- Test multiple power values for optimal scale-free topology
- Consider signed vs unsigned networks based on biological interpretation

### When Using Inter-tissue Analysis:
- Verify sample overlap across tissues for meaningful cross-tissue correlations
- Use lower cross-tissue power (typically 3) compared to within-tissue power (typically 6)
- Filter to genes expressed across all tissues for valid comparisons

### Statistical Considerations:
- Use permutation testing for robust significance assessment
- Apply FDR correction for multiple testing
- Consider effect size in addition to statistical significance
- Validate results with independent datasets when possible

---

This guide provides the complete workflow for both traditional and inter-tissue coexpression network analysis using the modernized Python pipeline. The methodology maintains compatibility with the original R-based approach while offering improved usability and additional features.
