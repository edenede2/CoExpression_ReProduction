# CoExpression Network Analysis Pipeline

A comprehensive Python-based pipeline for gene coexpression network analysis and Modular Differential Connectivity (MDC) analysis using GTEx data.

## 🔬 Overview

This pipeline modernizes gene coexpression network analysis by implementing WGCNA (Weighted Gene Co-expression Network Analysis) methodology in Python with a user-friendly Streamlit interface. It enables researchers to:

- Compare gene expression patterns between different groups (age, sex, disease, custom)
- Build coexpression networks and identify gene modules
- Calculate Modular Differential Connectivity (MDC) between conditions
- Perform statistical testing with permutation tests and FDR correction

## ✨ Features

- **Flexible Data Handling**: Optimized for large GTEx datasets (GCT format)
- **Robust Quality Control**: RIN filtering, PCA-based outlier detection, variance filtering
- **Modern Network Analysis**: WGCNA-style coexpression networks with soft thresholding
- **Advanced Module Detection**: Hierarchical clustering and topological overlap matrices
- **Statistical Rigor**: Permutation testing with FDR correction for multiple comparisons
- **Interactive Web Interface**: Streamlit app with step-by-step workflow
- **Flexible Group Comparisons**: Age groups, sex, disease status, or custom groupings
- **Rich Visualizations**: Interactive plotly graphics and network diagrams

## 📋 Requirements

- Python 3.8+
- 8GB+ RAM (for full GTEx dataset)
- Modern web browser (for Streamlit interface)

## 📦 Installation

### 1. Navigate to the Project Directory
```bash
cd /path/to/CoExpression_reProduction/new_script
```

### 2. Install Dependencies
```bash
pip install -r requirements.txt
```

### 3. Verify Data Files
Ensure GTEx data files are in the `../data` directory:
- `GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct`
- `GTEx_Analysis_v8_Annotations_SampleAttributesDS.tsv`
- `GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.tsv`
- `hgnc_complete_set.tsv`

## 🚀 Quick Start

### Automated Testing (Recommended)

From the project root directory:

```bash
./run_tests.sh
```

This script automatically:
- Finds the correct directories
- Checks for required data files
- Runs all tests in the proper order
- Provides clear next steps

### Manual Testing

Navigate to the script directory first:

```bash
cd new_script
```

Then run tests:

```bash
python minimal_test.py        # Quick data loading verification
python test_age_filtering.py  # Age filtering functionality test  
python test_pipeline.py       # Complete pipeline test (may take time)
```

### Launch Web Application

```bash
cd new_script
streamlit run streamlit_app.py
```

**Important Notes:**
- The app loads a subset of samples for demonstration (adjustable in sidebar)
- Select tissues with good sample counts (highlighted in green/yellow)
- Ensure sufficient samples per group for meaningful analysis
- See `STREAMLIT_FIX.md` for detailed usage instructions

## 📖 Usage Guide

### Command Line Interface

```python
from data_reader import GTExDataReader
from preprocessing import DataPreprocessor
from network_analysis import NetworkAnalyzer
from mdc_analysis import MDCAnalyzer

# 1. Load and filter data
reader = GTExDataReader()
expression_df, gene_metadata = reader.read_gct_file("expression_file.gct")
sample_attrs = reader.read_sample_attributes()
subject_phenos = reader.read_subject_phenotypes()

# Filter samples by age and tissue
young_samples = reader.filter_samples_by_metadata(
    sample_attrs, subject_phenos, 
    age_group='young', tissue='Brain - Cortex'
)

old_samples = reader.filter_samples_by_metadata(
    sample_attrs, subject_phenos, 
    age_group='old', tissue='Brain - Cortex'
)

# 2. Preprocess expression data
preprocessor = DataPreprocessor()
processed_df, info = preprocessor.preprocess_expression_data(
    expression_df[young_samples + old_samples]
)

# 3. Build coexpression network
network_analyzer = NetworkAnalyzer()
network_results = network_analyzer.construct_network(processed_df)

# 4. Analyze differential connectivity
mdc_analyzer = MDCAnalyzer()
mdc_results = mdc_analyzer.analyze_differential_connectivity(
    processed_df[young_samples], 
    processed_df[old_samples],
    network_results['modules'], 
    network_results['power'],
    "Young", "Old"
)
```

### Web Interface Workflow

1. **📂 Data Loading**: Load GTEx expression data and metadata
2. **🔧 Preprocessing**: Configure quality control and filtering parameters
3. **🕸️ Network Analysis**: Build coexpression networks and detect modules
4. **📊 MDC Analysis**: Compare connectivity between groups with statistical testing
5. **📈 Visualization**: Interactive plots and network visualizations
6. **💾 Export Results**: Download analysis results and plots

## 🏗️ Module Architecture

### `data_reader.py` - Data Loading & Filtering
- Efficient GTEx GCT file parsing
- Sample metadata integration
- Advanced filtering by tissue, age, sex, and quality metrics
- **Age Groups Available**:
  - `'young'`: Ages 20-29 and 30-39
  - `'old'`: Ages 60-69 and 70-79  
  - `'middle'`: Ages 40-49 and 50-59
  - Custom lists: e.g., `['20-29', '70-79']`

### `preprocessing.py` - Quality Control Pipeline
- Log transformation and quantile normalization
- Gene filtering (expression level, variance, protein-coding)
- PCA-based outlier detection with Mahalanobis distance
- Sample quality filtering (RIN scores)

### `network_analysis.py` - Coexpression Networks
- Correlation matrix calculation (Pearson, Spearman, Kendall)
- Soft threshold power optimization for scale-free topology
- Adjacency matrix and Topological Overlap Matrix (TOM) construction
- Module detection via hierarchical clustering
- Module eigengene calculation

### `mdc_analysis.py` - Differential Connectivity
- Modular Differential Connectivity (MDC) calculation
- Permutation testing for statistical significance
- False Discovery Rate (FDR) correction
- Cross-condition network comparison metrics

### `visualization.py` - Interactive Plots
- Power analysis and scale-free topology plots
- PCA plots with outlier highlighting
- Module size distributions and eigengene heatmaps
- MDC results visualization
- Network topology plots

### `streamlit_app.py` - Web Interface
- Multi-page application with guided workflow
- Real-time parameter adjustment
- Progress tracking and result caching
- Interactive result exploration and download

## ⚙️ Configuration Parameters

### Data Filtering
- **Minimum RIN**: RNA Integrity Number threshold (default: 6.0)
- **Expression Filter**: Minimum expression level and sample count
- **Variance Filter**: Minimum gene variance threshold
- **Protein Coding**: Filter for protein-coding genes only

### Network Construction
- **Correlation Method**: Pearson (default), Spearman, or Kendall
- **Soft Threshold Power**: Auto-determined (1-20) or manual selection
- **Minimum Module Size**: Minimum genes per module (default: 30)
- **Merge Threshold**: Module merging similarity threshold (default: 0.25)

### Statistical Testing
- **Permutation Count**: Number of permutation tests (default: 100)
- **P-value Threshold**: Significance threshold (default: 0.05)
- **Effect Size Threshold**: Minimum log MDC ratio (default: 0.5)
- **FDR Method**: Multiple testing correction method (default: 'fdr_bh')

### Quality Control
- **Outlier Detection**: Expected contamination fraction (default: 0.05)
- **Sample Filtering**: Minimum sample count per condition
- **Gene Filtering**: Expression and variance thresholds

## 📊 Example Applications

### 1. Age-Related Analysis
```python
# Compare young vs old brain samples
young_brain = reader.filter_samples_by_metadata(
    sample_attrs, subject_phenos,
    age_group='young', tissue='Brain - Cortex'
)

old_brain = reader.filter_samples_by_metadata(
    sample_attrs, subject_phenos,
    age_group='old', tissue='Brain - Cortex'
)
```

### 2. Sex-Based Comparison
```python
# Compare male vs female samples
male_samples = reader.filter_samples_by_metadata(
    sample_attrs, subject_phenos,
    sex=1, tissue='Heart - Left Ventricle'
)

female_samples = reader.filter_samples_by_metadata(
    sample_attrs, subject_phenos,
    sex=2, tissue='Heart - Left Ventricle'
)
```

### 3. Custom Group Analysis
```python
# Compare extreme age groups
very_young = reader.filter_samples_by_metadata(
    sample_attrs, subject_phenos,
    age_group=['20-29'], tissue='Liver'
)

very_old = reader.filter_samples_by_metadata(
    sample_attrs, subject_phenos,
    age_group=['70-79'], tissue='Liver'
)
```

## 📈 Output Files

### Analysis Results
- **MDC Results**: `mdc_results.csv` - Module connectivity metrics and p-values
- **Module Assignments**: `gene_modules.csv` - Gene-to-module mapping
- **Network Matrices**: `adjacency_matrix.csv`, `tom_matrix.csv`
- **Sample Information**: `processed_samples.csv` - QC and filtering details

### Visualizations
- **Power Analysis Plot**: Soft threshold selection diagnostics
- **Module Plots**: Size distribution, dendrograms, eigengene heatmaps
- **MDC Plots**: Differential connectivity results and significance
- **Quality Control Plots**: PCA with outliers, expression distributions

## 🔧 Troubleshooting

### Common Issues

**Memory Issues with Large Datasets**
```python
# Use sample limits for testing
expression_df, _ = reader.read_gct_file("file.gct", sample_limit=1000)
```

**No Samples Found After Filtering**
```python
# Check available age groups
age_groups = reader.get_available_age_groups(subject_phenos)
print("Available age groups:", age_groups)

# Check tissue sample counts
tissues = reader.get_available_tissues(sample_attrs)
for tissue in tissues[:5]:
    count = len(reader.get_tissue_samples(sample_attrs, tissue))
    print(f"{tissue}: {count} samples")
```

**Installation Issues**
```bash
# Update pip and install dependencies
pip install --upgrade pip
pip install -r requirements.txt --upgrade
```

## 🧪 Testing

### Comprehensive Pipeline Test
```bash
python test_pipeline.py
```

### Age Filtering Verification
```bash
python test_age_filtering.py
```

### Quick Functionality Check
```bash
python quick_test.py
```

## 📚 Scientific Background

This pipeline implements methods described in:

- **WGCNA**: Zhang & Horvath (2005), Langfelder & Horvath (2008)
- **MDC**: Oldham et al. (2008), Functional organization of the transcriptome
- **GTEx**: The GTEx Consortium (2020), The GTEx Consortium atlas of genetic regulatory effects

### Key Methodological Features

1. **Scale-Free Network Construction**: Uses soft thresholding to preserve continuous nature of correlations
2. **Topological Overlap Matrix**: Emphasizes shared network neighborhoods for robust module detection
3. **Permutation Testing**: Provides empirical p-values for differential connectivity
4. **Multiple Testing Correction**: Controls false discovery rate across modules

## 📄 License

This project is licensed under the MIT License - see the LICENSE file for details.

## 🤝 Contributing

Contributions are welcome! Please feel free to submit issues, feature requests, or pull requests.

## 📞 Support

For questions or issues:
1. Check the test script outputs for debugging information
2. Review the Streamlit app logs for detailed error messages
3. Examine the example usage in the test files

---

**Ready to explore gene coexpression networks? Start with `python test_pipeline.py` and then launch the web app with `streamlit run streamlit_app.py`!** 🚀
