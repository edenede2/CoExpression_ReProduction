"""
Streamlit Web Application for Coexpression Network Analysis
Main application file for running the analysis pipeline
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from pathlib import Path
from typing import List, Dict, Tuple
import warnings
warnings.filterwarnings('ignore')

# Import custom modules
from data_reader import GTExDataReader
from preprocessing import DataPreprocessor
from network_analysis import NetworkAnalysis
from mdc_analysis import MDCAnalyzer
from visualization import NetworkVisualizer

# Configure Streamlit page
st.set_page_config(
    page_title="CoExpression Network Analysis",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Initialize session state
if 'data_loaded' not in st.session_state:
    st.session_state.data_loaded = False
if 'sample_limit' not in st.session_state:
    st.session_state.sample_limit = 10000
if 'preprocessed' not in st.session_state:
    st.session_state.preprocessed = False
if 'network_built' not in st.session_state:
    st.session_state.network_built = False

def load_data():
    """Load and cache the GTEx data."""
    if not st.session_state.data_loaded:
        with st.spinner("Loading GTEx data..."):
            reader = GTExDataReader()
            
            # Allow user to choose sample limit
            sample_limit = st.sidebar.slider(
                "Sample limit (for demo)", 
                min_value=100, 
                max_value=5000, 
                value=1000, 
                step=100,
                help="Limit the number of samples loaded for demonstration. Use larger values for more comprehensive analysis."
            )

            st.session_state.sample_limit = sample_limit

            # Load expression data with sample limit
            expression_df, gene_metadata = reader.read_gct_file(
                "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct",
                sample_limit=st.session_state.sample_limit
            )

            print(f"Loaded {len(expression_df)} genes and {len(expression_df.columns)} samples")
            
            # Load metadata
            sample_attrs = reader.read_sample_attributes()
            subject_phenos = reader.read_subject_phenotypes()
            protein_coding = reader.read_protein_coding_genes()
            
            # Store in session state
            st.session_state.expression_df = expression_df
            st.session_state.gene_metadata = gene_metadata
            st.session_state.sample_attrs = sample_attrs
            st.session_state.subject_phenos = subject_phenos
            st.session_state.protein_coding = protein_coding
            st.session_state.data_loaded = True
            
            st.success("Data loaded successfully!")

def main():
    """Main application function."""
    
    # Title and description
    st.title("🧬 CoExpression Network Analysis")
    st.markdown("""
    This application performs gene coexpression network analysis and Modular Differential Connectivity (MDC) 
    analysis using GTEx data. You can compare different groups (age, sex, disease status, etc.) and identify 
    modules with significantly different connectivity patterns.
    """)
    
    # Sidebar for navigation
    st.sidebar.title("Navigation")
    page = st.sidebar.selectbox(
        "Choose a page:",
        ["Data Loading", "Data Preprocessing", "Network Analysis", "MDC Analysis", "Results Visualization"]
    )
    
    if page == "Data Loading":
        data_loading_page()
    elif page == "Data Preprocessing":
        preprocessing_page()
    elif page == "Network Analysis":
        network_analysis_page()
    elif page == "MDC Analysis":
        mdc_analysis_page()
    elif page == "Results Visualization":
        visualization_page()

def data_loading_page():
    """Data loading and exploration page."""
    st.header("📊 Data Loading and Exploration")
    # Allow user to choose sample limit
    sample_limit = st.sidebar.slider(
        "Sample limit (for demo)", 
        min_value=1000, 
        max_value=50000, 
        value=10000, 
        step=1000,
        help="Limit the number of samples loaded for demonstration. Use larger values for more comprehensive analysis."
    )

    st.session_state.sample_limit = sample_limit

    # Load data button
    if st.button("Load GTEx Data", type="primary"):
        load_data()
    
    if st.session_state.data_loaded:
        st.success("✅ Data loaded successfully!")
        
        # Display data information
        col1, col2, col3 = st.columns(3)
        
        with col1:
            st.metric("Genes", len(st.session_state.expression_df))
        with col2:
            st.metric("Samples", len(st.session_state.expression_df.columns))
        with col3:
            st.metric("Protein Coding Genes", len(st.session_state.protein_coding))
        
        # Sample metadata exploration
        st.subheader("Sample Metadata")
        
        # Get available tissues with sample counts in loaded data
        reader = GTExDataReader()
        available_tissues = reader.get_available_tissues(st.session_state.sample_attrs)
        
        # Calculate tissue sample counts for loaded expression data
        loaded_samples = set(st.session_state.expression_df.columns)
        tissue_counts = []
        for tissue in available_tissues:
            all_tissue_samples = reader.get_tissue_samples(st.session_state.sample_attrs, tissue)
            loaded_tissue_samples = [s for s in all_tissue_samples if s in loaded_samples]
            tissue_counts.append(len(loaded_tissue_samples))
        
        st.write(f"**Available tissues in loaded dataset:**")
        tissue_df = pd.DataFrame({
            'Tissue': available_tissues,
            'Loaded Samples': tissue_counts,
            'Total Samples': [
                len(reader.get_tissue_samples(st.session_state.sample_attrs, tissue))
                for tissue in available_tissues
            ]
        })
        
        # Sort by loaded samples count
        tissue_df = tissue_df.sort_values('Loaded Samples', ascending=False)
        
        # Highlight tissues with good sample counts
        def highlight_good_tissues(row):
            if row['Loaded Samples'] >= 20:
                return ['background-color: green'] * len(row)
            elif row['Loaded Samples'] >= 10:
                return ['background-color: yellow'] * len(row)
            else:
                return ['background-color: coral'] * len(row)
        
        styled_df = tissue_df.style.apply(highlight_good_tissues, axis=1)
        st.dataframe(styled_df, use_container_width=True)
        
        # Age and sex distribution
        st.subheader("Subject Demographics")
        col1, col2 = st.columns(2)
        
        with col1:
            # Age distribution
            age_counts = st.session_state.subject_phenos['AGE'].value_counts().sort_index()
            fig_age = px.bar(
                x=age_counts.index,
                y=age_counts.values,
                title="Age Distribution",
                labels={'x': 'Age Group', 'y': 'Count'}
            )
            st.plotly_chart(fig_age, use_container_width=True)
        
        with col2:
            # Sex distribution
            sex_counts = st.session_state.subject_phenos['SEX'].value_counts()
            fig_sex = px.pie(
                values=sex_counts.values,
                names=['Male' if x == 1 else 'Female' for x in sex_counts.index],
                title="Sex Distribution"
            )
            st.plotly_chart(fig_sex, use_container_width=True)
    
    else:
        st.info("Click 'Load GTEx Data' to begin analysis.")

def preprocessing_page():
    """Data preprocessing page."""
    st.header("🔧 Data Preprocessing")
    
    if not st.session_state.data_loaded:
        st.warning("Please load data first on the Data Loading page.")
        return
    
    st.subheader("Preprocessing Options")
    
    # Preprocessing parameters
    col1, col2 = st.columns(2)
    
    with col1:
        apply_log = st.checkbox("Apply log₂(TPM + 1) transformation", value=True)
        filter_low_expr = st.checkbox("Filter low expression genes", value=True)
        filter_low_var = st.checkbox("Filter low variance genes", value=True)
        
    with col2:
        detect_outliers = st.checkbox("Detect and remove outlier samples", value=True)
        quantile_norm = st.checkbox("Apply quantile normalization", value=True)
        regress_confounders = st.checkbox("Regress out confounding factors", value=True)
        min_rin = st.slider("Minimum RIN score", 0.0, 10.0, 6.0, 0.1)
    
    # Advanced parameters
    with st.expander("Advanced Parameters"):
        min_expression = st.number_input("Minimum expression threshold", 0.0, 10.0, 0.1, 0.1)
        min_samples_fraction = st.slider("Minimum sample fraction for expression filter", 0.0, 1.0, 0.2, 0.05, 
                                         help="Genes must be expressed above threshold in at least this fraction of samples")
        min_variance = st.number_input("Minimum variance threshold", 0.0, 10.0, 0.1, 0.1)
        outlier_contamination = st.slider("Expected outlier fraction", 0.01, 0.2, 0.01, 0.01)
    
    # Group selection for analysis
    st.subheader("Group Selection for Analysis")
    
    # Tissue selection
    reader = GTExDataReader()
    available_tissues = reader.get_available_tissues(st.session_state.sample_attrs)
    
    # Filter tissues with sufficient samples in loaded data
    loaded_samples = set(st.session_state.expression_df.columns)
    tissues_with_samples = []
    tissue_sample_counts = {}
    
    for tissue in available_tissues:
        all_tissue_samples = reader.get_tissue_samples(st.session_state.sample_attrs, tissue)
        loaded_tissue_samples = [s for s in all_tissue_samples if s in loaded_samples]
        if len(loaded_tissue_samples) >= 10:  # Minimum threshold
            tissues_with_samples.append(tissue)
            tissue_sample_counts[tissue] = len(loaded_tissue_samples)
    
    # Sort tissues by sample count
    tissues_with_samples.sort(key=lambda x: tissue_sample_counts[x], reverse=True)
    
    if not tissues_with_samples:
        st.error("No tissues have sufficient samples in the loaded dataset. Try increasing the sample limit.")
        return
    
    selected_tissue = st.selectbox(
        "Select tissue for analysis:",
        tissues_with_samples,
        format_func=lambda x: f"{x} ({tissue_sample_counts[x]} samples)"
    )
    
    # Group comparison type
    comparison_type = st.selectbox(
        "Select comparison type:",
        ["Age (Young vs Old)", "Sex (Male vs Female)", "Custom"]
    )
    
    if st.button("Run Preprocessing", type="primary"):
        with st.spinner("Running preprocessing pipeline..."):
            preprocessor = DataPreprocessor()
            
            # Get available samples from loaded expression data
            available_samples = set(st.session_state.expression_df.columns)
            
            # Filter samples based on criteria
            if comparison_type == "Age (Young vs Old)":
                young_samples_all = reader.filter_samples_by_metadata(
                    st.session_state.sample_attrs,
                    st.session_state.subject_phenos,
                    min_rin=min_rin,
                    age_group='young',
                    tissue=selected_tissue
                )
                old_samples_all = reader.filter_samples_by_metadata(
                    st.session_state.sample_attrs,
                    st.session_state.subject_phenos,
                    min_rin=min_rin,
                    age_group='old',
                    tissue=selected_tissue
                )
                
                # Keep only samples present in expression data
                young_samples = [s for s in young_samples_all if s in available_samples]
                old_samples = [s for s in old_samples_all if s in available_samples]
                all_samples = young_samples + old_samples
                group_1_name, group_2_name = "Young", "Old"
                
            elif comparison_type == "Sex (Male vs Female)":
                male_samples_all = reader.filter_samples_by_metadata(
                    st.session_state.sample_attrs,
                    st.session_state.subject_phenos,
                    min_rin=min_rin,
                    sex=1,
                    tissue=selected_tissue
                )
                female_samples_all = reader.filter_samples_by_metadata(
                    st.session_state.sample_attrs,
                    st.session_state.subject_phenos,
                    min_rin=min_rin,
                    sex=2,
                    tissue=selected_tissue
                )
                
                # Keep only samples present in expression data
                male_samples = [s for s in male_samples_all if s in available_samples]
                female_samples = [s for s in female_samples_all if s in available_samples]
                all_samples = male_samples + female_samples
                young_samples, old_samples = male_samples, female_samples
                group_1_name, group_2_name = "Male", "Female"
            
            # Check if we have enough samples
            if len(all_samples) < 10:
                st.error(f"Not enough samples found! Only {len(all_samples)} samples available.")
                st.info("Try selecting a different tissue or reducing quality thresholds.")
                return
            
            if len(young_samples) < 3 or len(old_samples) < 3:
                st.error(f"Insufficient samples in groups: {group_1_name}: {len(young_samples)}, {group_2_name}: {len(old_samples)}")
                st.info("Need at least 3 samples per group for analysis.")
                return
            
            # Filter expression data to selected samples
            filtered_expression = st.session_state.expression_df[all_samples]
            
            # Prepare confounders for regression
            confounders_df = None
            
            if regress_confounders:
                with st.spinner("Preparing confounders..."):
                    confounders_df = prepare_confounders(
                        st.session_state.sample_attrs, 
                        st.session_state.subject_phenos,
                        all_samples
                    )
                    st.info(f"Prepared {confounders_df.shape[1]} confounding factors for regression")
            
            # Run preprocessing
            processed_df, processing_info = preprocessor.preprocess_expression_data(
                filtered_expression,
                st.session_state.protein_coding,
                confounders_df=confounders_df,
                apply_log=apply_log,
                filter_low_expr=filter_low_expr,
                filter_low_var=filter_low_var,
                detect_outliers=detect_outliers,
                quantile_norm=quantile_norm,
                regress_confounders=regress_confounders,
                min_expression=min_expression,
                min_samples_fraction=min_samples_fraction,
                min_variance=min_variance,
                outlier_contamination=outlier_contamination
            )
            
            # Store results
            st.session_state.processed_df = processed_df
            st.session_state.processing_info = processing_info
            st.session_state.young_samples = young_samples
            st.session_state.old_samples = old_samples
            st.session_state.group_1_name = group_1_name
            st.session_state.group_2_name = group_2_name
            st.session_state.selected_tissue = selected_tissue
            st.session_state.preprocessed = True
            
            st.success("Preprocessing completed!")
            
            # Show sample counts
            st.info(f"Analysis will use {len(young_samples)} {group_1_name} samples and {len(old_samples)} {group_2_name} samples from {selected_tissue} tissue.")
    
    # Show preprocessing results
    if st.session_state.preprocessed:
        st.subheader("Preprocessing Results")
        
        info = st.session_state.processing_info
        
        col1, col2, col3 = st.columns(3)
        with col1:
            st.metric("Original genes", info['original_shape'][0])
        with col2:
            st.metric("Final genes", info['final_shape'][0])
        with col3:
            st.metric("Final samples", info['final_shape'][1])
        
        st.write("**Processing steps applied:**", info['steps_applied'])
        
        # Show PCA plot if outlier detection was performed
        if 'pca_results' in info:
            st.subheader("Sample Quality Control")
            visualizer = NetworkVisualizer()
            pca_fig = visualizer.plot_sample_pca(info['pca_results'])
            st.plotly_chart(pca_fig, use_container_width=True)

def network_analysis_page():
    """Network analysis page with traditional and inter-tissue options."""
    st.header("🕸️ Network Analysis")
    
    if not st.session_state.preprocessed:
        st.warning("Please complete data preprocessing first.")
        return
    
    st.subheader("Analysis Type Selection")
    
    # Analysis type selection
    analysis_type = st.radio(
        "Choose analysis type:",
        ["Traditional WGCNA", "Inter-tissue X-WGCNA"],
        help="Traditional: Single condition analysis. Inter-tissue: Cross-tissue analysis with separate TS/CT powers."
    )
    
    if analysis_type == "Traditional WGCNA":
        st.info("🔬 **Traditional WGCNA**: Analyzes coexpression within a single condition/tissue")
        _traditional_network_analysis()
    else:
        st.info("🌐 **Inter-tissue X-WGCNA**: Analyzes coexpression across multiple tissues with tissue-specific (TS) and cross-tissue (CT) connections")
        _inter_tissue_network_analysis()

def _traditional_network_analysis():
    """Traditional WGCNA network analysis."""
    st.subheader("Traditional Network Construction Parameters")
    
    col1, col2 = st.columns(2)
    
    with col1:
        correlation_method = st.selectbox(
            "Correlation method:",
            ["pearson", "spearman"]
        )
        min_module_size = st.number_input("Minimum module size", 10, 200, 30, 10)
        
    with col2:
        auto_power = st.checkbox("Auto-determine soft threshold power", value=True)
        if not auto_power:
            manual_power = st.slider("Soft threshold power", 1, 20, 6, 1)
    
    if st.button("Build Traditional Network", type="primary"):
        with st.spinner("Building coexpression network..."):
            network_analyzer = NetworkAnalysis()
            network_analyzer.set_analysis_type('traditional')
            
            power = None if auto_power else manual_power
            
            # Build network using all samples
            network_results = network_analyzer.construct_network(
                st.session_state.processed_df,
                power=power,
                min_module_size=min_module_size,
                correlation_method=correlation_method
            )
            
            # Calculate module eigengenes
            eigengenes = network_analyzer.calculate_module_eigengenes(
                st.session_state.processed_df,
                network_results['modules']
            )
            
            # Store results
            st.session_state.network_results = network_results
            st.session_state.eigengenes = eigengenes
            st.session_state.network_built = True
            st.session_state.analysis_type = 'traditional'
            
            st.success("Traditional network construction completed!")
    
    # Show results if available
    if (hasattr(st.session_state, 'network_built') and st.session_state.network_built and 
        hasattr(st.session_state, 'analysis_type') and st.session_state.analysis_type == 'traditional'):
        _display_network_results()

def _inter_tissue_network_analysis():
    """Inter-tissue X-WGCNA network analysis."""
    st.subheader("Inter-tissue Network Construction")
    
    # Step 1: Tissue grouping
    st.write("**Step 1: Tissue Grouping**")
    
    if 'filtered_metadata' not in st.session_state:
        st.error("No filtered metadata available. Please complete preprocessing first.")
        return
    
    # Get available tissues
    available_tissues = st.session_state.filtered_metadata['SMTSD'].unique()
    
    selected_tissues = st.multiselect(
        "Select tissues for inter-tissue analysis:",
        available_tissues,
        default=available_tissues[:3] if len(available_tissues) >= 3 else available_tissues,
        help="Select at least 2 tissues for inter-tissue analysis"
    )
    
    if len(selected_tissues) < 2:
        st.warning("Please select at least 2 tissues for inter-tissue analysis.")
        return
    
    # Step 2: Power determination
    st.write("**Step 2: Power Determination**")
    
    col1, col2 = st.columns(2)
    with col1:
        auto_determine_powers = st.checkbox("Auto-determine TS and CT powers", value=True)
        target_rsq = st.slider("Target R² for scale-free topology", 0.5, 0.95, 0.8, 0.05)
    
    with col2:
        if not auto_determine_powers:
            ts_power = st.slider("Tissue-Specific (TS) power", 1, 20, 6, 1)
            ct_power = st.slider("Cross-Tissue (CT) power", 1, 20, 3, 1)
    
    # Step 3: Network parameters
    st.write("**Step 3: Network Parameters**")
    
    col1, col2 = st.columns(2)
    with col1:
        correlation_method = st.selectbox("Correlation method:", ["pearson", "spearman"])
        min_module_size = st.number_input("Minimum module size", 10, 200, 30, 10)
    
    with col2:
        cluster_type_threshold = st.slider(
            "Cluster type threshold", 
            0.5, 1.0, 0.95, 0.05,
            help="Threshold for classifying modules as tissue-specific (TS) vs cross-tissue (CT)"
        )
    
    # Power determination button
    if auto_determine_powers and st.button("🔍 Determine Optimal Powers"):
        with st.spinner("Determining optimal powers for tissues..."):
            try:
                # Prepare tissue data
                tissue_data, tissue_sample_mapping = _prepare_tissue_data_for_inter_tissue(selected_tissues)
                
                # Initialize analyzer
                network_analyzer = NetworkAnalysis()
                network_analyzer.set_analysis_type('inter_tissue')
                network_analyzer.load_tissue_data_for_inter_tissue_analysis(tissue_data, tissue_sample_mapping)
                
                # Determine powers
                power_results = network_analyzer.determine_inter_tissue_powers(target_rsq=target_rsq)
                
                st.session_state.power_results = power_results
                st.session_state.tissue_data_prepared = True
                st.session_state.network_analyzer = network_analyzer
                
                st.success("Power determination completed!")
                
                # Display results
                st.write("**Power Analysis Results:**")
                col1, col2 = st.columns(2)
                with col1:
                    st.metric("Recommended TS Power", f"{power_results['recommended_ts_power']:.1f}")
                with col2:
                    st.metric("Recommended CT Power", f"{power_results['recommended_ct_power']:.1f}")
                
                # Show detailed results
                with st.expander("📊 Detailed Power Analysis"):
                    st.write("**Tissue-Specific Powers:**")
                    ts_results = []
                    for tissue, analysis in power_results['tissue_specific'].items():
                        ts_results.append({
                            'Tissue': tissue,
                            'Optimal_Power': analysis['optimal_power'],
                            'R²': f"{analysis['optimal_rsq']:.3f}"
                        })
                    st.dataframe(pd.DataFrame(ts_results))
                    
                    st.write("**Cross-Tissue Powers:**")
                    ct_results = []
                    for pair, analysis in power_results['cross_tissue'].items():
                        ct_results.append({
                            'Tissue_Pair': pair,
                            'Optimal_Power': analysis['optimal_power'],
                            'R²': f"{analysis['optimal_rsq']:.3f}"
                        })
                    st.dataframe(pd.DataFrame(ct_results))
                
            except Exception as e:
                st.error(f"Error in power determination: {str(e)}")
                return
    
    # Network construction button
    build_network = False
    if auto_determine_powers:
        if hasattr(st.session_state, 'power_results'):
            ts_power_to_use = st.session_state.power_results['recommended_ts_power']
            ct_power_to_use = st.session_state.power_results['recommended_ct_power']
            build_network = st.button("🕸️ Build Inter-tissue Network", type="primary")
    else:
        ts_power_to_use = ts_power
        ct_power_to_use = ct_power
        build_network = st.button("🕸️ Build Inter-tissue Network", type="primary")
    
    if build_network:
        with st.spinner("Building inter-tissue network..."):
            try:
                # Prepare data if not already done
                if not hasattr(st.session_state, 'tissue_data_prepared'):
                    tissue_data, tissue_sample_mapping = _prepare_tissue_data_for_inter_tissue(selected_tissues)
                    network_analyzer = NetworkAnalysis()
                    network_analyzer.set_analysis_type('inter_tissue')
                    network_analyzer.load_tissue_data_for_inter_tissue_analysis(tissue_data, tissue_sample_mapping)
                else:
                    network_analyzer = st.session_state.network_analyzer
                
                # Build network
                network_results = network_analyzer.construct_inter_tissue_network(
                    ts_power=ts_power_to_use,
                    ct_power=ct_power_to_use,
                    correlation_method=correlation_method,
                    min_module_size=min_module_size
                )
                
                # Store results
                st.session_state.network_results = network_results
                st.session_state.network_built = True
                st.session_state.analysis_type = 'inter_tissue'
                st.session_state.selected_tissues = selected_tissues
                
                st.success("Inter-tissue network construction completed!")
                
            except Exception as e:
                st.error(f"Error in network construction: {str(e)}")
                return
    
    # Show results if available
    if (hasattr(st.session_state, 'network_built') and st.session_state.network_built and 
        hasattr(st.session_state, 'analysis_type') and st.session_state.analysis_type == 'inter_tissue'):
        _display_inter_tissue_results()

def _prepare_tissue_data_for_inter_tissue(selected_tissues: List[str]) -> Tuple[Dict, Dict]:
    """Prepare tissue data for inter-tissue analysis."""
    tissue_data = {}
    tissue_sample_mapping = {}
    
    for tissue in selected_tissues:
        # Get samples for this tissue
        tissue_samples = st.session_state.filtered_metadata[
            st.session_state.filtered_metadata['SMTSD'] == tissue
        ]['SAMPID'].tolist()
        
        # Get expression data for these samples
        tissue_expr = st.session_state.processed_df[tissue_samples]
        
        tissue_data[tissue] = tissue_expr.T  # Samples x genes
        tissue_sample_mapping[tissue] = tissue_samples
    
    return tissue_data, tissue_sample_mapping

def _display_network_results():
    """Display traditional network analysis results."""
    st.subheader("Network Analysis Results")
    
    results = st.session_state.network_results
    
    col1, col2, col3 = st.columns(3)
    with col1:
        st.metric("Soft threshold power", results['power'])
    with col2:
        n_modules = len(set(results['modules'].values())) - (1 if 0 in results['modules'].values() else 0)
        st.metric("Number of modules", n_modules)
    with col3:
        grey_genes = sum(1 for m in results['modules'].values() if m == 0)
        st.metric("Unassigned genes", grey_genes)
    
    # Visualizations
    try:
        visualizer = NetworkVisualizer()
        
        # Power analysis plot
        if results['power_analysis'] is not None:
            st.subheader("Soft Threshold Power Selection")
            power_fig = visualizer.plot_power_analysis(results['power_analysis'])
            st.plotly_chart(power_fig, use_container_width=True)
        
        # Module sizes
        st.subheader("Module Size Distribution")
        module_fig = visualizer.plot_module_sizes(results['modules'])
        st.plotly_chart(module_fig, use_container_width=True)
    except Exception as e:
        st.warning(f"Visualization error: {e}")
        st.write("Network results available but visualization failed.")

def _display_inter_tissue_results():
    """Display inter-tissue network analysis results."""
    st.subheader("Inter-tissue Network Analysis Results")
    
    results = st.session_state.network_results
    
    # Basic metrics
    col1, col2, col3, col4 = st.columns(4)
    with col1:
        st.metric("TS Power", results['ts_power'])
    with col2:
        st.metric("CT Power", results['ct_power'])
    with col3:
        n_modules = len(results['module_details'])
        st.metric("Number of modules", n_modules)
    with col4:
        grey_genes = sum(1 for m in results['modules'].values() if m == 0)
        st.metric("Unassigned genes", grey_genes)
    
    # Module composition analysis
    st.subheader("Module Composition Analysis")
    
    module_details = results['module_details']
    
    # TS vs CT module classification
    ts_modules = len(module_details[module_details['Module_Type'] == 'TS'])
    ct_modules = len(module_details[module_details['Module_Type'] == 'CT'])
    
    col1, col2 = st.columns(2)
    with col1:
        st.metric("Tissue-Specific (TS) Modules", ts_modules)
    with col2:
        st.metric("Cross-Tissue (CT) Modules", ct_modules)
    
    # Display module details table
    st.write("**Module Details:**")
    st.dataframe(module_details)
    
    # Network validation
    if 'validation' in results:
        st.subheader("Network Validation")
        validation = results['validation']
        
        # Overall network properties
        st.write("**Overall Network Properties:**")
        col1, col2, col3 = st.columns(3)
        with col1:
            st.metric("Scale-free R²", f"{validation['overall_network']['rsq']:.3f}")
        with col2:
            st.metric("Network Slope", f"{validation['overall_network']['slope']:.3f}")
        with col3:
            st.metric("Mean Connectivity", f"{validation['overall_network']['mean_connectivity']:.2f}")
        
        # Tissue-specific validation
        with st.expander("📊 Detailed Validation Results"):
            st.write("**Tissue-Specific Network Properties:**")
            ts_validation = []
            for tissue, props in validation['tissue_specific'].items():
                ts_validation.append({
                    'Tissue': tissue,
                    'R²': f"{props['rsq']:.3f}",
                    'Slope': f"{props['slope']:.3f}",
                    'Mean_Connectivity': f"{props['mean_connectivity']:.2f}"
                })
            st.dataframe(pd.DataFrame(ts_validation))
            
            st.write("**Cross-Tissue Network Properties:**")
            ct_validation = []
            for pair, props in validation['cross_tissue'].items():
                ct_validation.append({
                    'Tissue_Pair': pair,
                    'Mean_Connectivity': f"{props['mean_connectivity']:.2f}",
                    'Median_Connectivity': f"{props['median_connectivity']:.2f}"
                })
            st.dataframe(pd.DataFrame(ct_validation))
        
        # Power analysis plot
        if results['power_analysis'] is not None:
            try:
                visualizer = NetworkVisualizer()
                st.subheader("Soft Threshold Power Selection")
                power_fig = visualizer.plot_power_analysis(results['power_analysis'])
                st.plotly_chart(power_fig, use_container_width=True)
            except Exception as e:
                st.warning(f"Power analysis visualization failed: {e}")
        
        # Module sizes
        try:
            visualizer = NetworkVisualizer()
            st.subheader("Module Size Distribution")
            module_fig = visualizer.plot_module_sizes(results['modules'])
            st.plotly_chart(module_fig, use_container_width=True)
        except Exception as e:
            st.warning(f"Module size visualization failed: {e}")
        
        # Module eigengenes
        if len(st.session_state.eigengenes) > 0:
            try:
                visualizer = NetworkVisualizer()
                st.subheader("Module Eigengenes")
                
                # Create sample groups for visualization
                sample_groups = {}
                for sample in st.session_state.eigengenes.index:
                    if sample in st.session_state.young_samples:
                        sample_groups[sample] = st.session_state.group_1_name
                    elif sample in st.session_state.old_samples:
                        sample_groups[sample] = st.session_state.group_2_name
                    else:
                        sample_groups[sample] = "Other"
                
                eigengene_fig = visualizer.plot_module_eigengenes(
                    st.session_state.eigengenes,
                    sample_groups
                )
                st.plotly_chart(eigengene_fig, use_container_width=True)
            except Exception as e:
                st.warning(f"Eigengene visualization failed: {e}")

def mdc_analysis_page():
    """MDC analysis page."""
    st.header("📊 Modular Differential Connectivity Analysis")
    
    if not st.session_state.network_built:
        st.warning("Please complete network analysis first.")
        return
    
    st.subheader("MDC Analysis Parameters")
    
    col1, col2 = st.columns(2)
    
    with col1:
        run_permutation = st.checkbox("Run permutation test", value=True)
        if run_permutation:
            n_permutations = st.number_input("Number of permutations", 10, 1000, 100, 10)
        
    with col2:
        p_threshold = st.number_input("P-value threshold", 0.001, 0.1, 0.05, 0.001, format="%.3f")
        log_ratio_threshold = st.number_input("Minimum log ratio threshold", 0.1, 2.0, 0.5, 0.1)
    
    if st.button("Run MDC Analysis", type="primary"):
        with st.spinner("Running MDC analysis..."):
            mdc_analyzer = MDCAnalyzer()
            
            # Split expression data by groups
            group_1_samples = [s for s in st.session_state.young_samples if s in st.session_state.processed_df.columns]
            group_2_samples = [s for s in st.session_state.old_samples if s in st.session_state.processed_df.columns]
            
            expr_group_1 = st.session_state.processed_df[group_1_samples]
            expr_group_2 = st.session_state.processed_df[group_2_samples]
            
            # Run MDC analysis
            mdc_results = mdc_analyzer.analyze_differential_connectivity(
                expr_group_1,
                expr_group_2,
                st.session_state.network_results['modules'],
                st.session_state.network_results['power'],
                st.session_state.group_1_name,
                st.session_state.group_2_name,
                run_permutation_test=run_permutation,
                n_permutations=n_permutations if run_permutation else 0
            )
            
            st.session_state.mdc_results = mdc_results
            st.session_state.mdc_analyzer = mdc_analyzer
            
            st.success("MDC analysis completed!")
    
    # Show MDC results
    if 'mdc_results' in st.session_state:
        st.subheader("MDC Analysis Results")
        
        mdc_results = st.session_state.mdc_results['mdc_results']
        
        # Summary statistics
        col1, col2, col3 = st.columns(3)
        with col1:
            st.metric("Modules analyzed", len(mdc_results))
        with col2:
            increased = len(mdc_results[mdc_results['log_mdc_ratio'] > 0])
            st.metric("Increased connectivity", increased)
        with col3:
            decreased = len(mdc_results[mdc_results['log_mdc_ratio'] < 0])
            st.metric("Decreased connectivity", decreased)
        
        # Visualizations
        visualizer = NetworkVisualizer()
        
        st.subheader("MDC Results Overview")
        mdc_fig = visualizer.plot_mdc_results(
            mdc_results,
            st.session_state.group_1_name,
            st.session_state.group_2_name
        )
        st.plotly_chart(mdc_fig, use_container_width=True)
        
        # Significant modules
        if 'permutation_results' in st.session_state.mdc_results:
            significant_modules = st.session_state.mdc_analyzer.get_significant_modules(
                p_threshold, log_ratio_threshold
            )
            
            if len(significant_modules) > 0:
                st.subheader("Significant Modules")
                
                sig_fig = visualizer.plot_significant_modules(significant_modules, top_n=10)
                st.plotly_chart(sig_fig, use_container_width=True)
                
                # Table of significant modules
                st.subheader("Significant Modules Table")
                display_cols = ['module_id', 'module_size', 'log_mdc_ratio', 'p_value_corrected']
                st.dataframe(
                    significant_modules[display_cols].round(4),
                    use_container_width=True
                )
            else:
                st.info("No significant modules found with current thresholds.")
        
        # Download results
        st.subheader("Download Results")
        
        # Prepare results for download
        results_data = mdc_results.copy()
        if 'permutation_results' in st.session_state.mdc_results:
            perm_results = st.session_state.mdc_results['permutation_results']
            results_data = results_data.merge(
                perm_results[['module_id', 'p_value', 'p_value_corrected', 'significant']],
                on='module_id',
                how='left'
            )
        
        csv = results_data.to_csv(index=False)
        st.download_button(
            label="Download MDC results as CSV",
            data=csv,
            file_name=f"mdc_results_{st.session_state.selected_tissue}_{st.session_state.group_1_name}_vs_{st.session_state.group_2_name}.csv",
            mime="text/csv"
        )

def visualization_page():
    """Results visualization page."""
    st.header("📈 Results Visualization")
    
    if 'mdc_results' not in st.session_state:
        st.warning("Please complete MDC analysis first.")
        return
    
    visualizer = NetworkVisualizer()
    
    # Network visualization
    st.subheader("Network Visualization")
    
    col1, col2 = st.columns(2)
    with col1:
        max_nodes = st.slider("Maximum nodes to display", 50, 500, 100, 50)
    with col2:
        layout_type = st.selectbox("Layout algorithm", ["spring", "circular"])
    
    if st.button("Generate Network Plot"):
        with st.spinner("Generating network visualization..."):
            network_fig = visualizer.create_network_plot(
                st.session_state.network_results['adjacency_matrix'],
                st.session_state.network_results['modules'],
                layout=layout_type,
                max_nodes=max_nodes
            )
            st.plotly_chart(network_fig, use_container_width=True)
    
    # Additional visualizations
    st.subheader("Correlation Analysis")
    
    if st.button("Show Correlation Heatmap"):
        with st.spinner("Generating correlation heatmap..."):
            corr_fig = visualizer.plot_correlation_heatmap(
                st.session_state.network_results['correlation_matrix'],
                max_genes=100
            )
            st.plotly_chart(corr_fig, use_container_width=True)
    
    # Summary report
    st.subheader("Analysis Summary")
    
    summary_data = {
        "Analysis Parameter": [
            "Tissue",
            "Comparison",
            "Total Genes",
            "Total Samples",
            "Modules Detected",
            "Soft Threshold Power",
            "Correlation Method"
        ],
        "Value": [
            st.session_state.selected_tissue,
            f"{st.session_state.group_1_name} vs {st.session_state.group_2_name}",
            st.session_state.processed_df.shape[0],
            st.session_state.processed_df.shape[1],
            len(set(st.session_state.network_results['modules'].values())),
            st.session_state.network_results['power'],
            "Pearson"  # Default
        ]
    }
    
    summary_df = pd.DataFrame(summary_data)
    st.table(summary_df)

def prepare_confounders(sample_attrs, subject_phenos, samples):
    """
    Prepare confounders DataFrame for regression, following the original script's logic.
    
    Args:
        sample_attrs: Sample attributes DataFrame
        subject_phenos: Subject phenotypes DataFrame
        samples: List of samples to include
        
    Returns:
        DataFrame of confounders ready for regression
    """
    # Filter sample attributes to selected samples
    filtered_attrs = sample_attrs.loc[list(set(samples) & set(sample_attrs.index))].copy()
    
    if len(filtered_attrs) == 0:
        st.error("No samples found in sample attributes")
        return pd.DataFrame(index=samples)
    
    # Extract subject IDs
    filtered_attrs['SUBJID'] = filtered_attrs.index.str.extract(r'(GTEX-[^-]+)', expand=False)
    
    # Merge with phenotypes to get AGE and DTHHRDY
    merged_data = filtered_attrs.merge(
        subject_phenos[['AGE', 'DTHHRDY']],
        left_on='SUBJID',
        right_index=True,
        how='left'
    )
    
    # Process batch information (SMGEBTCH) similar to original script
    if 'SMGEBTCH' in merged_data.columns:
        # Count batch occurrences
        batch_counts = merged_data['SMGEBTCH'].value_counts()
        valid_batches = batch_counts[batch_counts > 1].index.tolist()
        
        # Replace rare batches with singleton category
        merged_data['SMGEBTCH'] = merged_data['SMGEBTCH'].apply(
            lambda x: x if x in valid_batches else 'ASINGLETON_SMGEBTCH'
        )
    
    # Remove unnecessary columns
    cols_to_drop = ['SMTSD', 'SUBJID']
    merged_data = merged_data.drop([c for c in cols_to_drop if c in merged_data.columns], axis=1)
    
    # Ensure numeric columns are numeric and categorical are properly encoded
    numeric_cols = []
    categorical_cols = []
    
    for col in merged_data.columns:
        try:
            # Try to convert to numeric
            merged_data[col] = pd.to_numeric(merged_data[col])
            numeric_cols.append(col)
        except:
            # If it fails, it's categorical
            categorical_cols.append(col)
    
    print(f"Identified {len(numeric_cols)} numeric confounders and {len(categorical_cols)} categorical confounders")
    
    # Keep only columns that can be used for regression
    final_data = merged_data[numeric_cols + categorical_cols].copy()
    
    # Handle special case for AGE - convert from ranges to numeric values
    if 'AGE' in final_data.columns and isinstance(final_data['AGE'].iloc[0], str):
        try:
            # Extract lower bound of age range (e.g., "20-29" -> 20)
            final_data['AGE'] = final_data['AGE'].str.split('-').str[0].astype(int)
            print("Converted AGE to numeric")
        except:
            print("Could not convert AGE to numeric")
    
    # Final confounders DataFrame 
    # In the original script: x = sample.drop(['AGE', 'SMRIN', 'DTHHRDY'], axis=1)
    # We'll keep AGE, SMRIN, DTHHRDY in our confounders dataframe since they'll be 
    # handled when we create dummy variables in the regression function
    
    return final_data

if __name__ == "__main__":
    main()
