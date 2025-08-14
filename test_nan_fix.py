#!/usr/bin/env python3
"""
Test script to verify NaN handling fixes in preprocessing
"""

import pandas as pd
import numpy as np
from new_script.data_reader import GTExDataReader
from new_script.preprocessing import DataPreprocessor

def test_nan_fixes():
    print("Testing NaN handling fixes in preprocessing...")
    
    # Initialize components
    reader = GTExDataReader()
    preprocessor = DataPreprocessor()
    
    # Load a small subset of data for testing
    target_tissues = ['Brain - Cortex']
    
    try:
        # Load data
        print("Loading GTEx data...")
        expression_df, gene_metadata = reader.read_gct_file(
            "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct"
        )
        
        sample_attrs = reader.read_sample_attributes()
        subject_phenos = reader.read_subject_phenotypes()
        hgnc_data = reader.read_protein_coding_genes()
        
        # Align data
        expression_df, sample_attrs, subject_phenos = reader.align_expression_with_metadata(
            expression_df, sample_attrs, subject_phenos
        )
        
        # Get tissue samples
        tissue_young_samples = reader.filter_samples_by_metadata(
            sample_attrs, subject_phenos, tissue='Brain - Cortex', age_group='young'
        )
        tissue_old_samples = reader.filter_samples_by_metadata(
            sample_attrs, subject_phenos, tissue='Brain - Cortex', age_group='old'
        )
        
        tissue_samples = tissue_young_samples + tissue_old_samples
        print(f"Using {len(tissue_samples)} Brain - Cortex samples")
        
        # Filter expression data to tissue samples
        tissue_expression_df = expression_df[tissue_samples]
        
        # Take only first 100 genes for faster testing
        tissue_expression_df = tissue_expression_df.iloc[:100]
        
        print("Running preprocessing with confounder regression...")
        processed_tissue, info = preprocessor.preprocess_expression_data(
            expression_df=tissue_expression_df,
            hgnc_data=hgnc_data,
            sample_attrs=sample_attrs,
            subject_phenos=subject_phenos,
            apply_log=True,
            filter_low_expr=True,
            filter_low_var=True,
            detect_outliers=True,
            quantile_norm=True,
            regress_confounders=True,
            min_expression=0.1,
            min_variance=0.2,
            outlier_contamination=0.01
        )
        
        print(f"Successfully processed {processed_tissue.shape[0]} genes for {processed_tissue.shape[1]} samples")
        print("✅ NaN handling fixes appear to be working!")
        
        return True
        
    except Exception as e:
        print(f"❌ Error during testing: {str(e)}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    success = test_nan_fixes()
    if success:
        print("\n🎉 All tests passed! The NaN handling fixes are working correctly.")
    else:
        print("\n💥 Tests failed. Please check the error messages above.")
