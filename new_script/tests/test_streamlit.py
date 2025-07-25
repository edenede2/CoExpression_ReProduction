"""
Test script to verify the Streamlit app functionality
"""

import sys
import os
sys.path.append('.')

from data_reader import GTExDataReader
import pandas as pd

def test_streamlit_compatibility():
    """Test if the app can load data and filter samples correctly."""
    
    print("Testing Streamlit App Compatibility")
    print("=" * 40)
    
    try:
        # Initialize data reader
        reader = GTExDataReader()
        
        # Load small subset of expression data (like the app does)
        expression_df, gene_metadata = reader.read_gct_file(
            "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct",
            sample_limit=500  # Small subset for testing
        )
        
        sample_attrs = reader.read_sample_attributes()
        subject_phenos = reader.read_subject_phenotypes()
        
        print(f"✅ Loaded {expression_df.shape[1]} samples, {expression_df.shape[0]} genes")
        
        # Get available samples from loaded expression data
        available_samples = set(expression_df.columns)
        print(f"Available samples in expression data: {len(available_samples)}")
        
        # Test tissue filtering with sample intersection
        available_tissues = reader.get_available_tissues(sample_attrs)
        print(f"Total tissues available: {len(available_tissues)}")
        
        # Find tissues with good sample counts in loaded data
        good_tissues = []
        for tissue in available_tissues[:10]:  # Check first 10 tissues
            all_tissue_samples = reader.get_tissue_samples(sample_attrs, tissue)
            loaded_tissue_samples = [s for s in all_tissue_samples if s in available_samples]
            if len(loaded_tissue_samples) >= 10:
                good_tissues.append((tissue, len(loaded_tissue_samples)))
        
        print(f"Tissues with ≥10 samples in loaded data:")
        for tissue, count in good_tissues[:5]:
            print(f"  {tissue}: {count} samples")
        
        if good_tissues:
            # Test age filtering with the best tissue
            test_tissue = good_tissues[0][0]
            print(f"\nTesting age filtering with {test_tissue}...")
            
            # Get young samples
            young_samples_all = reader.filter_samples_by_metadata(
                sample_attrs, subject_phenos,
                age_group='young', tissue=test_tissue
            )
            young_samples = [s for s in young_samples_all if s in available_samples]
            
            # Get old samples
            old_samples_all = reader.filter_samples_by_metadata(
                sample_attrs, subject_phenos,
                age_group='old', tissue=test_tissue
            )
            old_samples = [s for s in old_samples_all if s in available_samples]
            
            print(f"Young samples found: {len(young_samples)}")
            print(f"Old samples found: {len(old_samples)}")
            
            if len(young_samples) >= 3 and len(old_samples) >= 3:
                # Test expression data access
                all_samples = young_samples + old_samples
                filtered_expression = expression_df[all_samples]
                print(f"✅ Successfully filtered expression data: {filtered_expression.shape}")
                print("🎉 Streamlit app should work correctly!")
                return True
            else:
                print("⚠️  Not enough samples for comparison in this tissue")
        else:
            print("⚠️  No tissues with sufficient samples found")
            
        print("💡 Tip: Try increasing sample_limit in load_data() function")
        return False
        
    except Exception as e:
        print(f"❌ Error: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    success = test_streamlit_compatibility()
    if success:
        print(f"\n🚀 Ready to run: streamlit run streamlit_app.py")
    else:
        print(f"\n🔧 App may have issues - check the problems above")
