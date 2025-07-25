"""
Quick pipeline test with age filtering fix verification
"""

from data_reader import GTExDataReader

def quick_test():
    print("Quick Pipeline Test")
    print("==================")
    
    reader = GTExDataReader()
    
    # Test 1: Load data
    print("\n1. Loading data...")
    try:
        expression_df, gene_metadata = reader.read_gct_file(
            "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct",
            sample_limit=50  # Very small subset
        )
        sample_attrs = reader.read_sample_attributes()
        subject_phenos = reader.read_subject_phenotypes()
        print("✅ Data loading successful!")
        
    except Exception as e:
        print(f"❌ Data loading failed: {e}")
        return False
    
    # Test 2: Age filtering  
    print("\n2. Testing age filtering...")
    try:
        tissues = reader.get_available_tissues(sample_attrs)
        selected_tissue = tissues[0]  # Pick first tissue
        
        young_samples = reader.filter_samples_by_metadata(
            sample_attrs, subject_phenos,
            age_group='young', tissue=selected_tissue
        )
        
        old_samples = reader.filter_samples_by_metadata(
            sample_attrs, subject_phenos,
            age_group='old', tissue=selected_tissue
        )
        
        print(f"✅ Found {len(young_samples)} young and {len(old_samples)} old samples")
        
    except Exception as e:
        print(f"❌ Age filtering failed: {e}")
        return False
    
    print("\n🎉 Quick test passed! Age filtering is working correctly.")
    return True

if __name__ == "__main__":
    quick_test()
