"""
Quick test script to verify age group filtering works correctly
"""

from data_reader import GTExDataReader
import pandas as pd

def test_age_filtering():
    """Test the corrected age group filtering functionality."""
    
    print("Testing Age Group Filtering")
    print("=" * 40)
    
    reader = GTExDataReader()
    
    try:
        # Load metadata
        sample_attrs = reader.read_sample_attributes()
        subject_phenos = reader.read_subject_phenotypes()
        
        print(f"Loaded {len(subject_phenos)} subjects")
        
        # Show available age groups
        age_groups = reader.get_available_age_groups(subject_phenos)
        age_counts = reader.get_age_group_counts(subject_phenos)
        
        print("\nAvailable age groups:")
        for age in age_groups:
            print(f"  {age}: {age_counts.get(age, 0)} subjects")
        
        # Show available sexes
        sexes = reader.get_available_sexes(subject_phenos)
        print(f"\nAvailable sexes: {sexes}")
        
        # Test young group filtering
        print("\nTesting 'young' group filtering...")
        young_samples = reader.filter_samples_by_metadata(
            sample_attrs, subject_phenos, age_group='young'
        )
        print(f"Found {len(young_samples)} young samples")
        
        # Test old group filtering  
        print("\nTesting 'old' group filtering...")
        old_samples = reader.filter_samples_by_metadata(
            sample_attrs, subject_phenos, age_group='old'
        )
        print(f"Found {len(old_samples)} old samples")
        
        # Test middle group filtering
        print("\nTesting 'middle' group filtering...")
        middle_samples = reader.filter_samples_by_metadata(
            sample_attrs, subject_phenos, age_group='middle'
        )
        print(f"Found {len(middle_samples)} middle-aged samples")
        
        # Test custom age group filtering
        print("\nTesting custom age group filtering ['20-29', '70-79']...")
        custom_samples = reader.filter_samples_by_metadata(
            sample_attrs, subject_phenos, age_group=['20-29', '70-79']
        )
        print(f"Found {len(custom_samples)} samples in custom age groups")
        
        # Test with tissue filtering
        tissues = reader.get_available_tissues(sample_attrs)
        if tissues:
            test_tissue = tissues[0]
            print(f"\nTesting with tissue filter ({test_tissue})...")
            
            tissue_young = reader.filter_samples_by_metadata(
                sample_attrs, subject_phenos, 
                age_group='young', tissue=test_tissue
            )
            
            tissue_old = reader.filter_samples_by_metadata(
                sample_attrs, subject_phenos, 
                age_group='old', tissue=test_tissue
            )
            
            print(f"  Young samples in {test_tissue}: {len(tissue_young)}")
            print(f"  Old samples in {test_tissue}: {len(tissue_old)}")
        
        print("\n✅ Age filtering test completed successfully!")
        
        return True
        
    except Exception as e:
        print(f"❌ Test failed: {e}")
        return False

if __name__ == "__main__":
    test_age_filtering()
