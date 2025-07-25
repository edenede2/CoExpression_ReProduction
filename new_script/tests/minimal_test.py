"""
Minimal test to verify data loading works
"""

from data_reader import GTExDataReader
import traceback

def minimal_test():
    print("Minimal Data Loading Test")
    print("=" * 30)
    
    try:
        reader = GTExDataReader()
        print(f"Data directory: {reader.data_dir.absolute()}")
        print(f"Directory exists: {reader.data_dir.exists()}")
        
        if reader.data_dir.exists():
            files = list(reader.data_dir.glob("*.gct"))
            print(f"GCT files found: {len(files)}")
            for f in files:
                print(f"  - {f.name}")
        
        # Try to load just sample attributes (smaller file)
        print("\nTrying to load sample attributes...")
        sample_attrs = reader.read_sample_attributes()
        print(f"✅ Success! Loaded {sample_attrs.shape[0]} samples")
        
        # Try subject phenotypes
        print("\nTrying to load subject phenotypes...")
        subject_phenos = reader.read_subject_phenotypes()
        print(f"✅ Success! Loaded {subject_phenos.shape[0]} subjects")
        
        return True
        
    except Exception as e:
        print(f"❌ Error: {e}")
        print("\nFull traceback:")
        traceback.print_exc()
        return False

if __name__ == "__main__":
    success = minimal_test()
    print(f"\nTest result: {'✅ PASSED' if success else '❌ FAILED'}")
