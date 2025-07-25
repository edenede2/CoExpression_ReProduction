#!/usr/bin/env python3
"""
Quick verification script to test the inter-tissue network analysis setup
"""

import os
import sys

def main():
    print("🔍 Verifying Inter-tissue Network Analysis Setup")
    print("=" * 50)
    
    # Check if we're in the right directory
    current_dir = os.getcwd()
    print(f"Current directory: {current_dir}")
    
    # Check if required files exist
    required_files = [
        "new_script/inter_tissue_network_analysis.py",
        "new_script/inter_tissue_mdc_analysis.py", 
        "new_script/streamlit_app.py",
        "data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct",
        "data/GTEx_Analysis_v8_Annotations_SampleAttributesDS.tsv"
    ]
    
    print("\n📁 Checking required files:")
    all_files_exist = True
    for file_path in required_files:
        exists = os.path.exists(file_path)
        status = "✅" if exists else "❌"
        print(f"  {status} {file_path}")
        if not exists:
            all_files_exist = False
    
    if not all_files_exist:
        print("\n⚠️  Some required files are missing!")
        return False
    
    # Test imports
    print("\n📦 Testing imports:")
    try:
        # Add current directory to path
        sys.path.insert(0, '.')
        
        # Test core imports
        import pandas as pd
        import numpy as np
        print("  ✅ pandas and numpy")
        
        import networkx as nx
        print("  ✅ networkx")
        
        import sklearn
        print("  ✅ sklearn")
        
        # Test custom modules
        from new_script.inter_tissue_network_analysis import InterTissueNetworkAnalysis
        print("  ✅ InterTissueNetworkAnalysis")
        
        from new_script.inter_tissue_mdc_analysis import InterTissueMDCAnalysis
        print("  ✅ InterTissueMDCAnalysis")
        
        from new_script.network_analysis import NetworkAnalysis
        print("  ✅ NetworkAnalysis")
        
        print("\n🎉 All imports successful!")
        
    except ImportError as e:
        print(f"\n❌ Import error: {e}")
        return False
    except Exception as e:
        print(f"\n❌ Unexpected error: {e}")
        return False
    
    # Test basic functionality
    print("\n⚙️  Testing basic functionality:")
    try:
        # Initialize inter-tissue analyzer
        analyzer = InterTissueNetworkAnalysis()
        print("  ✅ InterTissueNetworkAnalysis initialization")
        
        # Test data file reading capability
        sample_metadata_path = "data/GTEx_Analysis_v8_Annotations_SampleAttributesDS.tsv"
        if os.path.exists(sample_metadata_path):
            # Read just first few lines to test
            with open(sample_metadata_path, 'r') as f:
                header = f.readline()
                if 'SAMPID' in header:
                    print("  ✅ GTEx sample metadata format verified")
                else:
                    print("  ⚠️  GTEx sample metadata format unexpected")
        
        print("\n🎯 Setup verification complete!")
        print("\nNext steps:")
        print("1. Run: streamlit run new_script/streamlit_app.py")
        print("2. Select 'Inter-tissue Analysis' mode")
        print("3. Upload your GTEx data and test the pipeline")
        
        return True
        
    except Exception as e:
        print(f"\n❌ Functionality test error: {e}")
        return False

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
