#!/usr/bin/env python3
"""
Test inter-tissue network analysis functionality
"""

import sys
import os
sys.path.append('/Users/edeneldar/CoExpression_reProduction/new_script')

def test_imports():
    """Test that all modules can be imported successfully."""
    print("Testing imports...")
    
    try:
        from network_analysis import NetworkAnalysis
        print("✅ NetworkAnalysis imported successfully")
        
        from inter_tissue_network_analysis import InterTissueNetworkAnalysis
        print("✅ InterTissueNetworkAnalysis imported successfully")
        
        from inter_tissue_mdc_analysis import InterTissueMDCAnalysis
        print("✅ InterTissueMDCAnalysis imported successfully")
        
        # Test initialization
        network_analyzer = NetworkAnalysis()
        print("✅ NetworkAnalysis initialized successfully")
        
        network_analyzer.set_analysis_type('inter_tissue')
        print("✅ Inter-tissue analysis type set successfully")
        
        network_analyzer.set_analysis_type('traditional')
        print("✅ Traditional analysis type set successfully")
        
        return True
        
    except Exception as e:
        print(f"❌ Import error: {e}")
        return False

def test_functionality():
    """Test basic functionality."""
    print("\nTesting basic functionality...")
    
    try:
        import numpy as np
        import pandas as pd
        from network_analysis import NetworkAnalysis
        
        # Create simple test data
        n_samples = 50
        n_genes = 100
        
        # Generate synthetic expression data
        np.random.seed(42)
        expression_data = np.random.randn(n_genes, n_samples)
        
        gene_names = [f"Gene_{i:03d}" for i in range(n_genes)]
        sample_names = [f"Sample_{i:03d}" for i in range(n_samples)]
        
        expr_df = pd.DataFrame(expression_data, index=gene_names, columns=sample_names)
        
        print(f"Created test expression data: {expr_df.shape}")
        
        # Test traditional network analysis
        network_analyzer = NetworkAnalysis()
        network_analyzer.set_analysis_type('traditional')
        
        # Test correlation calculation
        corr_matrix = network_analyzer.calculate_correlation_matrix(expr_df)
        print(f"✅ Correlation matrix calculated: {corr_matrix.shape}")
        
        # Test power determination (with small range for speed)
        optimal_power, power_results = network_analyzer.determine_soft_threshold_power(
            corr_matrix, powers=[1, 2, 3], target_rsq=0.5
        )
        print(f"✅ Optimal power determined: {optimal_power}")
        
        return True
        
    except Exception as e:
        print(f"❌ Functionality test error: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    print("🧪 Testing Inter-tissue Network Analysis Implementation")
    print("=" * 60)
    
    # Test imports
    import_success = test_imports()
    
    if import_success:
        # Test functionality
        func_success = test_functionality()
        
        if func_success:
            print("\n🎉 All tests passed! Inter-tissue network analysis is ready.")
            print("\nNext steps:")
            print("1. Run the Streamlit app: streamlit run streamlit_app.py")
            print("2. Load your GTEx data")
            print("3. Try both Traditional and Inter-tissue analysis modes")
        else:
            print("\n❌ Functionality tests failed.")
    else:
        print("\n❌ Import tests failed.")
    
    print("\n" + "=" * 60)
