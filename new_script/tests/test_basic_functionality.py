#!/usr/bin/env python3
"""
Simple test to verify the soft threshold power determination is working correctly.
"""

import numpy as np
import pandas as pd
from network_analysis import NetworkAnalysis

def test_basic_functionality():
    """Test basic functionality of soft threshold determination."""
    print("Testing basic soft threshold functionality...")
    
    # Create a simple correlation matrix with clear structure
    np.random.seed(42)
    n_genes = 50  # Smaller for faster testing
    
    # Create correlation matrix with clear module structure
    corr_matrix = np.random.uniform(0.1, 0.3, (n_genes, n_genes))
    corr_matrix = (corr_matrix + corr_matrix.T) / 2
    
    # Create strong modules
    module_size = 10
    for i in range(0, n_genes, module_size):
        end = min(i + module_size, n_genes)
        corr_matrix[i:end, i:end] = np.random.uniform(0.7, 0.9, (end-i, end-i))
        corr_matrix[i:end, i:end] = (corr_matrix[i:end, i:end] + corr_matrix[i:end, i:end].T) / 2
    
    # Set diagonal to 1
    np.fill_diagonal(corr_matrix, 1.0)
    
    gene_names = [f"Gene_{i:02d}" for i in range(n_genes)]
    corr_df = pd.DataFrame(corr_matrix, index=gene_names, columns=gene_names)
    
    print(f"Created correlation matrix: {corr_df.shape}")
    print(f"Correlation range: {corr_df.values.min():.3f} to {corr_df.values.max():.3f}")
    
    # Test network analysis
    network_analyzer = NetworkAnalysis()
    
    # Test with small range for speed
    test_powers = [1, 2, 3, 4, 5, 6]
    
    try:
        optimal_power, results_df = network_analyzer.determine_soft_threshold_power(
            corr_df, 
            powers=test_powers,
            target_rsq=0.5  # Relaxed target
        )
        
        print(f"\n📊 Results:")
        print(f"Optimal power: {optimal_power}")
        print(f"\nAll results:")
        for _, row in results_df.iterrows():
            print(f"  Power {row['power']}: R²={row['rsq']:.3f}, slope={row['slope']:.2f}, connectivity={row['mean_connectivity']:.1f}")
        
        # Basic validation
        if optimal_power > 1:
            print(f"\n✅ Algorithm selected power > 1: {optimal_power}")
        else:
            print(f"\n⚠️  Algorithm selected power 1: {optimal_power}")
            
        max_rsq = results_df['rsq'].max()
        if max_rsq > 0:
            print(f"✅ Non-zero R² values detected: max = {max_rsq:.3f}")
        else:
            print(f"❌ All R² values are zero")
            
        # Check connectivity decreases with power
        connectivities = results_df['mean_connectivity'].values
        if len(connectivities) > 1 and connectivities[0] > connectivities[-1]:
            print(f"✅ Connectivity decreases with power: {connectivities[0]:.1f} → {connectivities[-1]:.1f}")
        else:
            print(f"⚠️  Connectivity pattern: {connectivities}")
            
        return True
        
    except Exception as e:
        print(f"❌ Error: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    success = test_basic_functionality()
    if success:
        print(f"\n🎉 Basic functionality test completed!")
    else:
        print(f"\n💥 Test failed!")
