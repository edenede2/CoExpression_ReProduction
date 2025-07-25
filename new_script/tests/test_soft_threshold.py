#!/usr/bin/env python3
"""
Quick test for soft threshold power determination fix.
"""

import numpy as np
import pandas as pd
from network_analysis import NetworkAnalysis

def test_soft_threshold():
    """Test the soft threshold power determination with synthetic data."""
    print("Testing soft threshold power determination...")
    
    # Create synthetic correlation matrix (100 genes)
    np.random.seed(42)
    n_genes = 100
    
    # Generate more realistic gene expression correlation structure
    # Create hub genes with many connections (scale-free-like)
    base_corr = np.random.rand(n_genes, n_genes) * 0.2
    base_corr = (base_corr + base_corr.T) / 2
    
    # Create hub structure - some genes are highly connected
    n_hubs = 10
    hub_indices = np.random.choice(n_genes, n_hubs, replace=False)
    
    # Hubs have stronger correlations
    for hub in hub_indices:
        # Connect hub to many other genes
        n_connections = np.random.randint(20, 40)
        connected_genes = np.random.choice(n_genes, n_connections, replace=False)
        for gene in connected_genes:
            if gene != hub:
                base_corr[hub, gene] = base_corr[gene, hub] = np.random.uniform(0.6, 0.9)
    
    # Add some module structure
    for i in range(0, n_genes, 25):
        end = min(i + 25, n_genes)
        module_corr = np.random.uniform(0.4, 0.7)
        base_corr[i:end, i:end] = module_corr
    
    # Ensure diagonal is 1
    np.fill_diagonal(base_corr, 1.0)
    
    # Clip to valid correlation range
    corr_matrix = np.clip(base_corr, -1, 1)
    
    # Convert to DataFrame
    gene_names = [f"Gene_{i:03d}" for i in range(n_genes)]
    corr_df = pd.DataFrame(corr_matrix, index=gene_names, columns=gene_names)
    
    print(f"Created synthetic correlation matrix: {corr_df.shape}")
    print(f"Correlation range: {corr_df.values.min():.3f} to {corr_df.values.max():.3f}")
    
    # Test soft threshold determination
    network_analyzer = NetworkAnalysis()
    
    # Test with a smaller range of powers for speed
    test_powers = list(range(1, 11))
    
    try:
        optimal_power, results_df = network_analyzer.determine_soft_threshold_power(
            corr_df, 
            powers=test_powers,
            target_rsq=0.8
        )
        
        print(f"\nResults:")
        print(f"Optimal power: {optimal_power}")
        print(f"\nDetailed results:")
        print(results_df[['power', 'rsq', 'slope', 'mean_connectivity']].round(3))
        
        # Check if we got reasonable results
        max_rsq = results_df['rsq'].max()
        if max_rsq > 0.1:  # Should be better than 0.1
            print(f"\n✅ Success! Maximum R² achieved: {max_rsq:.3f}")
        else:
            print(f"\n❌ Low R² values. Maximum R²: {max_rsq:.3f}")
            
    except Exception as e:
        print(f"❌ Error during testing: {e}")
        raise

if __name__ == "__main__":
    test_soft_threshold()
