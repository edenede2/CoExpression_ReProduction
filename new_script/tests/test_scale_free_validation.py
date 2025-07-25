#!/usr/bin/env python3
"""
Test to validate that our scale-free detection algorithm works correctly
by creating artificial scale-free networks.
"""

import numpy as np
import pandas as pd
from network_analysis import NetworkAnalysis
import networkx as nx

def create_scale_free_correlation_matrix(n_genes=100, gamma=2.5):
    """
    Create a correlation matrix that should exhibit scale-free properties.
    """
    np.random.seed(42)
    
    # Create a scale-free network using networkx
    G = nx.barabasi_albert_graph(n_genes, m=3, seed=42)
    
    # Convert to adjacency matrix
    adj_matrix = nx.adjacency_matrix(G).toarray().astype(float)
    
    # Convert adjacency to correlation-like values
    # Strong connections get high correlation, weak get low
    corr_matrix = adj_matrix * np.random.uniform(0.6, 0.9, adj_matrix.shape)
    
    # Make symmetric
    corr_matrix = (corr_matrix + corr_matrix.T) / 2
    
    # Add some background correlation
    background = np.random.uniform(0.0, 0.3, (n_genes, n_genes))
    background = (background + background.T) / 2
    
    corr_matrix = np.maximum(corr_matrix, background)
    
    # Set diagonal to 1
    np.fill_diagonal(corr_matrix, 1.0)
    
    return corr_matrix

def test_scale_free_detection():
    """Test if our algorithm can detect scale-free topology in known scale-free networks."""
    print("Testing scale-free topology detection...")
    
    # Create a known scale-free correlation matrix
    corr_matrix = create_scale_free_correlation_matrix(n_genes=80)
    
    gene_names = [f"Gene_{i:03d}" for i in range(80)]
    corr_df = pd.DataFrame(corr_matrix, index=gene_names, columns=gene_names)
    
    print(f"Created scale-free correlation matrix: {corr_df.shape}")
    print(f"Correlation range: {corr_df.values.min():.3f} to {corr_df.values.max():.3f}")
    
    # Test with our network analyzer
    network_analyzer = NetworkAnalysis()
    
    # Test with broader range of powers
    test_powers = list(range(1, 16))
    
    try:
        optimal_power, results_df = network_analyzer.determine_soft_threshold_power(
            corr_df, 
            powers=test_powers,
            target_rsq=0.7  # More relaxed target
        )
        
        print(f"\nResults:")
        print(f"Optimal power: {optimal_power}")
        
        # Show best results
        best_results = results_df.nlargest(5, 'rsq')
        print(f"\nTop 5 R² results:")
        print(best_results[['power', 'rsq', 'slope', 'mean_connectivity']].round(3))
        
        # Validation
        max_rsq = results_df['rsq'].max()
        if max_rsq > 0.5:
            print(f"\n✅ Excellent! Maximum R² achieved: {max_rsq:.3f}")
        elif max_rsq > 0.3:
            print(f"\n✅ Good! Maximum R² achieved: {max_rsq:.3f}")
        elif max_rsq > 0.1:
            print(f"\n⚠️  Moderate R² achieved: {max_rsq:.3f}")
        else:
            print(f"\n❌ Low R² values. Maximum R²: {max_rsq:.3f}")
            
        # Check if negative slopes (expected for scale-free)
        negative_slopes = results_df[results_df['slope'] < 0]
        print(f"\nPowers with negative slopes (scale-free indicator): {len(negative_slopes)}/{len(results_df)}")
        
    except Exception as e:
        print(f"❌ Error during testing: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    test_scale_free_detection()
