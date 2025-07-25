"""
Modular Differential Connectivity (MDC) Analysis Module
Implements MDC calculations for comparing network connectivity between conditions
"""

import pandas as pd
import numpy as np
from scipy import stats
from scipy.stats import chi2
from typing import Dict, List, Tuple, Optional
import warnings

class MDCAnalyzer:
    """Class for Modular Differential Connectivity analysis."""
    
    def __init__(self):
        """Initialize the MDC analyzer."""
        self.mdc_results = None
        
    def calculate_module_connectivity(self, 
                                      adjacency_matrix: pd.DataFrame,
                                      modules: Dict[str, int],
                                      module_id: int) -> Dict[str, float]:
        """
        Calculate connectivity metrics for a specific module.
        
        Args:
            adjacency_matrix: Gene adjacency matrix
            modules: Gene-to-module mapping
            module_id: Module ID to analyze
            
        Returns:
            Dictionary with connectivity metrics
        """
        # Get genes in this module
        module_genes = [gene for gene, mod in modules.items() if mod == module_id]
        
        if len(module_genes) < 2:
            return {'intra_connectivity': 0, 'total_connectivity': 0, 'module_size': 0}
        
        # Filter adjacency matrix for module genes
        module_adj = adjacency_matrix.loc[module_genes, module_genes]
        full_adj = adjacency_matrix.loc[module_genes]
        
        # Calculate intra-module connectivity (within module)
        # Sum of all connections within module (excluding diagonal)
        intra_conn = (module_adj.values.sum() - np.diag(module_adj).sum()) / 2
        
        # Calculate total connectivity (all connections from module genes)
        total_conn = (full_adj.values.sum() - np.diag(module_adj).sum()) / 2
        
        # Normalize by number of possible connections
        n_genes = len(module_genes)
        max_intra_conn = n_genes * (n_genes - 1) / 2
        
        return {
            'intra_connectivity': intra_conn,
            'total_connectivity': total_conn,
            'module_size': n_genes,
            'max_intra_connectivity': max_intra_conn,
            'intra_density': intra_conn / max_intra_conn if max_intra_conn > 0 else 0
        }
    
    def calculate_mdc_for_modules(self,
                                  adjacency_matrix_1: pd.DataFrame,
                                  adjacency_matrix_2: pd.DataFrame,
                                  modules: Dict[str, int],
                                  condition_1_name: str = "Condition1",
                                  condition_2_name: str = "Condition2") -> pd.DataFrame:
        """
        Calculate MDC for all modules comparing two conditions.
        
        Args:
            adjacency_matrix_1: Adjacency matrix for condition 1
            adjacency_matrix_2: Adjacency matrix for condition 2
            modules: Gene-to-module mapping
            condition_1_name: Name for condition 1
            condition_2_name: Name for condition 2
            
        Returns:
            DataFrame with MDC results for each module
        """
        print(f"Calculating MDC between {condition_1_name} and {condition_2_name}...")
        
        results = []
        unique_modules = set(modules.values())
        
        for module_id in unique_modules:
            if module_id == 0:  # Skip grey module
                continue
                
            print(f"Processing module {module_id}...")
            
            # Calculate connectivity for both conditions
            conn_1 = self.calculate_module_connectivity(adjacency_matrix_1, modules, module_id)
            conn_2 = self.calculate_module_connectivity(adjacency_matrix_2, modules, module_id)
            
            if conn_1['module_size'] < 3 or conn_2['module_size'] < 3:
                continue
            
            # Calculate MDC ratio
            if conn_1['intra_connectivity'] > 0 and conn_2['intra_connectivity'] > 0:
                mdc_ratio = conn_2['intra_connectivity'] / conn_1['intra_connectivity']
                log_mdc_ratio = np.log2(mdc_ratio)
            else:
                mdc_ratio = np.nan
                log_mdc_ratio = np.nan
            
            # Calculate density difference
            density_diff = conn_2['intra_density'] - conn_1['intra_density']
            
            results.append({
                'module_id': module_id,
                'module_size': conn_1['module_size'],
                f'{condition_1_name}_connectivity': conn_1['intra_connectivity'],
                f'{condition_2_name}_connectivity': conn_2['intra_connectivity'],
                f'{condition_1_name}_density': conn_1['intra_density'],
                f'{condition_2_name}_density': conn_2['intra_density'],
                'mdc_ratio': mdc_ratio,
                'log_mdc_ratio': log_mdc_ratio,
                'density_difference': density_diff
            })
        
        mdc_df = pd.DataFrame(results)
        
        print(f"Calculated MDC for {len(mdc_df)} modules")
        return mdc_df
    
    def permutation_test_mdc(self,
                            expression_df_1: pd.DataFrame,
                            expression_df_2: pd.DataFrame,
                            modules: Dict[str, int],
                            power: int,
                            n_permutations: int = 100,
                            correlation_method: str = 'pearson') -> pd.DataFrame:
        """
        Perform permutation test for MDC significance.
        
        Args:
            expression_df_1: Expression data for condition 1
            expression_df_2: Expression data for condition 2
            modules: Gene-to-module mapping
            power: Soft threshold power
            n_permutations: Number of permutations
            correlation_method: Correlation method
            
        Returns:
            DataFrame with permutation test results
        """
        print(f"Performing permutation test with {n_permutations} permutations...")
        
        from network_analysis import NetworkAnalysis
        
        # Combine expression data
        combined_expr = pd.concat([expression_df_1, expression_df_2], axis=1)
        n_samples_1 = expression_df_1.shape[1]
        n_samples_2 = expression_df_2.shape[1]
        total_samples = n_samples_1 + n_samples_2
        
        # Calculate observed MDC
        print("Calculating observed MDC...")
        network_analyzer = NetworkAnalysis()
        
        # Build networks for both conditions
        corr_1 = network_analyzer.calculate_correlation_matrix(expression_df_1, correlation_method)
        corr_2 = network_analyzer.calculate_correlation_matrix(expression_df_2, correlation_method)
        
        adj_1 = network_analyzer._calculate_adjacency_matrix(corr_1, power)
        adj_2 = network_analyzer._calculate_adjacency_matrix(corr_2, power)
        
        observed_mdc = self.calculate_mdc_for_modules(adj_1, adj_2, modules)
        
        # Permutation testing
        print("Running permutations...")
        permuted_mdc_ratios = {module_id: [] for module_id in observed_mdc['module_id']}
        
        for perm in range(n_permutations):
            if perm % 20 == 0:
                print(f"Permutation {perm}/{n_permutations}")
            
            # Randomly shuffle sample assignments
            shuffled_indices = np.random.permutation(total_samples)
            perm_samples_1 = shuffled_indices[:n_samples_1]
            perm_samples_2 = shuffled_indices[n_samples_1:]
            
            # Create permuted expression matrices
            perm_expr_1 = combined_expr.iloc[:, perm_samples_1]
            perm_expr_2 = combined_expr.iloc[:, perm_samples_2]
            
            # Calculate permuted networks
            perm_corr_1 = network_analyzer.calculate_correlation_matrix(perm_expr_1, correlation_method)
            perm_corr_2 = network_analyzer.calculate_correlation_matrix(perm_expr_2, correlation_method)
            
            perm_adj_1 = network_analyzer._calculate_adjacency_matrix(perm_corr_1, power)
            perm_adj_2 = network_analyzer._calculate_adjacency_matrix(perm_corr_2, power)
            
            # Calculate permuted MDC
            perm_mdc = self.calculate_mdc_for_modules(perm_adj_1, perm_adj_2, modules)
            
            # Store permuted ratios
            for _, row in perm_mdc.iterrows():
                module_id = row['module_id']
                if not np.isnan(row['mdc_ratio']):
                    permuted_mdc_ratios[module_id].append(row['log_mdc_ratio'])
        
        # Calculate p-values
        results = []
        for _, row in observed_mdc.iterrows():
            module_id = row['module_id']
            observed_log_ratio = row['log_mdc_ratio']
            
            if np.isnan(observed_log_ratio) or len(permuted_mdc_ratios[module_id]) == 0:
                p_value = np.nan
            else:
                # Two-tailed test
                permuted_ratios = np.array(permuted_mdc_ratios[module_id])
                p_value = np.mean(np.abs(permuted_ratios) >= np.abs(observed_log_ratio))
            
            results.append({
                'module_id': module_id,
                'module_size': row['module_size'],
                'mdc_ratio': row['mdc_ratio'],
                'log_mdc_ratio': row['log_mdc_ratio'],
                'p_value': p_value,
                'n_permutations': len(permuted_mdc_ratios[module_id])
            })
        
        results_df = pd.DataFrame(results)
        
        # Multiple testing correction (FDR)
        valid_p_values = results_df['p_value'].dropna()
        if len(valid_p_values) > 0:
            from statsmodels.stats.multitest import multipletests
            reject, corrected_p, _, _ = multipletests(valid_p_values, method='fdr_bh')
            
            results_df['p_value_corrected'] = np.nan
            results_df.loc[valid_p_values.index, 'p_value_corrected'] = corrected_p
            results_df.loc[valid_p_values.index, 'significant'] = reject
        else:
            results_df['p_value_corrected'] = np.nan
            results_df['significant'] = False
        
        print("Permutation test completed")
        return results_df
    
    def analyze_differential_connectivity(self,
                                          expression_df_1: pd.DataFrame,
                                          expression_df_2: pd.DataFrame,
                                          modules: Dict[str, int],
                                          power: int,
                                          condition_1_name: str = "Condition1",
                                          condition_2_name: str = "Condition2",
                                          run_permutation_test: bool = True,
                                          n_permutations: int = 100) -> Dict:
        """
        Complete differential connectivity analysis pipeline.
        
        Args:
            expression_df_1: Expression data for condition 1
            expression_df_2: Expression data for condition 2
            modules: Gene-to-module mapping
            power: Soft threshold power
            condition_1_name: Name for condition 1
            condition_2_name: Name for condition 2
            run_permutation_test: Whether to run permutation test
            n_permutations: Number of permutations
            
        Returns:
            Dictionary with analysis results
        """
        print("Starting differential connectivity analysis...")
        
        from network_analysis import NetworkAnalysis
        
        # Build networks for both conditions
        network_analyzer = NetworkAnalysis()
        
        print(f"Building network for {condition_1_name}...")
        network_1 = network_analyzer.construct_network(
            expression_df_1,
            power=power, 
            min_module_size=1  # Use existing modules
        )
        
        print(f"Building network for {condition_2_name}...")
        network_2 = network_analyzer.construct_network(
            expression_df_2, 
            power=power, 
            min_module_size=1  # Use existing modules
        )
        
        # Calculate MDC
        mdc_results = self.calculate_mdc_for_modules(
            network_1['adjacency_matrix'],
            network_2['adjacency_matrix'],
            modules,
            condition_1_name,
            condition_2_name
        )
        
        results = {
            'mdc_results': mdc_results,
            'network_1': network_1,
            'network_2': network_2,
            'condition_1_name': condition_1_name,
            'condition_2_name': condition_2_name
        }
        
        # Permutation test
        if run_permutation_test:
            permutation_results = self.permutation_test_mdc(
                expression_df_1,
                expression_df_2,
                modules,
                power,
                n_permutations
            )
            results['permutation_results'] = permutation_results
        
        self.mdc_results = results
        
        print("Differential connectivity analysis completed!")
        return results
    
    def get_significant_modules(self, 
                                p_threshold: float = 0.05,
                                log_ratio_threshold: float = 0.5) -> pd.DataFrame:
        """
        Get modules with significant differential connectivity.
        
        Args:
            p_threshold: P-value threshold for significance
            log_ratio_threshold: Minimum absolute log ratio threshold
            
        Returns:
            DataFrame with significant modules
        """
        if self.mdc_results is None or 'permutation_results' not in self.mdc_results:
            print("No permutation test results available")
            return pd.DataFrame()
        
        perm_results = self.mdc_results['permutation_results']
        
        # Filter significant modules
        significant = perm_results[
            (perm_results['p_value_corrected'] <= p_threshold) &
            (np.abs(perm_results['log_mdc_ratio']) >= log_ratio_threshold)
        ].copy()
        
        # Sort by absolute log ratio
        significant = significant.sort_values('log_mdc_ratio', key=abs, ascending=False)
        
        print(f"Found {len(significant)} significant modules")
        return significant
