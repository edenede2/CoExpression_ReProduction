"""
Inter-tissue Modular Differential Connectivity (MDC) Analysis
Implementation based on the original R script methodology.

MDC computes the ratio of mean intra-module connectivity in one network 
to that of the same set of genes in another network.
"""

import numpy as np
import pandas as pd
from typing import Dict, List, Tuple, Optional
from scipy import stats
import warnings

class InterTissueMDCAnalysis:
    """
    Inter-tissue Modular Differential Connectivity analysis for comparing networks between conditions.
    """
    
    def __init__(self):
        self.modules = None
        self.networks = {}
        self.mdc_results = None
    
    def load_modules(self, gene_modules: Dict[str, int]) -> None:
        """
        Load module assignments.
        
        Args:
            gene_modules: Dictionary mapping gene names to module IDs
        """
        self.modules = gene_modules
        print(f"Loaded {len(gene_modules)} gene module assignments")
        
        # Get module sizes
        module_sizes = {}
        for gene, module_id in gene_modules.items():
            if module_id != 0:  # Exclude grey module
                module_sizes[module_id] = module_sizes.get(module_id, 0) + 1
        
        print(f"Module sizes: {dict(sorted(module_sizes.items()))}")
    
    def load_condition_networks(self, 
                              condition_networks: Dict[str, pd.DataFrame]) -> None:
        """
        Load adjacency matrices for different conditions.
        
        Args:
            condition_networks: Dict mapping condition names to adjacency matrices
        """
        self.networks = condition_networks
        print(f"Loaded networks for conditions: {list(condition_networks.keys())}")
        
        for condition, network in condition_networks.items():
            print(f"  {condition}: {network.shape}")
    
    def calculate_mdc(self, 
                     condition1: str, 
                     condition2: str,
                     min_module_size: int = 10,
                     n_permutations: int = 1000) -> pd.DataFrame:
        """
        Calculate MDC between two conditions.
        
        Args:
            condition1: Name of first condition
            condition2: Name of second condition  
            min_module_size: Minimum module size to analyze
            n_permutations: Number of permutations for FDR calculation
            
        Returns:
            DataFrame with MDC results
        """
        if self.modules is None:
            raise ValueError("Modules not loaded")
        if condition1 not in self.networks or condition2 not in self.networks:
            raise ValueError(f"Networks not found for conditions: {condition1}, {condition2}")
        
        print(f"Calculating MDC: {condition1} vs {condition2}")
        
        network1 = self.networks[condition1]
        network2 = self.networks[condition2]
        
        # Get common genes
        common_genes = list(set(network1.index) & set(network2.index))
        if len(common_genes) == 0:
            raise ValueError("No common genes between networks")
        
        # Filter networks to common genes
        net1_filtered = network1.loc[common_genes, common_genes]
        net2_filtered = network2.loc[common_genes, common_genes]
        
        print(f"Analyzing {len(common_genes)} common genes")
        
        # Get modules for common genes
        gene_modules_filtered = {gene: self.modules[gene] for gene in common_genes 
                               if gene in self.modules}
        
        # Group genes by module
        modules_dict = {}
        for gene, module_id in gene_modules_filtered.items():
            if module_id != 0:  # Exclude grey module
                if module_id not in modules_dict:
                    modules_dict[module_id] = []
                modules_dict[module_id].append(gene)
        
        # Filter modules by size
        valid_modules = {mod_id: genes for mod_id, genes in modules_dict.items() 
                        if len(genes) >= min_module_size}
        
        print(f"Analyzing {len(valid_modules)} modules (min size: {min_module_size})")
        
        # Calculate MDC for each module
        mdc_results = []
        
        for module_id, module_genes in valid_modules.items():
            print(f"Processing module {module_id} ({len(module_genes)} genes)")
            
            # Calculate intra-module connectivity for both conditions
            conn1 = self._calculate_intra_module_connectivity(net1_filtered, module_genes)
            conn2 = self._calculate_intra_module_connectivity(net2_filtered, module_genes)
            
            if conn1 == 0 or conn2 == 0:
                print(f"  Warning: Zero connectivity in module {module_id}")
                continue
            
            # Calculate MDC
            mdc = conn1 / conn2
            
            # Calculate permutation-based FDR
            perm_mdcs = self._calculate_permutation_mdc(
                net1_filtered, net2_filtered, module_genes, n_permutations
            )
            
            fdr = self._calculate_fdr(mdc, perm_mdcs)
            
            mdc_results.append({
                'Module_ID': module_id,
                'Module_Size': len(module_genes),
                'Connectivity_Condition1': conn1,
                'Connectivity_Condition2': conn2,
                'MDC': mdc,
                'Permutation_Mean': np.mean(perm_mdcs),
                'Permutation_Std': np.std(perm_mdcs),
                'FDR': fdr,
                'Log2_MDC': np.log2(mdc) if mdc > 0 else np.nan,
                'Direction': 'Up' if mdc > 1 else 'Down'
            })
        
        self.mdc_results = pd.DataFrame(mdc_results)
        self.mdc_results = self.mdc_results.sort_values('FDR')
        
        print(f"MDC analysis completed for {len(self.mdc_results)} modules")
        return self.mdc_results
    
    def _calculate_intra_module_connectivity(self, 
                                           network: pd.DataFrame, 
                                           module_genes: List[str]) -> float:
        """
        Calculate mean intra-module connectivity.
        """
        # Get subnetwork for module genes
        module_network = network.loc[module_genes, module_genes]
        
        # Remove diagonal (self-connections)
        np.fill_diagonal(module_network.values, 0)
        
        # Calculate mean connectivity
        if module_network.shape[0] <= 1:
            return 0.0
        
        # Sum of all connections divided by number of possible connections
        total_connections = module_network.sum().sum()
        possible_connections = len(module_genes) * (len(module_genes) - 1)
        
        return total_connections / possible_connections if possible_connections > 0 else 0.0
    
    def _calculate_permutation_mdc(self, 
                                 network1: pd.DataFrame,
                                 network2: pd.DataFrame,
                                 module_genes: List[str],
                                 n_permutations: int) -> List[float]:
        """
        Calculate permutation-based MDC distribution.
        """
        all_genes = list(network1.index)
        module_size = len(module_genes)
        perm_mdcs = []
        
        for _ in range(n_permutations):
            # Randomly sample genes of same size
            random_genes = np.random.choice(all_genes, module_size, replace=False)
            
            # Calculate connectivity for random module
            conn1 = self._calculate_intra_module_connectivity(network1, random_genes)
            conn2 = self._calculate_intra_module_connectivity(network2, random_genes)
            
            if conn1 > 0 and conn2 > 0:
                perm_mdc = conn1 / conn2
                perm_mdcs.append(perm_mdc)
        
        return perm_mdcs
    
    def _calculate_fdr(self, observed_mdc: float, permutation_mdcs: List[float]) -> float:
        """
        Calculate False Discovery Rate based on permutation test.
        """
        if len(permutation_mdcs) == 0:
            return 1.0
        
        if observed_mdc < 1:
            # For MDC < 1, count how many permutations are less than observed
            fdr = sum(1 for perm_mdc in permutation_mdcs if perm_mdc < observed_mdc)
        else:
            # For MDC > 1, count how many permutations are greater than observed
            fdr = sum(1 for perm_mdc in permutation_mdcs if perm_mdc > observed_mdc)
        
        return fdr / len(permutation_mdcs)
    
    def get_significant_modules(self, fdr_threshold: float = 0.05) -> pd.DataFrame:
        """
        Get modules with significant MDC changes.
        
        Args:
            fdr_threshold: FDR threshold for significance
            
        Returns:
            DataFrame with significant modules
        """
        if self.mdc_results is None:
            raise ValueError("MDC results not calculated yet")
        
        significant = self.mdc_results[self.mdc_results['FDR'] <= fdr_threshold]
        
        print(f"Found {len(significant)} significant modules (FDR <= {fdr_threshold})")
        return significant
    
    def summarize_mdc_results(self) -> Dict:
        """
        Summarize MDC analysis results.
        """
        if self.mdc_results is None:
            raise ValueError("MDC results not calculated yet")
        
        summary = {
            'total_modules': len(self.mdc_results),
            'significant_modules_005': len(self.mdc_results[self.mdc_results['FDR'] <= 0.05]),
            'significant_modules_01': len(self.mdc_results[self.mdc_results['FDR'] <= 0.1]),
            'upregulated_modules': len(self.mdc_results[
                (self.mdc_results['MDC'] > 1) & (self.mdc_results['FDR'] <= 0.05)
            ]),
            'downregulated_modules': len(self.mdc_results[
                (self.mdc_results['MDC'] < 1) & (self.mdc_results['FDR'] <= 0.05)
            ]),
            'mean_mdc': self.mdc_results['MDC'].mean(),
            'median_mdc': self.mdc_results['MDC'].median()
        }
        
        return summary
    
    def export_results(self, output_dir: str) -> None:
        """
        Export MDC results to files.
        """
        if self.mdc_results is None:
            raise ValueError("MDC results not calculated yet")
        
        # Export full results
        self.mdc_results.to_csv(f"{output_dir}/mdc_results.csv", index=False)
        
        # Export significant results
        significant = self.get_significant_modules(0.05)
        significant.to_csv(f"{output_dir}/mdc_significant.csv", index=False)
        
        # Export summary
        summary = self.summarize_mdc_results()
        summary_df = pd.DataFrame([summary])
        summary_df.to_csv(f"{output_dir}/mdc_summary.csv", index=False)
        
        print(f"MDC results exported to {output_dir}")
