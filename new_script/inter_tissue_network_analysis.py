"""
Inter-tissue Coexpression Network Analysis
Implementation based on X-WGCNA methodology from original R scripts.

This module implements cross-tissue network analysis where:
- TS (Tissue-Specific) connections: within the same tissue
- CT (Cross-Tissue) connections: between different tissues
- Different soft threshold powers for TS vs CT connections
- Proper scale-free topology validation for inter-tissue networks
"""

import numpy as np
import pandas as pd
import scipy.stats as stats
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from sklearn.decomposition import PCA
from typing import Dict, List, Tuple, Optional, Union
import warnings
from itertools import combinations

class InterTissueNetworkAnalysis:
    """
    Inter-tissue coexpression network analysis using X-WGCNA methodology.
    """
    
    def __init__(self):
        self.tissue_data = {}
        self.tissue_indices = {}
        self.adjacency_matrix = None
        self.tom_matrix = None
        self.clusters = None
        
    def load_tissue_data(self, 
                        tissue_expression_data: Dict[str, pd.DataFrame],
                        tissue_sample_mapping: Dict[str, List[str]]) -> None:
        """
        Load expression data for multiple tissues.
        
        Args:
            tissue_expression_data: Dict mapping tissue names to expression matrices
                                  (samples x genes)
            tissue_sample_mapping: Dict mapping tissue names to sample IDs
        """
        print("Loading tissue expression data...")
        
        self.tissue_data = {}
        self.tissue_indices = {}
        
        # Track gene positions in the combined matrix
        current_position = 0
        all_gene_names = []
        
        for tissue_name, expr_data in tissue_expression_data.items():
            print(f"Processing {tissue_name}: {expr_data.shape}")
            
            # Store tissue data
            tissue_genes = [f"{tissue_name}_{gene}" for gene in expr_data.columns]
            self.tissue_data[tissue_name] = {
                'expression': expr_data,
                'genes': tissue_genes,
                'samples': tissue_sample_mapping[tissue_name]
            }
            
            # Track positions for adjacency matrix construction
            start_idx = current_position
            end_idx = current_position + len(tissue_genes)
            self.tissue_indices[tissue_name] = (start_idx, end_idx)
            current_position = end_idx
            
            all_gene_names.extend(tissue_genes)
        
        self.all_gene_names = all_gene_names
        print(f"Total genes across all tissues: {len(self.all_gene_names)}")
        
    def determine_optimal_powers(self,
                               powers_to_test: List[float] = None,
                               target_rsq: float = 0.8) -> Dict[str, Dict]:
        """
        Determine optimal soft threshold powers for tissue-specific and cross-tissue connections.
        
        Args:
            powers_to_test: List of powers to test
            target_rsq: Target R² for scale-free topology
            
        Returns:
            Dictionary with power analysis results
        """
        if powers_to_test is None:
            powers_to_test = np.arange(0.5, 20.5, 1.0)
        
        print("Determining optimal soft threshold powers...")
        power_results = {}
        
        # 1. Determine powers for individual tissues (TS connections)
        print("\n=== Tissue-Specific Power Analysis ===")
        ts_powers = {}
        
        for tissue_name, tissue_info in self.tissue_data.items():
            print(f"\nAnalyzing {tissue_name}...")
            expr_data = tissue_info['expression']
            
            power_analysis = self._analyze_power_for_expression(
                expr_data, powers_to_test, target_rsq
            )
            
            ts_powers[tissue_name] = power_analysis
            print(f"{tissue_name} optimal TS power: {power_analysis['optimal_power']} "
                  f"(R² = {power_analysis['optimal_rsq']:.3f})")
        
        # 2. Determine powers for tissue pairs (CT connections)
        print("\n=== Cross-Tissue Power Analysis ===")
        ct_powers = {}
        
        tissue_names = list(self.tissue_data.keys())
        for i, tissue1 in enumerate(tissue_names):
            for j, tissue2 in enumerate(tissue_names[i+1:], i+1):
                pair_name = f"{tissue1}-{tissue2}"
                print(f"\nAnalyzing cross-tissue pair: {pair_name}")
                
                # Find common samples
                samples1 = set(self.tissue_data[tissue1]['samples'])
                samples2 = set(self.tissue_data[tissue2]['samples'])
                common_samples = list(samples1.intersection(samples2))
                
                if len(common_samples) < 10:
                    print(f"Warning: Only {len(common_samples)} common samples for {pair_name}")
                    continue
                
                # Combine expression data for common samples
                expr1 = self.tissue_data[tissue1]['expression'].loc[common_samples]
                expr2 = self.tissue_data[tissue2]['expression'].loc[common_samples]
                combined_expr = pd.concat([expr1, expr2], axis=1)
                
                power_analysis = self._analyze_power_for_expression(
                    combined_expr, powers_to_test, target_rsq
                )
                
                ct_powers[pair_name] = power_analysis
                print(f"{pair_name} optimal CT power: {power_analysis['optimal_power']} "
                      f"(R² = {power_analysis['optimal_rsq']:.3f})")
        
        power_results = {
            'tissue_specific': ts_powers,
            'cross_tissue': ct_powers,
            'recommended_ts_power': self._get_recommended_power(ts_powers),
            'recommended_ct_power': self._get_recommended_power(ct_powers)
        }
        
        print(f"\n=== Final Recommendations ===")
        print(f"Recommended TS power: {power_results['recommended_ts_power']}")
        print(f"Recommended CT power: {power_results['recommended_ct_power']}")
        
        return power_results
    
    def _analyze_power_for_expression(self,
                                    expr_data: pd.DataFrame,
                                    powers_to_test: List[float],
                                    target_rsq: float) -> Dict:
        """
        Analyze scale-free topology for given expression data across powers.
        """
        results = []
        
        # Calculate correlation matrix
        corr_matrix = expr_data.T.corr(method='pearson')
        
        for power in powers_to_test:
            # Calculate adjacency matrix
            adj_matrix = np.abs(corr_matrix.values) ** power
            np.fill_diagonal(adj_matrix, 0)
            
            # Calculate connectivity
            connectivity = adj_matrix.sum(axis=1)
            
            # Scale-free topology analysis
            rsq, slope = self._calculate_scale_free_fit(connectivity)
            mean_connectivity = np.mean(connectivity)
            
            results.append({
                'power': power,
                'rsq': rsq,
                'slope': slope,
                'mean_connectivity': mean_connectivity
            })
        
        results_df = pd.DataFrame(results)
        
        # Find optimal power
        high_rsq = results_df[results_df['rsq'] >= target_rsq]
        if len(high_rsq) > 0:
            optimal_power = high_rsq.iloc[0]['power']
            optimal_rsq = high_rsq.iloc[0]['rsq']
        else:
            # Choose power with highest R²
            optimal_idx = results_df['rsq'].idxmax()
            optimal_power = results_df.loc[optimal_idx, 'power']
            optimal_rsq = results_df.loc[optimal_idx, 'rsq']
        
        return {
            'results': results_df,
            'optimal_power': optimal_power,
            'optimal_rsq': optimal_rsq,
            'target_rsq': target_rsq
        }
    
    def _calculate_scale_free_fit(self, connectivity: np.ndarray) -> Tuple[float, float]:
        """
        Calculate scale-free topology fit using WGCNA methodology.
        """
        # Remove zero connectivity
        nonzero_k = connectivity[connectivity > 0]
        
        if len(nonzero_k) < 10:
            return 0.0, 0.0
        
        # Bin connectivity values
        n_breaks = min(50, int(np.sqrt(len(nonzero_k))))
        hist, bin_edges = np.histogram(nonzero_k, bins=n_breaks)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        
        # Remove empty bins
        valid_bins = hist > 0
        if valid_bins.sum() < 5:
            return 0.0, 0.0
        
        k_values = bin_centers[valid_bins]
        p_k_values = hist[valid_bins] / len(nonzero_k)
        
        # Log-log regression
        try:
            log_k = np.log10(k_values)
            log_p_k = np.log10(p_k_values + 1e-9)
            
            slope, intercept, r_value, p_value, std_err = stats.linregress(log_k, log_p_k)
            rsq = r_value ** 2
            
            # For scale-free networks, we want negative slope
            if slope > 0:
                rsq = 0.0
                
            return rsq, slope
        except:
            return 0.0, 0.0
    
    def _get_recommended_power(self, power_dict: Dict) -> float:
        """
        Get recommended power from multiple power analyses.
        """
        powers = []
        rsqs = []
        
        for analysis in power_dict.values():
            if analysis['optimal_rsq'] > 0.5:  # Only consider good fits
                powers.append(analysis['optimal_power'])
                rsqs.append(analysis['optimal_rsq'])
        
        if not powers:
            # Fallback to any available power
            powers = [analysis['optimal_power'] for analysis in power_dict.values()]
            if powers:
                return np.median(powers)
            return 6.0  # Default
        
        # Weight by R² and take median
        weighted_powers = np.array(powers) * np.array(rsqs)
        return np.median(weighted_powers) / np.median(rsqs) if np.median(rsqs) > 0 else np.median(powers)
    
    def construct_inter_tissue_adjacency(self,
                                       ts_power: float = 6.0,
                                       ct_power: float = 3.0,
                                       correlation_method: str = 'pearson') -> pd.DataFrame:
        """
        Construct inter-tissue adjacency matrix using X-WGCNA methodology.
        
        Args:
            ts_power: Soft threshold power for tissue-specific connections
            ct_power: Soft threshold power for cross-tissue connections
            correlation_method: Correlation method ('pearson' or 'spearman')
            
        Returns:
            Inter-tissue adjacency matrix
        """
        print(f"Constructing inter-tissue adjacency matrix...")
        print(f"TS power: {ts_power}, CT power: {ct_power}")
        
        total_genes = len(self.all_gene_names)
        adj_matrix = np.zeros((total_genes, total_genes))
        
        tissue_names = list(self.tissue_data.keys())
        
        # 1. Calculate tissue-specific (TS) adjacencies
        print("Calculating tissue-specific adjacencies...")
        for tissue_name, (start_idx, end_idx) in self.tissue_indices.items():
            expr_data = self.tissue_data[tissue_name]['expression']
            
            # Calculate correlation matrix for this tissue
            corr_matrix = expr_data.T.corr(method=correlation_method)
            
            # Apply soft thresholding with TS power
            ts_adj = np.abs(corr_matrix.values) ** ts_power
            np.fill_diagonal(ts_adj, 0)
            
            # Place in the adjacency matrix
            adj_matrix[start_idx:end_idx, start_idx:end_idx] = ts_adj
            
            print(f"  {tissue_name}: {ts_adj.shape} -> positions [{start_idx}:{end_idx}]")
        
        # 2. Calculate cross-tissue (CT) adjacencies
        print("Calculating cross-tissue adjacencies...")
        for i, tissue1 in enumerate(tissue_names):
            for j, tissue2 in enumerate(tissue_names[i+1:], i+1):
                print(f"  Processing {tissue1} - {tissue2}")
                
                # Get tissue indices
                start1, end1 = self.tissue_indices[tissue1]
                start2, end2 = self.tissue_indices[tissue2]
                
                # Find common samples
                samples1 = set(self.tissue_data[tissue1]['samples'])
                samples2 = set(self.tissue_data[tissue2]['samples'])
                common_samples = list(samples1.intersection(samples2))
                
                if len(common_samples) < 5:
                    print(f"    Warning: Only {len(common_samples)} common samples")
                    continue
                
                # Get expression data for common samples
                expr1 = self.tissue_data[tissue1]['expression'].loc[common_samples]
                expr2 = self.tissue_data[tissue2]['expression'].loc[common_samples]
                
                # Calculate cross-correlation
                cross_corr = np.corrcoef(expr1.T, expr2.T)
                cross_corr = cross_corr[0:expr1.shape[1], expr1.shape[1]:]
                
                # Apply soft thresholding with CT power
                ct_adj = np.abs(cross_corr) ** ct_power
                
                # Place in adjacency matrix (symmetric)
                adj_matrix[start1:end1, start2:end2] = ct_adj
                adj_matrix[start2:end2, start1:end1] = ct_adj.T
                
                print(f"    Cross-correlation shape: {cross_corr.shape}")
        
        # Create DataFrame
        self.adjacency_matrix = pd.DataFrame(
            adj_matrix,
            index=self.all_gene_names,
            columns=self.all_gene_names
        )
        
        print(f"Inter-tissue adjacency matrix created: {self.adjacency_matrix.shape}")
        return self.adjacency_matrix
    
    def calculate_topological_overlap(self) -> pd.DataFrame:
        """
        Calculate Topological Overlap Matrix (TOM) for inter-tissue network.
        """
        if self.adjacency_matrix is None:
            raise ValueError("Adjacency matrix not constructed yet")
        
        print("Calculating Topological Overlap Matrix...")
        
        adj_np = self.adjacency_matrix.values
        
        # TOM calculation: TOM_ij = (a_ij + sum_u(a_iu * a_uj)) / (min(k_i, k_j) + 1 - a_ij)
        # where k_i is the connectivity of node i
        
        # Calculate connectivity
        connectivity = adj_np.sum(axis=1)
        
        # Matrix multiplication for numerator
        numerator = adj_np + np.dot(adj_np, adj_np)
        
        # Calculate denominator
        k_matrix = np.outer(connectivity, np.ones_like(connectivity))
        denominator = np.minimum(k_matrix, k_matrix.T) + 1 - adj_np
        
        # Avoid division by zero
        denominator[denominator <= 0] = 1e-10
        
        # Calculate TOM
        tom_matrix = numerator / denominator
        
        # Set diagonal to 0
        np.fill_diagonal(tom_matrix, 0)
        
        self.tom_matrix = pd.DataFrame(
            tom_matrix,
            index=self.adjacency_matrix.index,
            columns=self.adjacency_matrix.columns
        )
        
        print(f"TOM matrix calculated: {self.tom_matrix.shape}")
        return self.tom_matrix
    
    def identify_inter_tissue_modules(self,
                                    min_module_size: int = 30,
                                    cluster_type_threshold: float = 0.95) -> Dict:
        """
        Identify inter-tissue modules using hierarchical clustering.
        
        Args:
            min_module_size: Minimum number of genes in a module
            cluster_type_threshold: Threshold to classify modules as TS or CT
            
        Returns:
            Dictionary with clustering results
        """
        if self.tom_matrix is None:
            raise ValueError("TOM matrix not calculated yet")
        
        print("Identifying inter-tissue modules...")
        
        # Convert TOM to distance matrix
        distance_matrix = 1 - self.tom_matrix.values
        
        # Hierarchical clustering
        linkage_matrix = linkage(squareform(distance_matrix), method='average')
        
        # Dynamic tree cut equivalent (simplified)
        # Use different distance thresholds to find optimal clustering
        best_clustering = None
        best_n_modules = 0
        
        for t in np.arange(0.1, 0.9, 0.05):
            clusters = fcluster(linkage_matrix, t=t, criterion='distance')
            
            # Filter small clusters
            cluster_sizes = pd.Series(clusters).value_counts()
            valid_clusters = cluster_sizes[cluster_sizes >= min_module_size]
            
            if len(valid_clusters) > best_n_modules and len(valid_clusters) < 50:
                best_clustering = clusters
                best_n_modules = len(valid_clusters)
        
        if best_clustering is None:
            best_clustering = fcluster(linkage_matrix, t=0.5, criterion='distance')
        
        # Create module assignments
        gene_modules = {}
        for i, gene in enumerate(self.all_gene_names):
            module_id = best_clustering[i]
            # Filter out small modules
            if pd.Series(best_clustering).value_counts()[module_id] >= min_module_size:
                gene_modules[gene] = module_id
            else:
                gene_modules[gene] = 0  # Grey module (unassigned)
        
        # Analyze module composition
        module_details = self._analyze_module_composition(gene_modules, cluster_type_threshold)
        
        self.clusters = {
            'gene_modules': gene_modules,
            'module_details': module_details,
            'linkage_matrix': linkage_matrix
        }
        
        print(f"Identified {len(module_details)} modules")
        return self.clusters
    
    def _analyze_module_composition(self,
                                  gene_modules: Dict[str, int],
                                  cluster_type_threshold: float) -> pd.DataFrame:
        """
        Analyze the tissue composition of modules.
        """
        module_ids = set(gene_modules.values())
        module_ids.discard(0)  # Remove grey module
        
        tissue_names = list(self.tissue_data.keys())
        
        details = []
        for module_id in module_ids:
            module_genes = [gene for gene, mod in gene_modules.items() if mod == module_id]
            
            # Count genes per tissue
            tissue_counts = {}
            for tissue in tissue_names:
                tissue_counts[tissue] = sum(1 for gene in module_genes if gene.startswith(f"{tissue}_"))
            
            total_genes = len(module_genes)
            dominant_tissue = max(tissue_counts, key=tissue_counts.get)
            dominant_fraction = tissue_counts[dominant_tissue] / total_genes
            
            # Classify module type
            module_type = 'TS' if dominant_fraction >= cluster_type_threshold else 'CT'
            
            details.append({
                'Module_ID': module_id,
                'Module_Size': total_genes,
                'Module_Type': module_type,
                'Dominant_Tissue': dominant_tissue,
                'Dominant_Fraction': dominant_fraction,
                **{f'{tissue}_Count': tissue_counts[tissue] for tissue in tissue_names}
            })
        
        return pd.DataFrame(details)
    
    def validate_network_properties(self) -> Dict:
        """
        Validate inter-tissue network properties including scale-free topology.
        """
        if self.adjacency_matrix is None:
            raise ValueError("Adjacency matrix not constructed yet")
        
        print("Validating network properties...")
        
        # Overall network connectivity
        connectivity = self.adjacency_matrix.sum(axis=1)
        overall_rsq, overall_slope = self._calculate_scale_free_fit(connectivity.values)
        
        # Tissue-specific connectivity analysis
        tissue_connectivity = {}
        for tissue_name, (start_idx, end_idx) in self.tissue_indices.items():
            ts_adj = self.adjacency_matrix.iloc[start_idx:end_idx, start_idx:end_idx]
            ts_connectivity = ts_adj.sum(axis=1)
            ts_rsq, ts_slope = self._calculate_scale_free_fit(ts_connectivity.values)
            
            tissue_connectivity[tissue_name] = {
                'mean_connectivity': ts_connectivity.mean(),
                'median_connectivity': ts_connectivity.median(),
                'max_connectivity': ts_connectivity.max(),
                'rsq': ts_rsq,
                'slope': ts_slope
            }
        
        # Cross-tissue connectivity analysis
        ct_connectivity = {}
        tissue_names = list(self.tissue_data.keys())
        for i, tissue1 in enumerate(tissue_names):
            for j, tissue2 in enumerate(tissue_names[i+1:], i+1):
                start1, end1 = self.tissue_indices[tissue1]
                start2, end2 = self.tissue_indices[tissue2]
                
                ct_adj = self.adjacency_matrix.iloc[start1:end1, start2:end2]
                ct_conn = ct_adj.sum(axis=1)
                
                ct_connectivity[f"{tissue1}-{tissue2}"] = {
                    'mean_connectivity': ct_conn.mean(),
                    'median_connectivity': ct_conn.median(),
                    'max_connectivity': ct_conn.max()
                }
        
        validation_results = {
            'overall_network': {
                'rsq': overall_rsq,
                'slope': overall_slope,
                'mean_connectivity': connectivity.mean(),
                'median_connectivity': connectivity.median()
            },
            'tissue_specific': tissue_connectivity,
            'cross_tissue': ct_connectivity
        }
        
        return validation_results
