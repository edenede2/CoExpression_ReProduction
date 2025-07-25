"""
Network Analysis Module
Implements WGCNA-style network construction and analysis
"""

import pandas as pd
import numpy as np
import networkx as nx
from scipy import stats
from scipy.cluster.hierarchy import linkage, fcluster
from sklearn.metrics import adjusted_rand_score
from sklearn.cluster import KMeans
from typing import Dict, List, Tuple, Optional, Union
import warnings

class NetworkAnalysis:
    """Class for constructing and analyzing coexpression networks."""
    
    def __init__(self):
        """Initialize the network analyzer."""
        self.adjacency_matrix = None
        self.tom_matrix = None
        self.modules = None
        self.power = None
        self.analysis_type = 'traditional'  # 'traditional' or 'inter_tissue'
        self.inter_tissue_analyzer = None
        
    def set_analysis_type(self, analysis_type: str) -> None:
        """
        Set the type of network analysis to perform.
        
        Args:
            analysis_type: 'traditional' for single-condition WGCNA or 
                          'inter_tissue' for X-WGCNA cross-tissue analysis
        """
        if analysis_type not in ['traditional', 'inter_tissue']:
            raise ValueError("analysis_type must be 'traditional' or 'inter_tissue'")
        
        self.analysis_type = analysis_type
        
        if analysis_type == 'inter_tissue':
            from inter_tissue_network_analysis import InterTissueNetworkAnalysis
            self.inter_tissue_analyzer = InterTissueNetworkAnalysis()
            print("Initialized inter-tissue network analysis")
        else:
            print("Using traditional WGCNA analysis")
    
    def load_tissue_data_for_inter_tissue_analysis(self,
                                                  tissue_expression_data: Dict[str, pd.DataFrame],
                                                  tissue_sample_mapping: Dict[str, List[str]]) -> None:
        """
        Load tissue data for inter-tissue analysis.
        
        Args:
            tissue_expression_data: Dict mapping tissue names to expression matrices
            tissue_sample_mapping: Dict mapping tissue names to sample IDs
        """
        if self.analysis_type != 'inter_tissue':
            raise ValueError("Must set analysis_type to 'inter_tissue' first")
        
        self.inter_tissue_analyzer.load_tissue_data(
            tissue_expression_data, tissue_sample_mapping
        )
    
    def determine_inter_tissue_powers(self,
                                    powers_to_test: List[float] = None,
                                    target_rsq: float = 0.8) -> Dict:
        """
        Determine optimal powers for inter-tissue analysis.
        
        Args:
            powers_to_test: List of powers to test
            target_rsq: Target R² for scale-free topology
            
        Returns:
            Power analysis results
        """
        if self.analysis_type != 'inter_tissue':
            raise ValueError("Must set analysis_type to 'inter_tissue' first")
        
        return self.inter_tissue_analyzer.determine_optimal_powers(
            powers_to_test, target_rsq
        )
    
    def construct_inter_tissue_network(self,
                                     ts_power: float = 6.0,
                                     ct_power: float = 3.0,
                                     correlation_method: str = 'pearson',
                                     min_module_size: int = 30) -> Dict:
        """
        Construct inter-tissue network using X-WGCNA methodology.
        
        Args:
            ts_power: Soft threshold power for tissue-specific connections
            ct_power: Soft threshold power for cross-tissue connections
            correlation_method: Correlation method
            min_module_size: Minimum module size
            
        Returns:
            Network construction results
        """
        if self.analysis_type != 'inter_tissue':
            raise ValueError("Must set analysis_type to 'inter_tissue' first")
        
        print("Constructing inter-tissue network...")
        
        # Construct adjacency matrix
        self.adjacency_matrix = self.inter_tissue_analyzer.construct_inter_tissue_adjacency(
            ts_power, ct_power, correlation_method
        )
        
        # Calculate TOM
        self.tom_matrix = self.inter_tissue_analyzer.calculate_topological_overlap()
        
        # Identify modules
        clusters = self.inter_tissue_analyzer.identify_inter_tissue_modules(min_module_size)
        self.modules = clusters['gene_modules']
        
        # Validate network properties
        validation = self.inter_tissue_analyzer.validate_network_properties()
        
        results = {
            'adjacency_matrix': self.adjacency_matrix,
            'tom_matrix': self.tom_matrix,
            'modules': self.modules,
            'module_details': clusters['module_details'],
            'ts_power': ts_power,
            'ct_power': ct_power,
            'validation': validation,
            'analysis_type': 'inter_tissue'
        }
        
        print(f"Inter-tissue network constructed with {len(clusters['module_details'])} modules")
        return results
    
    def calculate_correlation_matrix(self, 
                                     expression_df: pd.DataFrame,
                                     method: str = 'pearson') -> pd.DataFrame:
        """
        Calculate gene-gene correlation matrix.
        
        Args:
            expression_df: Expression data (genes x samples)
            method: Correlation method ('pearson', 'spearman', 'kendall')
            
        Returns:
            Correlation matrix DataFrame
        """
        print(f"Calculating {method} correlation matrix...")
        
        if method == 'pearson':
            corr_matrix = expression_df.T.corr()
        elif method == 'spearman':
            corr_matrix = expression_df.T.corr(method='spearman')
        elif method == 'kendall':
            corr_matrix = expression_df.T.corr(method='kendall')
        else:
            raise ValueError(f"Unknown correlation method: {method}")
        
        print(f"Correlation matrix shape: {corr_matrix.shape}")
        return corr_matrix
    
    def determine_soft_threshold_power(self, 
                                       corr_matrix: pd.DataFrame,
                                       powers: List[int] = None,
                                       target_rsq: float = 0.85) -> Tuple[int, pd.DataFrame]:
        """
        Determine optimal soft threshold power for scale-free topology.
        Implementation based on WGCNA's pickSoftThreshold function.
        
        Args:
            corr_matrix: Gene correlation matrix
            powers: List of powers to test
            target_rsq: Target R-squared for scale-free topology
            
        Returns:
            tuple: (optimal_power, power_analysis_results)
        """
        if powers is None:
            powers = list(range(1, 21))
        
        print(f"Determining optimal soft threshold power (target R² = {target_rsq})...")
        
        results = []
        
        for power in powers:
            print(f"Testing power {power}...")
            
            # Calculate adjacency matrix
            adj_matrix = self._calculate_adjacency_matrix(corr_matrix, power)
            
            # Calculate connectivity (sum of adjacencies for each gene)
            connectivity = adj_matrix.sum(axis=1).values
            
            # WGCNA-style scale-free analysis
            # Remove genes with zero connectivity
            nonzero_k = connectivity[connectivity > 0]
            
            if len(nonzero_k) < 10:
                rsq = 0.0
                slope = 0.0
                mean_connectivity = 0.0
            else:
                # Bin the connectivity values and calculate frequency
                # Use log-spaced bins for better scale-free analysis
                min_k = max(1, np.min(nonzero_k))
                max_k = np.max(nonzero_k)
                
                if max_k > min_k:
                    # Create bins in log space
                    n_bins = min(50, int(np.sqrt(len(nonzero_k))))
                    bins = np.logspace(np.log10(min_k), np.log10(max_k), n_bins)
                    
                    # Calculate histogram
                    hist, bin_edges = np.histogram(nonzero_k, bins=bins)
                    
                    # Use bin centers
                    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
                    
                    # Remove empty bins
                    valid_bins = hist > 0
                    if valid_bins.sum() >= 5:  # Need at least 5 points
                        k_values = bin_centers[valid_bins]
                        p_k_values = hist[valid_bins] / len(nonzero_k)
                        
                        # Convert to log scale
                        log_k = np.log10(k_values)
                        log_p_k = np.log10(p_k_values)
                        
                        # Linear regression
                        slope, intercept, r_value, p_value, std_err = stats.linregress(log_k, log_p_k)
                        rsq = r_value ** 2
                        
                        # For scale-free networks, we want negative slope
                        if slope > 0:
                            rsq = 0.0
                    else:
                        rsq = 0.0
                        slope = 0.0
                else:
                    rsq = 0.0
                    slope = 0.0
            
            # Calculate mean connectivity
            mean_connectivity = np.mean(connectivity)
            
            results.append({
                'power': power,
                'rsq': rsq,
                'slope': slope,
                'mean_connectivity': mean_connectivity,
                'median_connectivity': np.median(connectivity),
                'max_connectivity': np.max(connectivity)
            })
        
        results_df = pd.DataFrame(results)
        
        # Find optimal power
        # 1. First try to find power with R² >= target
        high_rsq = results_df[results_df['rsq'] >= target_rsq]
        
        if len(high_rsq) > 0:
            # Choose the lowest power that meets the criterion
            optimal_power = high_rsq.iloc[0]['power']
        else:
            # If no power meets target, use more relaxed criteria
            # Look for R² >= 0.7 or the highest R²
            relaxed_rsq = results_df[results_df['rsq'] >= 0.7]
            
            if len(relaxed_rsq) > 0:
                optimal_power = relaxed_rsq.iloc[0]['power']
            else:
                # Choose power with highest R², but avoid very high powers
                # that might lead to disconnected networks
                reasonable_powers = results_df[results_df['mean_connectivity'] >= 1]
                if len(reasonable_powers) > 0:
                    optimal_power = reasonable_powers.loc[reasonable_powers['rsq'].idxmax(), 'power']
                else:
                    optimal_power = results_df.loc[results_df['rsq'].idxmax(), 'power']
            
            optimal_rsq = results_df[results_df['power'] == optimal_power]['rsq'].iloc[0]
            if optimal_rsq < target_rsq:
                print(f"Warning: No power achieved target R² of {target_rsq}")
        
        optimal_rsq = results_df[results_df['power'] == optimal_power]['rsq'].iloc[0]
        print(f"Optimal power: {optimal_power} (R² = {optimal_rsq:.3f})")
        
        return optimal_power, results_df
    
    def _calculate_adjacency_matrix(self, 
                                    corr_matrix: pd.DataFrame,
                                    power: int) -> pd.DataFrame:
        """
        Calculate adjacency matrix using soft thresholding.
        
        Args:
            corr_matrix: Gene correlation matrix
            power: Soft threshold power
            
        Returns:
            Adjacency matrix
        """
        # Take absolute values and apply power
        adj_matrix = np.abs(corr_matrix) ** power
        
        # Set diagonal to 0 (no self-connections)
        np.fill_diagonal(adj_matrix.values, 0)
        
        return adj_matrix
    
    def calculate_topological_overlap(self, adj_matrix: pd.DataFrame) -> pd.DataFrame:
        """
        Calculate Topological Overlap Matrix (TOM).
        
        Args:
            adj_matrix: Adjacency matrix
            
        Returns:
            TOM matrix
        """
        print("Calculating Topological Overlap Matrix...")
        
        # Convert to numpy for faster computation
        adj_np = adj_matrix.values
        n_genes = adj_np.shape[0]
        
        # Calculate TOM
        tom_np = np.zeros_like(adj_np)
        
        for i in range(n_genes):
            if i % 1000 == 0:
                print(f"Processing gene {i}/{n_genes}")
            
            for j in range(i+1, n_genes):
                # Numerator: l_ij + a_ij
                l_ij = np.sum(adj_np[i, :] * adj_np[j, :])
                numerator = l_ij + adj_np[i, j]
                
                # Denominator: min(k_i, k_j) + 1 - a_ij
                k_i = np.sum(adj_np[i, :])
                k_j = np.sum(adj_np[j, :])
                denominator = min(k_i, k_j) + 1 - adj_np[i, j]
                
                if denominator > 0:
                    tom_value = numerator / denominator
                else:
                    tom_value = 0
                
                tom_np[i, j] = tom_value
                tom_np[j, i] = tom_value
        
        # Set diagonal to 1
        np.fill_diagonal(tom_np, 1)
        
        tom_matrix = pd.DataFrame(tom_np, index=adj_matrix.index, columns=adj_matrix.columns)
        
        print("TOM calculation completed")
        return tom_matrix
    
    def detect_modules(self, 
                       tom_matrix: pd.DataFrame,
                       min_module_size: int = 30,
                       method: str = 'hierarchical') -> Dict[str, int]:
        """
        Detect gene modules using clustering.
        
        Args:
            tom_matrix: Topological overlap matrix
            min_module_size: Minimum size for modules
            method: Clustering method ('hierarchical' or 'kmeans')
            
        Returns:
            Dictionary mapping gene IDs to module numbers
        """
        print(f"Detecting modules using {method} clustering...")
        
        if method == 'hierarchical':
            modules = self._hierarchical_clustering(tom_matrix, min_module_size)
        elif method == 'kmeans':
            modules = self._kmeans_clustering(tom_matrix, min_module_size)
        else:
            raise ValueError(f"Unknown clustering method: {method}")
        
        module_counts = pd.Series(modules).value_counts().sort_index()
        print(f"Detected {len(module_counts)} modules")
        print(f"Module sizes: {dict(module_counts)}")
        
        return modules
    
    def _hierarchical_clustering(self, 
                                 tom_matrix: pd.DataFrame,
                                 min_module_size: int) -> Dict[str, int]:
        """Perform hierarchical clustering to detect modules."""
        # Calculate dissimilarity (1 - TOM)
        dissim_matrix = 1 - tom_matrix.values
        
        # Perform hierarchical clustering
        # Convert to condensed distance matrix for linkage
        n_genes = dissim_matrix.shape[0]
        condensed_dist = []
        
        for i in range(n_genes):
            for j in range(i+1, n_genes):
                condensed_dist.append(dissim_matrix[i, j])
        
        condensed_dist = np.array(condensed_dist)
        
        # Perform linkage
        linkage_matrix = linkage(condensed_dist, method='average')
        
        # Dynamic tree cutting - try different cut heights
        best_modules = None
        best_n_modules = 0
        
        for cut_height in np.arange(0.1, 0.9, 0.1):
            clusters = fcluster(linkage_matrix, cut_height, criterion='distance')
            
            # Filter small modules
            cluster_counts = pd.Series(clusters).value_counts()
            valid_clusters = cluster_counts[cluster_counts >= min_module_size].index
            
            if len(valid_clusters) > best_n_modules:
                best_n_modules = len(valid_clusters)
                best_modules = clusters.copy()
                
                # Reassign small modules to module 0 (grey module)
                for i, cluster in enumerate(clusters):
                    if cluster not in valid_clusters:
                        best_modules[i] = 0
        
        # If no valid clusters found, assign all genes to module 0
        if best_modules is None:
            best_modules = np.zeros(n_genes, dtype=int)
            print(f"Warning: No modules found with minimum size {min_module_size}. All genes assigned to module 0.")
        
        # Create gene-to-module mapping
        modules = dict(zip(tom_matrix.index, best_modules))
        
        return modules
    
    def _kmeans_clustering(self, 
                           tom_matrix: pd.DataFrame,
                           min_module_size: int) -> Dict[str, int]:
        """Perform K-means clustering to detect modules."""
        # Estimate number of clusters
        n_genes = tom_matrix.shape[0]
        max_k = min(20, n_genes // min_module_size)
        
        best_k = 5  # Default
        best_score = -1
        
        for k in range(2, max_k + 1):
            kmeans = KMeans(n_clusters=k, random_state=42, n_init=10)
            clusters = kmeans.fit_predict(tom_matrix.values)
            
            # Check if all clusters meet minimum size
            cluster_counts = pd.Series(clusters).value_counts()
            if all(count >= min_module_size for count in cluster_counts):
                # Use silhouette score or other metric
                score = k  # Simple heuristic - prefer more clusters if size requirements met
                if score > best_score:
                    best_score = score
                    best_k = k
        
        # Final clustering with best k
        kmeans = KMeans(n_clusters=best_k, random_state=42, n_init=10)
        clusters = kmeans.fit_predict(tom_matrix.values)
        
        # Create gene-to-module mapping
        modules = dict(zip(tom_matrix.index, clusters + 1))  # +1 to start from module 1
        
        return modules
    
    def construct_network(self, 
                          expression_df: pd.DataFrame,
                          power: Optional[int] = None,
                          min_module_size: int = 30,
                          correlation_method: str = 'pearson') -> Dict:
        """
        Complete network construction pipeline.
        
        Args:
            expression_df: Expression data (genes x samples)
            power: Soft threshold power (auto-determined if None)
            min_module_size: Minimum module size
            correlation_method: Correlation method
            
        Returns:
            Dictionary with network analysis results
        """
        print("Starting network construction pipeline...")
        
        # Calculate correlation matrix
        corr_matrix = self.calculate_correlation_matrix(expression_df, correlation_method)
        
        # Determine optimal power if not provided
        if power is None:
            power, power_analysis = self.determine_soft_threshold_power(corr_matrix)
        else:
            power_analysis = None
        
        self.power = power
        
        # Calculate adjacency matrix
        print(f"Calculating adjacency matrix with power {power}...")
        self.adjacency_matrix = self._calculate_adjacency_matrix(corr_matrix, power)
        
        # Calculate TOM
        self.tom_matrix = self.calculate_topological_overlap(self.adjacency_matrix)
        
        # Detect modules
        self.modules = self.detect_modules(self.tom_matrix, min_module_size)
        
        # Create module assignment DataFrame
        module_df = pd.DataFrame({
            'gene_id': list(self.modules.keys()),
            'module': list(self.modules.values())
        })
        
        results = {
            'correlation_matrix': corr_matrix,
            'adjacency_matrix': self.adjacency_matrix,
            'tom_matrix': self.tom_matrix,
            'modules': self.modules,
            'module_df': module_df,
            'power': power,
            'power_analysis': power_analysis
        }
        
        print("Network construction completed!")
        return results
    
    def calculate_module_eigengenes(self, 
                                    expression_df: pd.DataFrame,
                                    modules: Dict[str, int]) -> pd.DataFrame:
        """
        Calculate module eigengenes (first principal component of each module).
        
        Args:
            expression_df: Expression data (genes x samples)
            modules: Gene-to-module mapping
            
        Returns:
            Module eigengenes DataFrame (modules x samples)
        """
        print("Calculating module eigengenes...")
        
        from sklearn.decomposition import PCA
        
        module_eigengenes = {}
        
        for module_id in set(modules.values()):
            if module_id == 0:  # Skip grey module
                continue
                
            # Get genes in this module
            module_genes = [gene for gene, mod in modules.items() if mod == module_id]
            
            if len(module_genes) < 3:  # Need at least 3 genes for PCA
                continue
                
            # Get expression data for module genes
            module_expr = expression_df.loc[module_genes].T  # Samples x genes
            
            # Calculate first principal component
            pca = PCA(n_components=1)
            eigengene = pca.fit_transform(module_expr).flatten()
            
            module_eigengenes[f'ME{module_id}'] = eigengene
        
        eigengenes_df = pd.DataFrame(module_eigengenes, index=expression_df.columns)
        
        print(f"Calculated eigengenes for {len(module_eigengenes)} modules")
        return eigengenes_df
