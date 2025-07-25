"""
Data Preprocessing Module
Handles filtering, normalization, and quality control of expression data
"""

import pandas as pd
import numpy as np
from sklearn.decomposition import TruncatedSVD
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LinearRegression
from scipy import stats
from scipy.stats import chi2
from typing import List, Tuple, Dict, Optional

class DataPreprocessor:
    """Class for preprocessing expression data."""
    
    def __init__(self):
        """Initialize the preprocessor."""
        self.scaler = StandardScaler()
        self.pca = None
        
    def filter_protein_coding_genes(self, 
                                    expression_df: pd.DataFrame,
                                    hgnc_data: pd.DataFrame) -> pd.DataFrame:
        """
        Filter expression data to keep only protein coding genes.
        
        Args:
            expression_df: Expression data DataFrame
            hgnc_data: HGNC protein coding genes DataFrame
            
        Returns:
            Filtered expression DataFrame
        """
        # Extract gene symbols from Ensembl IDs
        expression_df['gene_symbol'] = expression_df.index.str.split('.').str[0]
        
        # Get protein coding gene symbols
        protein_coding_symbols = set(hgnc_data['ensembl_gene_id'].dropna())
        
        # Filter expression data
        protein_coding_mask = expression_df['gene_symbol'].isin(protein_coding_symbols)
        filtered_df = expression_df[protein_coding_mask].copy()
        
        # Remove the temporary gene_symbol column
        filtered_df = filtered_df.drop('gene_symbol', axis=1)
        
        print(f"Filtered to {len(filtered_df)} protein coding genes (from {len(expression_df)})")
        return filtered_df
    
    def log_transform(self, expression_df: pd.DataFrame, pseudocount: float = 1.0) -> pd.DataFrame:
        """
        Apply log2 transformation to expression data.
        
        Args:
            expression_df: Expression data DataFrame
            pseudocount: Pseudocount to add before log transformation
            
        Returns:
            Log-transformed expression DataFrame
        """
        print("Applying log2(TPM + 1) transformation...")
        log_df = np.log2(expression_df + pseudocount)
        return log_df
    
    def filter_low_variance_genes(self, 
                                  expression_df: pd.DataFrame,
                                  min_variance: float = 0.02) -> pd.DataFrame:
        """
        Filter out genes with low variance across samples.
        
        Args:
            expression_df: Expression data DataFrame
            min_variance: Minimum variance threshold
            
        Returns:
            Variance-filtered expression DataFrame
        """
        gene_variance = expression_df.var(axis=1)
        high_var_genes = gene_variance[gene_variance >= min_variance].index
        
        filtered_df = expression_df.loc[high_var_genes]
        
        print(f"Filtered to {len(filtered_df)} genes with variance >= {min_variance}")
        return filtered_df
    
    def filter_low_expression_genes(self, 
                                    expression_df: pd.DataFrame,
                                    min_expression: float = 0.1,
                                    min_samples_fraction: float = 0.2) -> pd.DataFrame:
        """
        Filter genes with low expression across samples.
        
        Args:
            expression_df: Expression data DataFrame
            min_expression: Minimum expression threshold (log2(TPM+1) > min_expression)
            min_samples_fraction: Minimum fraction of samples that must meet threshold
            
        Returns:
            Expression-filtered DataFrame
        """
        # In the original script:
        # tissue_matrix[(tissue_matrix.T < np.log2(0.1+1)).sum() > 0.2*tissue_matrix.shape[1]]
        # This filters out genes with expression < log2(0.1+1) in more than 20% of samples
        # So we keep genes with low expression in less than 20% of samples
        
        # Convert fraction to number of samples
        min_samples = int(min_samples_fraction * expression_df.shape[1])
        
        # Count samples with expression below threshold for each gene
        low_expressed_samples = (expression_df < min_expression).sum(axis=1)
        
        # Keep genes with low expression in fewer than min_samples
        keep_genes = low_expressed_samples <= min_samples
        filtered_df = expression_df[keep_genes]
        
        print(f"Filtered to {len(filtered_df)} genes (removed {len(expression_df) - len(filtered_df)} genes with low expression in > {min_samples_fraction*100}% of samples)")
        return filtered_df
    
    def detect_outliers_pca(self, 
                           expression_df: pd.DataFrame,
                           n_components: int = 20,
                           contamination: float = 0.01) -> Tuple[List[str], pd.DataFrame]:
        """
        Detect outlier samples using TruncatedSVD and Mahalanobis distance.
        
        Args:
            expression_df: Expression data DataFrame (genes x samples)
            n_components: Number of SVD components
            contamination: Expected fraction of outliers (defaults to 0.01 which corresponds to 99% confidence level)
            
        Returns:
            tuple: (outlier_sample_ids, pca_results_df)
        """
        print(f"Detecting outliers using TruncatedSVD with {n_components} components...")
        
        # Ensure there are no missing values in the data
        matrix_val = expression_df.copy()
        if matrix_val.isna().sum().sum() > 0:
            print("Filling NA values with small value (0.0000001)")
            matrix_val = matrix_val.fillna(0.0000001)
        
        # Save sample IDs for reference
        samples = matrix_val.columns.tolist()
        
        # Apply TruncatedSVD as in the original script
        svd = TruncatedSVD(n_components=n_components, random_state=1001)
        svd.fit(matrix_val)
        components = svd.components_.T  # Transpose to match the original script
        self.pca = svd  # Store for later use
        
        # Calculate covariance matrix and its inverse
        covariance = np.cov(components, rowvar=False)
        inv_covariance = np.linalg.matrix_power(covariance, -1)
        
        # Calculate center point
        center_point = np.mean(components, axis=0)
        
        # Calculate Mahalanobis distances
        distances = []
        for i, val in enumerate(components):
            distance = (val - center_point).T.dot(inv_covariance).dot(val - center_point)
            distances.append(distance)
        distances = np.array(distances)
        
        # Determine threshold using chi-square distribution with 99% confidence level
        # In original code: cutoff = chi2.ppf(0.99, components.shape[1]*2)
        threshold = chi2.ppf(1 - contamination, components.shape[1]*2)
        outliers = distances > threshold
        
        # Get outlier sample IDs
        outlier_indices = np.where(outliers)[0]
        outlier_samples = [samples[i] for i in outlier_indices]
        
        # Create results DataFrame
        pca_df = pd.DataFrame(
            components,
            index=samples,
            columns=[f'PC{i+1}' for i in range(n_components)]
        )
        pca_df['mahalanobis_distance'] = distances
        pca_df['is_outlier'] = outliers
        
        print(f"Detected {len(outlier_samples)} outlier samples")
        if hasattr(svd, 'explained_variance_ratio_'):
            print(f"SVD explained variance ratio: {svd.explained_variance_ratio_[:5]}")
        
        return outlier_samples, pca_df
    
    def remove_outlier_samples(self, 
                              expression_df: pd.DataFrame,
                              outlier_samples: List[str]) -> pd.DataFrame:
        """
        Remove outlier samples from expression data.
        
        Args:
            expression_df: Expression data DataFrame
            outlier_samples: List of outlier sample IDs
            
        Returns:
            Expression DataFrame with outliers removed
        """
        clean_samples = [s for s in expression_df.columns if s not in outlier_samples]
        clean_df = expression_df[clean_samples].copy()
        
        print(f"Removed {len(outlier_samples)} outlier samples")
        print(f"Remaining samples: {len(clean_df.columns)}")
        
        return clean_df
    
    def quantile_normalize(self, expression_df: pd.DataFrame) -> pd.DataFrame:
        """
        Apply quantile normalization to expression data.
        
        Args:
            expression_df: Expression data DataFrame
            
        Returns:
            Quantile-normalized expression DataFrame
        """
        print("Applying quantile normalization...")
        
        # Following the old script's implementation
        df = expression_df
        
        # Sort each column (sample)
        df_sorted = pd.DataFrame(
            np.sort(df.values, axis=0),
            index=df.index,
            columns=df.columns
        )
        
        # Calculate row means (quantile means)
        df_mean = df_sorted.mean(axis=1)
        
        # Assign ranks 1 to n
        df_mean.index = np.arange(1, len(df_mean) + 1)
        
        # Use ranks to map to quantile means
        # This is a more efficient implementation of the same logic in the old script
        df_qn = df.rank(method="min").stack().astype(int).map(df_mean).unstack()
        
        print("Quantile normalization completed")
        return df_qn
    
    def regress_confounders(self, expression_df: pd.DataFrame, 
                           confounders_df: pd.DataFrame) -> pd.DataFrame:
        """
        Regress out confounding factors from expression data.
        
        Args:
            expression_df: Expression data DataFrame (genes x samples)
            confounders_df: Confounding factors DataFrame (samples x factors)
            
        Returns:
            Residualized expression DataFrame
        """
        print("Regressing out confounding factors...")
        
        # Transpose for regression (samples x genes) as in the original script
        expr_t = expression_df.T  # samples x genes
        
        # Create DataFrame to hold residuals with same shape as expr_t
        sk_resid_mat = expr_t.copy(deep=True)
        for col in sk_resid_mat.columns:
            sk_resid_mat[col].values[:] = 0
        
        # Make sure we have common samples between expression and confounders
        common_samples = sorted(set(expr_t.index) & set(confounders_df.index))
        if len(common_samples) == 0:
            raise ValueError("No common samples between expression data and confounders")
        
        # Filter data to common samples
        expr_t_common = expr_t.loc[common_samples]
        confounders_common = confounders_df.loc[common_samples]
        
        # Prepare confounders for regression - convert all to numeric
        X = pd.get_dummies(confounders_common, drop_first=True)
        
        # Check for non-numeric columns and remove them
        numeric_cols = X.select_dtypes(include=['number']).columns
        if len(numeric_cols) < X.shape[1]:
            print(f"Warning: Removed {X.shape[1] - len(numeric_cols)} non-numeric columns from confounders")
            X = X[numeric_cols]
        
        print(f"Using {X.shape[1]} confounders for regression with {len(common_samples)} samples")
        
        # Fit linear model for each gene
        for i in range(expr_t_common.shape[1]):  # loop through all genes
            gene = expr_t_common.columns[i]
            y = expr_t_common[gene].values  # all samples for one gene
            
            # Fit linear regression model and get predictions
            try:
                reg = LinearRegression()
                reg.fit(X, y)
                y_pred = reg.predict(X)
                
                # Store residuals for common samples
                sk_resid_mat.loc[common_samples, gene] = y - y_pred
            except Exception as e:
                print(f"Warning: Failed to regress gene {gene}: {str(e)}")
                # Keep original values for this gene
                sk_resid_mat.loc[common_samples, gene] = expr_t_common[gene].values
        
        # Transpose back to genes x samples and apply quantile normalization
        # as in the old script
        sk_resid_mat_norm = self.quantile_normalize(sk_resid_mat.T)
        
        print("Confounding factors regression completed")
        return sk_resid_mat_norm
        

    
    def preprocess_expression_data(self,
                                   expression_df: pd.DataFrame,
                                   hgnc_data: Optional[pd.DataFrame] = None,
                                   confounders_df: Optional[pd.DataFrame] = None,
                                   apply_log: bool = True,
                                   filter_low_expr: bool = True,
                                   filter_low_var: bool = True,
                                   detect_outliers: bool = True,
                                   quantile_norm: bool = True, 
                                   regress_confounders: bool = False,
                                   min_expression: float = 0.1,
                                   min_samples_fraction: float = 0.2,  
                                   min_variance: float = 0.02,
                                   outlier_contamination: float = 0.01) -> Tuple[pd.DataFrame, Dict]:  
        """
        Complete preprocessing pipeline for expression data.
        
        Args:
            expression_df: Raw expression data DataFrame
            hgnc_data: HGNC protein coding genes (optional)
            confounders_df: Confounding factors DataFrame (optional)
            apply_log: Whether to apply log transformation
            filter_low_expr: Whether to filter low expression genes
            filter_low_var: Whether to filter low variance genes
            detect_outliers: Whether to detect and remove outliers
            quantile_norm: Whether to apply quantile normalization (default: True)
            regress_confounders: Whether to regress out confounding factors
            min_expression: Minimum expression threshold
            min_samples_fraction: Minimum fraction of samples where genes should be expressed
            min_variance: Minimum variance threshold
            outlier_contamination: Expected outlier fraction (default: 0.01 for 99% confidence)
            
        Returns:
            tuple: (processed_expression_df, processing_info)
        """
        print("Starting expression data preprocessing pipeline...")
        
        processed_df = expression_df.copy()
        processing_info = {
            'original_shape': expression_df.shape,
            'steps_applied': []
        }
        
        # Filter protein coding genes
        if hgnc_data is not None:
            processed_df = self.filter_protein_coding_genes(processed_df, hgnc_data)
            processing_info['steps_applied'].append('protein_coding_filter')
            processing_info['after_protein_coding'] = processed_df.shape
        
        # Log transformation
        if apply_log:
            processed_df = self.log_transform(processed_df)
            processing_info['steps_applied'].append('log_transform')
        
        # Filter low expression genes
        if filter_low_expr:
            processed_df = self.filter_low_expression_genes(
                processed_df, min_expression, min_samples_fraction
            )
            processing_info['steps_applied'].append('low_expression_filter')
            processing_info['after_expression_filter'] = processed_df.shape
        
        # Filter low variance genes
        if filter_low_var:
            processed_df = self.filter_low_variance_genes(processed_df, min_variance)
            processing_info['steps_applied'].append('low_variance_filter')
            processing_info['after_variance_filter'] = processed_df.shape
        
        # Detect and remove outliers
        if detect_outliers:
            outlier_samples, pca_results = self.detect_outliers_pca(
                processed_df, contamination=outlier_contamination
            )
            processed_df = self.remove_outlier_samples(processed_df, outlier_samples)
            processing_info['steps_applied'].append('outlier_removal')
            processing_info['outlier_samples'] = outlier_samples
            processing_info['pca_results'] = pca_results
            processing_info['after_outlier_removal'] = processed_df.shape
        
        # Quantile normalization
        if quantile_norm:
            processed_df = self.quantile_normalize(processed_df)
            processing_info['steps_applied'].append('quantile_normalization')
        
        # Regress out confounding factors
        if regress_confounders and confounders_df is not None:
            processed_df = self.regress_confounders(processed_df, confounders_df)
            processing_info['steps_applied'].append('confounder_regression')
            processing_info['after_confounder_regression'] = processed_df.shape
        
        processing_info['final_shape'] = processed_df.shape
        
        print(f"Preprocessing completed. Final shape: {processed_df.shape}")
        print(f"Steps applied: {processing_info['steps_applied']}")
        
        return processed_df, processing_info
