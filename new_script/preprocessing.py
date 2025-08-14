"""
Data Preprocessing Module
Handles filtering, normalization, and quality control of expression data
"""

import pandas as pd
import numpy as np
from sklearn.decomposition import TruncatedSVD
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LinearRegression
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
    
    def generate_confounders_from_metadata(self, 
                                           samples: List[str],
                                           sample_attrs: pd.DataFrame,
                                           subject_phenos: pd.DataFrame) -> pd.DataFrame:
        """
        Generate confounders DataFrame from GTEx metadata, following the original script methodology exactly.
        
        Based on the original script process_tissue function:
        1. Merge sample attributes with subject phenotypes 
        2. Process SMGEBTCH (batch) by grouping rare batches into singletons
        3. Convert AGE from ranges to numeric 
        4. Apply get_dummies to SMGEBTCH 
        5. Drop ['SMTSD','SMGEBTCH','SUBJID'] 
        6. Use ALL remaining columns as confounders EXCEPT ['AGE', 'SMRIN', 'DTHHRDY']
        
        Args:
            samples: List of sample IDs
            sample_attrs: Sample attributes DataFrame  
            subject_phenos: Subject phenotypes DataFrame
            
        Returns:
            Confounders DataFrame with samples as index
        """
        print("Generating confounders DataFrame from GTEx metadata following original script...")
        
        # Filter to available samples
        filtered_attrs = sample_attrs.loc[list(set(samples) & set(sample_attrs.index))].copy()
        
        if len(filtered_attrs) == 0:
            raise ValueError("No samples found in sample attributes")
        
        # Extract subject IDs from sample IDs (format: GTEX-XXXXX-XXXX-...)
        filtered_attrs['SUBJID'] = filtered_attrs.index.str.extract(r'(GTEX-[^-]+)', expand=False)
        
        # Merge with phenotypes to get AGE and DTHHRDY (following original script)
        merged_data = filtered_attrs.merge(
            subject_phenos[['AGE', 'DTHHRDY']],
            left_on='SUBJID',
            right_index=True,
            how='left'
        )
        
        # Process batch information (SMGEBTCH) exactly as in original script
        if 'SMGEBTCH' in merged_data.columns:
            # Count batch occurrences - group rare batches (count <= 1) into singletons
            batch_counts = merged_data['SMGEBTCH'].value_counts(normalize=False, sort=True)
            batch_names = batch_counts.index.tolist()
            batch_values = batch_counts.tolist()
            valid_batches = []
            
            for i in range(len(batch_names)):
                if batch_values[i] > 1:
                    valid_batches.append(batch_names[i])
            
            # Replace rare batches with singleton category
            batches = []
            for i in range(merged_data.shape[0]):
                if merged_data['SMGEBTCH'].iloc[i] in valid_batches:
                    batches.append(merged_data['SMGEBTCH'].iloc[i])
                else:
                    batches.append('ASINGLETON_SMGEBTCH')
            
            merged_data['SMGEBTCH'] = batches
            print(f"Processed {len(valid_batches)} valid batches, {len(batch_names) - len(valid_batches)} rare batches grouped as singletons")
        
        # Handle AGE conversion from ranges to numeric (e.g., "20-29" -> 20)
        # Following original script: sample3['AGE'] = [elem[:2] for elem in agegroup]
        if 'AGE' in merged_data.columns:
            try:
                # Original script extracts first 2 characters, not split by '-'
                agegroup = merged_data['AGE'].astype(str).tolist()
                merged_data['AGE'] = [elem[:2] for elem in agegroup]
                merged_data['AGE'] = pd.to_numeric(merged_data['AGE'], errors='coerce')
                print("Converted AGE from ranges to numeric values (first 2 characters)")
            except (ValueError, TypeError) as e:
                print(f"Warning: Could not convert AGE to numeric: {e}")
        
        # Apply get_dummies to SMGEBTCH (as in original script)
        if 'SMGEBTCH' in merged_data.columns:
            smgebtch_dummies = pd.get_dummies(merged_data['SMGEBTCH'], drop_first=True)
        else:
            smgebtch_dummies = pd.DataFrame(index=merged_data.index)
        
        # Drop columns as in original script: ['SMTSD','SMGEBTCH','SUBJID']
        columns_to_drop = ['SMTSD', 'SMGEBTCH', 'SUBJID']
        sample_data = merged_data.drop(columns=[col for col in columns_to_drop if col in merged_data.columns], axis=1)
        
        # Concatenate sample data with SMGEBTCH dummies (as in original script)
        if not smgebtch_dummies.empty:
            final_confounders = pd.concat([sample_data, smgebtch_dummies], join='inner', axis=1)
        else:
            final_confounders = sample_data
        
        # Now drop ['AGE', 'SMRIN', 'DTHHRDY'] as in original script
        # Original: x = sample.drop(['AGE', 'SMRIN', 'DTHHRDY'], axis=1)
        confounders_to_exclude = ['AGE', 'SMRIN', 'DTHHRDY']
        final_confounders = final_confounders.drop(
            columns=[col for col in confounders_to_exclude if col in final_confounders.columns], 
            axis=1
        )
        
        # Ensure all columns are numeric after dummy encoding
        for col in final_confounders.columns:
            if final_confounders[col].dtype == 'object':
                try:
                    final_confounders[col] = pd.to_numeric(final_confounders[col], errors='coerce')
                    if final_confounders[col].isna().any():
                        print(f"Warning: Column {col} contains non-numeric values, dropping it")
                        final_confounders = final_confounders.drop(col, axis=1)
                except (ValueError, TypeError):
                    print(f"Warning: Could not convert column {col} to numeric, dropping it")
                    final_confounders = final_confounders.drop(col, axis=1)
        
        print(f"DEBUG: Final confounders shape before cleaning: {final_confounders.shape}")
        print(f"DEBUG: Final confounders columns: {list(final_confounders.columns)}")
        print(f"DEBUG: Final confounders dtypes: {final_confounders.dtypes.to_dict()}")
        
        # Check each column for problems
        problematic_cols = []
        for col in final_confounders.columns:
            nan_count = final_confounders[col].isna().sum()
            inf_count = np.isinf(final_confounders[col]).sum() if pd.api.types.is_numeric_dtype(final_confounders[col]) else 0
            print(f"DEBUG: Column {col}: NaN={nan_count}, Inf={inf_count}, dtype={final_confounders[col].dtype}")
            
            if nan_count > 0 or inf_count > 0:
                problematic_cols.append(col)
        
        print(f"DEBUG: Problematic columns: {problematic_cols}")
        
        # For now, let's try a different approach - fill NaN with median for numeric columns
        for col in final_confounders.columns:
            if final_confounders[col].dtype in ['object', 'category']:
                print(f"DEBUG: Dropping non-numeric column {col}")
                final_confounders = final_confounders.drop(col, axis=1)
            elif final_confounders[col].isna().any():
                if pd.api.types.is_numeric_dtype(final_confounders[col]):
                    median_val = final_confounders[col].median()
                    print(f"DEBUG: Filling {final_confounders[col].isna().sum()} NaN values in {col} with median {median_val}")
                    final_confounders[col] = final_confounders[col].fillna(median_val)
                else:
                    print(f"DEBUG: Dropping non-numeric column with NaN: {col}")
                    final_confounders = final_confounders.drop(col, axis=1)
        
        print(f"DEBUG: After NaN handling - shape: {final_confounders.shape}")
        
        # Check for remaining infinite values
        for col in final_confounders.columns:
            if np.isinf(final_confounders[col]).any():
                print(f"DEBUG: Replacing infinite values in {col}")
                final_confounders[col] = final_confounders[col].replace([np.inf, -np.inf], final_confounders[col].median())
        
        # Final check - ensure we have data left
        if final_confounders.empty:
            raise ValueError("No valid confounders remaining after cleaning")
        if final_confounders.shape[0] == 0:
            raise ValueError("No samples remaining after confounder cleaning")
        
        print(f"Generated confounders DataFrame with {final_confounders.shape[1]} variables")
        print("Following original script methodology - using ALL metadata EXCEPT AGE, SMRIN, DTHHRDY")
        print(f"Confounder columns: {list(final_confounders.columns)}")
        
        return final_confounders
    
    # def regress_confounders(self, expression_df: pd.DataFrame, 
    #                        confounders_df: Optional[pd.DataFrame] = None,
    #                        sample_attrs: Optional[pd.DataFrame] = None,
    #                        subject_phenos: Optional[pd.DataFrame] = None) -> pd.DataFrame:
    #     """
    #     Regress out confounding factors from expression data.
        
    #     Args:
    #         expression_df: Expression data DataFrame (genes x samples)
    #         confounders_df: Confounding factors DataFrame (samples x factors).
    #                        If None, will generate from sample_attrs and subject_phenos
    #         sample_attrs: Sample attributes DataFrame (needed if confounders_df is None)
    #         subject_phenos: Subject phenotypes DataFrame (needed if confounders_df is None)
            
    #     Returns:
    #         Residualized expression DataFrame
    #     """
    #     print("Regressing out confounding factors...")
        
    #     # Generate confounders if not provided
    #     if confounders_df is None:
    #         if sample_attrs is None or subject_phenos is None:
    #             raise ValueError("If confounders_df is None, both sample_attrs and subject_phenos must be provided")
            
    #         samples = expression_df.columns.tolist()
    #         confounders_df = self.generate_confounders_from_metadata(
    #             samples, sample_attrs, subject_phenos
    #         )
        
    #     # Align samples between expression and confounders
    #     common_samples = expression_df.columns.intersection(confounders_df.index)
    #     print(f"Using {len(common_samples)} common samples for regression")
        
    #     # Filter to common samples
    #     expr_aligned = expression_df[common_samples]  # genes × samples
    #     conf_aligned = confounders_df.loc[common_samples]  # samples × confounders
        
    #     # Prepare matrices for vectorized computation
    #     Y = expr_aligned.values.T  # samples × genes (transpose!)
    #     X = conf_aligned.values    # samples × confounders
        
    #     print(f"Matrix dimensions:")
    #     print(f"  Y (expression): {Y.shape} (samples × genes)")
    #     print(f"  X (confounders): {X.shape} (samples × confounders)")
        
    #     # Add intercept term to confounders
    #     X = np.column_stack([np.ones(X.shape[0]), X])
    #     print(f"  X with intercept: {X.shape}")
        
    #     # Check for perfect multicollinearity
    #     try:
    #         # Compute X'X once for all genes
    #         XtX = X.T @ X
    #         print(f"Computing (X'X)^-1 matrix...")
            
    #         # Check condition number
    #         cond_num = np.linalg.cond(XtX)
    #         print(f"Condition number: {cond_num:.2e}")
            
    #         if cond_num > 1e12:
    #             print("Warning: High condition number detected. Adding ridge regularization.")
    #             # Add small ridge regularization to diagonal
    #             ridge_param = 1e-6
    #             XtX += ridge_param * np.eye(XtX.shape[0])
            
    #         # Solve normal equations: β = (X'X)^-1 * X'Y
    #         print("Solving normal equations for all genes simultaneously...")
    #         XtY = X.T @ Y  # confounders × genes
            
    #         # This is the key step - solve for ALL genes at once
    #         beta = np.linalg.solve(XtX, XtY)  # confounders × genes
    #         print(f"Regression coefficients shape: {beta.shape}")
            
    #         # Compute predictions and residuals for all genes
    #         print("Computing residuals...")
    #         Y_pred = X @ beta  # samples × genes
    #         residuals = Y - Y_pred  # samples × genes
            
    #         print(f"Residuals shape: {residuals.shape}")
            
    #         # Convert back to original format (genes × samples)
    #         residual_df = pd.DataFrame(
    #             residuals.T,  # Transpose back to genes × samples
    #             index=expr_aligned.index,
    #             columns=common_samples
    #         )
            
    #         # Compute R² for quality check
    #         ss_res = np.sum(residuals**2, axis=0)  # Sum of squares residual per gene
    #         ss_tot = np.sum((Y - np.mean(Y, axis=0))**2, axis=0)  # Total sum of squares per gene
    #         r2_scores = 1 - (ss_res / ss_tot)
            
    #         print(f"Regression R² statistics:")
    #         print(f"  Mean R²: {np.mean(r2_scores):.4f}")
    #         print(f"  Median R²: {np.median(r2_scores):.4f}")
    #         print(f"  Min R²: {np.min(r2_scores):.4f}")
    #         print(f"  Max R²: {np.max(r2_scores):.4f}")
            
    #         return residual_df
            
    #     except np.linalg.LinAlgError as e:
    #         print(f"Linear algebra error: {e}")
    #         print("Falling back to ridge regression...")
            
    #         # Fallback to ridge regression
    #         from sklearn.linear_model import Ridge
    #         ridge = Ridge(alpha=1.0)
            
    #         residual_matrix = np.zeros_like(Y)
    #         for i in range(Y.shape[1]):
    #             ridge.fit(X, Y[:, i])
    #             residual_matrix[:, i] = Y[:, i] - ridge.predict(X)
            
    #         residual_df = pd.DataFrame(
    #             residual_matrix.T,
    #             index=expr_aligned.index,
    #             columns=common_samples
    #         )
            
    #         return residual_df
    #     # # Transpose for regression (samples x genes) as in the original script
    #     # expr_t = expression_df.T  # samples x genes
        
    #     # # Create DataFrame to hold residuals with same shape as expr_t
    #     # sk_resid_mat = expr_t.copy(deep=True)
    #     # for col in sk_resid_mat.columns:
    #     #     sk_resid_mat[col].values[:] = 0
        
    #     # # Make sure we have common samples between expression and confounders
    #     # common_samples = sorted(set(expr_t.index) & set(confounders_df.index))
    #     # if len(common_samples) == 0:
    #     #     raise ValueError("No common samples between expression data and confounders")
        
    #     # # Filter data to common samples
    #     # expr_t_common = expr_t.loc[common_samples]
    #     # confounders_common = confounders_df.loc[common_samples]
        
    #     # # Following original script approach: confounders should already be processed 
    #     # # (dummy encoded, NaNs handled) in generate_confounders_from_metadata
    #     # X = confounders_common
        
    #     # # Additional safety checks to ensure no NaNs or non-numeric data
    #     # # Check for non-numeric columns and remove them
    #     # numeric_cols = X.select_dtypes(include=['number']).columns
    #     # if len(numeric_cols) < X.shape[1]:
    #     #     print(f"Warning: Removed {X.shape[1] - len(numeric_cols)} non-numeric columns from confounders")
    #     #     X = X[numeric_cols]
        
    #     # # Final comprehensive check for NaN/inf values
    #     # # Check for any NaN or inf values in confounders
    #     # nan_mask = X.isna() | np.isinf(X.replace([np.inf, -np.inf], np.nan))
    #     # if nan_mask.any().any():
    #     #     print("Warning: Found NaN/inf values in confounders.")
    #     #     # Count problematic values per sample and column
    #     #     problematic_samples = nan_mask.any(axis=1)
    #     #     problematic_columns = nan_mask.any(axis=0)
            
    #     #     print(f"Samples with missing values: {problematic_samples.sum()}")
    #     #     print(f"Columns with missing values: {problematic_columns.sum()}")
            
    #     #     # Drop columns with too many missing values (>50% missing)
    #     #     col_missing_pct = nan_mask.mean(axis=0)
    #     #     bad_cols = col_missing_pct[col_missing_pct > 0.5].index
    #     #     if len(bad_cols) > 0:
    #     #         print(f"Dropping {len(bad_cols)} columns with >50% missing values: {list(bad_cols)}")
    #     #         X = X.drop(columns=bad_cols)
    #     #         nan_mask = X.isna() | np.isinf(X.replace([np.inf, -np.inf], np.nan))
            
    #     #     # Drop remaining samples with any missing values
    #     #     if nan_mask.any().any():
    #     #         valid_samples = X.index[~nan_mask.any(axis=1)]
    #     #         print(f"Dropping {len(X) - len(valid_samples)} samples with remaining missing values.")
    #     #         X = X.loc[valid_samples]
    #     #         expr_t_common = expr_t_common.loc[valid_samples]
    #     #         common_samples = valid_samples.tolist()
        
    #     # # Final validation - ensure absolutely no NaN/inf values remain
    #     # assert not X.isna().any().any(), "NaN values still present in confounders after cleaning"
    #     # assert not np.isinf(X).any().any(), "Infinite values present in confounders"
    #     # assert X.shape[0] > 0, "No samples remaining after cleaning confounders"
    #     # assert X.shape[1] > 0, "No confounders remaining after cleaning"
        
    #     # print(f"Using {X.shape[1]} confounders for regression with {len(common_samples)} samples")
    #     # print("Confounders validation passed - no NaN/inf values present")
        
    #     # # Fit linear model for each gene
    #     # for i in range(expr_t_common.shape[1]):  # loop through all genes
    #     #     gene = expr_t_common.columns[i]
    #     #     y = expr_t_common[gene].values  # all samples for one gene
            
    #     #     # Check for NaN/inf in expression values for this gene
    #     #     if np.isnan(y).any() or np.isinf(y).any():
    #     #         print(f"Warning: Gene {gene} has NaN/inf values, skipping regression")
    #     #         sk_resid_mat.loc[common_samples, gene] = y
    #     #         continue
            
    #     #     # Fit linear regression model and get predictions
    #     #     try:
    #     #         reg = LinearRegression()
    #     #         reg.fit(X, y)
    #     #         y_pred = reg.predict(X)
                
    #     #         # Store residuals for common samples
    #     #         sk_resid_mat.loc[common_samples, gene] = y - y_pred
    #     #     except (ValueError, np.linalg.LinAlgError) as e:
    #     #         print(f"Warning: Failed to regress gene {gene}: {str(e)}")
    #     #         # Keep original values for this gene
    #     #         sk_resid_mat.loc[common_samples, gene] = y
        
    #     # # Transpose back to genes x samples and apply quantile normalization
    #     # # as in the old script
    #     # sk_resid_mat_norm = self.quantile_normalize(sk_resid_mat.T)
        
    #     # print("Confounding factors regression completed")
    #     # return sk_resid_mat_norm

    def regress_confounders(self, expression_df, confounders_df=None, 
                        sample_attrs=None, subject_phenos=None):
        """
        Regress out confounding factors using vectorized approach.
        """
        print("Regressing out confounding factors using vectorized approach...")
        
        # Use the vectorized implementation
        return self.regress_confounders_vectorized(
            expression_df, 
            confounders_df=confounders_df,
            sample_attrs=sample_attrs,
            subject_phenos=subject_phenos
        )


    def regress_confounders_vectorized(self, expression_df, confounders_df=None, 
                                    sample_attrs=None, subject_phenos=None):
        """
        Vectorized confounding regression with robust data type handling.
        """
        
        print("Starting vectorized confounding regression...")
        
        # Generate or use provided confounders
        if confounders_df is None:
            if sample_attrs is None or subject_phenos is None:
                raise ValueError("Must provide confounders or metadata")
            
            samples = expression_df.columns.tolist()
            confounders_df = self.generate_confounders_from_metadata(
                samples, sample_attrs, subject_phenos
            )
        
        # Align samples between expression and confounders
        common_samples = expression_df.columns.intersection(confounders_df.index)
        print(f"Using {len(common_samples)} common samples for regression")
        
        if len(common_samples) < 10:
            raise ValueError(f"Too few common samples ({len(common_samples)}) for regression")
        
        # Filter to common samples
        expr_aligned = expression_df[common_samples]  # genes × samples
        conf_aligned = confounders_df.loc[common_samples]  # samples × confounders
        
        print(f"Before cleaning - Confounders shape: {conf_aligned.shape}")
        print(f"Confounders dtypes: {conf_aligned.dtypes.value_counts()}")
        
        # CRITICAL: Ensure ALL columns are numeric
        print("Converting all confounder columns to numeric...")
        numeric_confounders = pd.DataFrame(index=conf_aligned.index)
        
        for col in conf_aligned.columns:
            try:
                # Convert to numeric, coercing errors to NaN
                numeric_col = pd.to_numeric(conf_aligned[col], errors='coerce')
                
                # Check if conversion was successful
                nan_count = numeric_col.isna().sum()
                total_count = len(numeric_col)
                
                if nan_count == total_count:
                    print(f"Dropping column {col}: all values are non-numeric")
                    continue
                elif nan_count > total_count * 0.5:
                    print(f"Dropping column {col}: >50% non-numeric values ({nan_count}/{total_count})")
                    continue
                else:
                    if nan_count > 0:
                        # Fill NaN values with median
                        median_val = numeric_col.median()
                        numeric_col = numeric_col.fillna(median_val)
                        print(f"Column {col}: filled {nan_count} NaN values with median {median_val:.4f}")
                    
                    numeric_confounders[col] = numeric_col
                    
            except Exception as e:
                print(f"Error processing column {col}: {e}")
                continue
        
        print(f"After cleaning - Confounders shape: {numeric_confounders.shape}")
        
        if numeric_confounders.shape[1] == 0:
            raise ValueError("No valid numeric confounders remaining after cleaning")
        
        # Check for any remaining non-numeric data
        non_numeric_cols = numeric_confounders.select_dtypes(include=['object']).columns
        if len(non_numeric_cols) > 0:
            print(f"Warning: Found remaining non-numeric columns: {list(non_numeric_cols)}")
            numeric_confounders = numeric_confounders.drop(columns=non_numeric_cols)
        
        # Final validation - ensure all values are finite
        print("Checking for infinite values...")
        for col in numeric_confounders.columns:
            inf_mask = np.isinf(numeric_confounders[col])
            if inf_mask.any():
                print(f"Column {col}: replacing {inf_mask.sum()} infinite values with median")
                median_val = numeric_confounders[col][~inf_mask].median()
                numeric_confounders.loc[inf_mask, col] = median_val
        
        # Prepare matrices for vectorized computation
        Y = expr_aligned.values.T.astype(np.float64)  # samples × genes, ensure float64
        X = numeric_confounders.values.astype(np.float64)  # samples × confounders, ensure float64
        
        print(f"Final matrix dimensions:")
        print(f"  Y (expression): {Y.shape} (samples × genes)")
        print(f"  X (confounders): {X.shape} (samples × confounders)")
        print(f"  Y dtype: {Y.dtype}")
        print(f"  X dtype: {X.dtype}")
        
        # Check for any remaining NaN values
        if np.isnan(Y).any():
            print("Warning: Found NaN values in expression matrix")
            raise ValueError("Expression matrix contains NaN values")
        
        if np.isnan(X).any():
            print("Warning: Found NaN values in confounders matrix")
            raise ValueError("Confounders matrix contains NaN values")
        
        # Add intercept term to confounders
        X = np.column_stack([np.ones(X.shape[0], dtype=np.float64), X])
        print(f"  X with intercept: {X.shape}")
        
        try:
            # Compute X'X once for all genes
            XtX = X.T @ X
            print(f"Computing (X'X)^-1 matrix...")
            
            # Check condition number
            cond_num = np.linalg.cond(XtX)
            print(f"Condition number: {cond_num:.2e}")
            
            if cond_num > 1e12:
                print("Warning: High condition number detected. Adding ridge regularization.")
                ridge_param = 1e-6
                XtX += ridge_param * np.eye(XtX.shape[0])
            
            # Solve normal equations: β = (X'X)^-1 * X'Y
            print("Solving normal equations for all genes simultaneously...")
            XtY = X.T @ Y  # confounders × genes
            
            # This is the key step - solve for ALL genes at once
            beta = np.linalg.solve(XtX, XtY)  # confounders × genes
            print(f"Regression coefficients shape: {beta.shape}")
            
            # Compute predictions and residuals for all genes
            print("Computing residuals...")
            Y_pred = X @ beta  # samples × genes
            residuals = Y - Y_pred  # samples × genes
            
            print(f"Residuals shape: {residuals.shape}")
            
            # Convert back to original format (genes × samples)
            residual_df = pd.DataFrame(
                residuals.T,  # Transpose back to genes × samples
                index=expr_aligned.index,
                columns=common_samples
            )
            
            # Compute R² for quality check
            ss_res = np.sum(residuals**2, axis=0)  # Sum of squares residual per gene
            ss_tot = np.sum((Y - np.mean(Y, axis=0))**2, axis=0)  # Total sum of squares per gene
            r2_scores = 1 - (ss_res / ss_tot)
            
            print(f"Regression R² statistics:")
            print(f"  Mean R²: {np.mean(r2_scores):.4f}")
            print(f"  Median R²: {np.median(r2_scores):.4f}")
            print(f"  Min R²: {np.min(r2_scores):.4f}")
            print(f"  Max R²: {np.max(r2_scores):.4f}")
            
            print("Vectorized confounding regression completed successfully!")
            return residual_df
            
        except np.linalg.LinAlgError as e:
            print(f"Linear algebra error: {e}")
            print("Falling back to ridge regression...")
            
            # Fallback to ridge regression
            from sklearn.linear_model import Ridge
            ridge = Ridge(alpha=1.0)
            
            residual_matrix = np.zeros_like(Y)
            for i in range(Y.shape[1]):
                ridge.fit(X, Y[:, i])
                residual_matrix[:, i] = Y[:, i] - ridge.predict(X)
            
            residual_df = pd.DataFrame(
                residual_matrix.T,
                index=expr_aligned.index,
                columns=common_samples
            )
            
            return residual_df

    def preprocess_expression_data(self,
                                    expression_df: pd.DataFrame,
                                    hgnc_data: Optional[pd.DataFrame] = None,
                                    confounders_df: Optional[pd.DataFrame] = None,
                                    sample_attrs: Optional[pd.DataFrame] = None,
                                    subject_phenos: Optional[pd.DataFrame] = None,
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
            confounders_df: Confounding factors DataFrame (optional, will be generated if None)
            sample_attrs: Sample attributes DataFrame (required if confounders_df is None and regress_confounders is True)
            subject_phenos: Subject phenotypes DataFrame (required if confounders_df is None and regress_confounders is True)
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
        if regress_confounders:
            if confounders_df is None:
                if sample_attrs is None or subject_phenos is None:
                    raise ValueError("If confounders_df is None and regress_confounders is True, both sample_attrs and subject_phenos must be provided")
                print("No confounders provided - generating from GTEx metadata following original script methodology")
            
            processed_df = self.regress_confounders(
                processed_df, 
                confounders_df=confounders_df,
                sample_attrs=sample_attrs,
                subject_phenos=subject_phenos
            )
            processing_info['steps_applied'].append('confounder_regression')
            processing_info['after_confounder_regression'] = processed_df.shape
        
        processing_info['final_shape'] = processed_df.shape
        
        print(f"Preprocessing completed. Final shape: {processed_df.shape}")
        print(f"Steps applied: {processing_info['steps_applied']}")
        
        return processed_df, processing_info
