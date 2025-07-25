"""
GTEx Data Reader Module
Handles reading and parsing of GTEx GCT files and sample metadata
"""

import pandas as pd
import numpy as np
import io
from typing import Tuple, Dict, List, Optional
from pathlib import Path

class GTExDataReader:
    """Class for reading and processing GTEx expression data and metadata."""
    
    def __init__(self, data_dir: str = None):
        """
        Initialize the GTEx data reader.
        
        Args:
            data_dir: Path to directory containing GTEx data files. 
                     If None, will search for data in common locations.
        """
        if data_dir is None:
            # Try to find data directory automatically
            possible_paths = [
                "../data",           # When running from new_script/
                "data",              # When running from project root
                "./data",            # Explicit current directory
                "new_script/../data" # When running from parent directory
            ]
            
            for path in possible_paths:
                test_path = Path(path)
                if test_path.exists() and test_path.is_dir():
                    # Check if it contains GTEx files
                    gct_files = list(test_path.glob("*gene_tpm.gct"))
                    if gct_files:
                        self.data_dir = test_path
                        print(f"Found data directory: {self.data_dir.absolute()}")
                        break
            else:
                # Default fallback
                self.data_dir = Path("../data")
                print(f"Warning: Using default data directory: {self.data_dir.absolute()}")
        else:
            self.data_dir = Path(data_dir)
        
    def read_gct_file(self, gct_file: str, sample_limit: Optional[int] = None) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        Read GTEx GCT format file efficiently.
        
        Args:
            gct_file: Path to GCT file
            sample_limit: Optional limit on number of samples to read (for testing)
            
        Returns:
            tuple: (expression_data, gene_metadata)
        """
        gct_path = self.data_dir / gct_file
        
        print(f"Reading GCT file: {gct_path}")
        
        # First pass: read just the header to get dimensions
        with open(gct_path, 'r') as f:
            version = f.readline().strip()  # #1.2
            dimensions = f.readline().strip().split('\t')
            n_genes = int(dimensions[0])
            n_samples = int(dimensions[1])
            
            # Read header line with sample IDs
            header_line = f.readline().strip()
            columns = header_line.split('\t')
            
            # First two columns are gene ID and description
            gene_id_col = columns[0]
            gene_desc_col = columns[1]
            sample_ids = columns[2:]
            
            if sample_limit:
                sample_ids = sample_ids[:sample_limit]
                n_samples = len(sample_ids)
        
        print(f"GCT version: {version}")
        print(f"Dimensions: {n_genes} genes, {n_samples} samples")
        print(f"Preparing to read {n_samples} samples")
        
        # Preallocate numpy array for better memory efficiency
        expr_matrix = np.zeros((n_genes, n_samples))
        gene_ids = []
        gene_descriptions = []
        
        # Second pass: read the data efficiently
        with open(gct_path, 'r') as f:
            # Skip header lines
            for _ in range(3):
                f.readline()
            
            # Read data in chunks
            chunk_size = 1000
            genes_read = 0
            
            while genes_read < n_genes:
                chunk_lines = []
                for _ in range(min(chunk_size, n_genes - genes_read)):
                    line = f.readline()
                    if not line:
                        break
                    chunk_lines.append(line)
                
                if not chunk_lines:
                    break
                    
                # Process the chunk
                for i, line in enumerate(chunk_lines):
                    parts = line.strip().split('\t')
                    gene_id = parts[0]
                    gene_desc = parts[1]
                    expression_values = parts[2:2+n_samples] if not sample_limit else parts[2:2+n_samples]
                    
                    # Convert to float
                    try:
                        expr_numeric = [float(x) for x in expression_values]
                        gene_ids.append(gene_id)
                        gene_descriptions.append(gene_desc)
                        expr_matrix[genes_read + i, :] = expr_numeric
                    except ValueError:
                        print(f"Warning: Could not parse gene {gene_id}, skipping...")
                        # Fill with zeros
                        gene_ids.append(gene_id)
                        gene_descriptions.append(gene_desc)
                        # We keep the zeros that were preallocated
                
                genes_read += len(chunk_lines)
                print(f"Processed {genes_read}/{n_genes} genes ({genes_read/n_genes*100:.1f}%)...")
        
        print(f"Finished reading {len(gene_ids)} genes.")
        
        # Create DataFrames - only use valid rows
        expression_df = pd.DataFrame(
            expr_matrix[:len(gene_ids)],
            index=gene_ids,
            columns=sample_ids
        )
        
        gene_metadata = pd.DataFrame({
            'gene_id': gene_ids,
            'description': gene_descriptions
        }, index=gene_ids)
        
        print(f"Loaded expression data: {expression_df.shape}")
        return expression_df, gene_metadata
    
    def read_sample_attributes(self, attr_file: str = "GTEx_Analysis_v8_Annotations_SampleAttributesDS.tsv") -> pd.DataFrame:
        """Read GTEx sample attributes file."""
        attr_path = self.data_dir / attr_file
        
        print(f"Reading sample attributes: {attr_path}")
        sample_attrs = pd.read_csv(attr_path, sep='\t', low_memory=False)
        
        # Set SAMPID as index for easy lookup
        sample_attrs.set_index('SAMPID', inplace=True)
        
        print(f"Loaded sample attributes: {sample_attrs.shape}")
        return sample_attrs
    
    def read_subject_phenotypes(self, pheno_file: str = "GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.tsv") -> pd.DataFrame:
        """Read GTEx subject phenotypes file."""
        pheno_path = self.data_dir / pheno_file
        
        print(f"Reading subject phenotypes: {pheno_path}")
        subject_phenos = pd.read_csv(pheno_path, sep='\t', low_memory=False)
        
        # Set SUBJID as index
        subject_phenos.set_index('SUBJID', inplace=True)
        
        print(f"Loaded subject phenotypes: {subject_phenos.shape}")
        return subject_phenos
    
    def read_protein_coding_genes(self, hgnc_file: str = "hgnc_complete_set.tsv") -> pd.DataFrame:
        """Read HGNC protein coding genes file."""
        hgnc_path = self.data_dir / hgnc_file
        
        print(f"Reading HGNC data: {hgnc_path}")
        hgnc_data = pd.read_csv(hgnc_path, sep='\t', low_memory=False)
        
        # Filter for protein coding genes
        protein_coding = hgnc_data[
            (hgnc_data['locus_type'] == 'gene with protein product') &
            (hgnc_data['status'] == 'Approved')
        ].copy()
        
        print(f"Found {len(protein_coding)} protein coding genes")
        return protein_coding
    
    def get_tissue_samples(self, sample_attrs: pd.DataFrame, tissue_name: str) -> List[str]:
        """
        Get sample IDs for a specific tissue.
        
        Args:
            sample_attrs: Sample attributes DataFrame
            tissue_name: Name of tissue (SMTSD column)
            
        Returns:
            List of sample IDs for the tissue
        """
        tissue_samples = sample_attrs[sample_attrs['SMTSD'] == tissue_name].index.tolist()
        print(f"Found {len(tissue_samples)} samples for tissue: {tissue_name}")
        return tissue_samples
    
    def get_available_tissues(self, sample_attrs: pd.DataFrame) -> List[str]:
        """Get list of available tissues."""
        tissues = sample_attrs['SMTSD'].unique().tolist()
        tissues.sort()
        return tissues
    
    def filter_samples_by_metadata(self, 
                                   sample_attrs: pd.DataFrame,
                                   subject_phenos: pd.DataFrame,
                                   min_rin: float = 5.7,
                                   age_group: Optional[str] = None,
                                   sex: Optional[str] = None,
                                   tissue: Optional[str] = None) -> List[str]:
        """
        Filter samples based on quality and phenotype criteria.
        
        Args:
            sample_attrs: Sample attributes DataFrame
            subject_phenos: Subject phenotypes DataFrame  
            min_rin: Minimum RIN score
            age_group: Age group filter options:
                      - 'young': 20-29 and 30-39 age groups
                      - 'old': 60-69 and 70-79 age groups  
                      - list: Custom list of age ranges (e.g., ['20-29', '50-59'])
                      - None: No age filtering
            sex: Sex filter (1=male, 2=female, or None)
            tissue: Tissue type filter
            
        Returns:
            List of filtered sample IDs
        """
        # Start with all samples
        filtered_samples = sample_attrs.copy()
        
        # Filter by RIN score
        if 'SMRIN' in filtered_samples.columns:
            filtered_samples = filtered_samples[filtered_samples['SMRIN'] >= min_rin]
            print(f"After RIN >= {min_rin}: {len(filtered_samples)} samples")
        
        # Filter by tissue
        if tissue:
            filtered_samples = filtered_samples[filtered_samples['SMTSD'] == tissue]
            print(f"After tissue filter ({tissue}): {len(filtered_samples)} samples")
        
        # Filter by subject phenotypes
        if age_group or sex:
            # Extract subject IDs from sample IDs (GTEX-XXXXX from GTEX-XXXXX-XXXX-XX-XXXXX)
            # Use str.extract with expand=False to get a Series instead of DataFrame
            filtered_samples = filtered_samples.copy()  # Avoid SettingWithCopyWarning
            filtered_samples['SUBJID'] = filtered_samples.index.str.extract(r'(GTEX-[^-]+)', expand=False)
            
            # Debug: Check if extraction worked
            extracted_count = filtered_samples['SUBJID'].notna().sum()
            print(f"   Extracted {extracted_count} subject IDs from {len(filtered_samples)} samples")
            
            if extracted_count == 0:
                print(f"   WARNING: No subject IDs extracted. Sample format may be different.")
                print(f"   Sample ID examples: {list(filtered_samples.index[:3])}")
                return []
            
            # Remove samples without valid subject IDs
            filtered_samples = filtered_samples.dropna(subset=['SUBJID'])
            print(f"   Samples with valid subject IDs: {len(filtered_samples)}")
            
            # Merge with phenotypes
            filtered_with_pheno = filtered_samples.merge(
                subject_phenos, 
                left_on='SUBJID', 
                right_index=True, 
                how='inner'
            )
            
            print(f"   Samples with phenotype data: {len(filtered_with_pheno)}")
            
            # Age group filter
            if age_group:
                if age_group == 'young':
                    # Include young adults: 20-29 and 30-39
                    filtered_with_pheno = filtered_with_pheno[filtered_with_pheno['AGE'].isin(['20-29', '30-39','40-49', '50-59'])]
                elif age_group == 'old':
                    # Include older adults: 60-69 and 70-79
                    filtered_with_pheno = filtered_with_pheno[filtered_with_pheno['AGE'].isin(['60-69', '70-79'])]
                elif isinstance(age_group, list):
                    # Allow custom age group list
                    filtered_with_pheno = filtered_with_pheno[filtered_with_pheno['AGE'].isin(age_group)]
                print(f"After age group filter ({age_group}): {len(filtered_with_pheno)} samples")
            
            # Sex filter
            if sex:
                filtered_with_pheno = filtered_with_pheno[filtered_with_pheno['SEX'] == sex]
                print(f"After sex filter ({sex}): {len(filtered_with_pheno)} samples")
            
            return filtered_with_pheno.index.tolist()
        
        return filtered_samples.index.tolist()
    
    def get_available_age_groups(self, subject_phenos: pd.DataFrame) -> List[str]:
        """Get list of available age groups in the dataset."""
        age_groups = subject_phenos['AGE'].dropna().unique().tolist()
        age_groups.sort()
        return age_groups
    
    def get_age_group_counts(self, subject_phenos: pd.DataFrame) -> Dict[str, int]:
        """Get counts of subjects in each age group."""
        return subject_phenos['AGE'].value_counts().to_dict()
    
    def get_available_sexes(self, subject_phenos: pd.DataFrame) -> Dict[int, str]:
        """Get available sex categories with their meanings."""
        sex_mapping = {1: 'Male', 2: 'Female'}
        available_sexes = {}
        for sex_code in subject_phenos['SEX'].dropna().unique():
            if sex_code in sex_mapping:
                available_sexes[sex_code] = sex_mapping[sex_code]
        return available_sexes
