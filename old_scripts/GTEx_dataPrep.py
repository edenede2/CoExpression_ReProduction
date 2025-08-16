#imports
import pandas as pd
import os 
from pathlib import Path
import re
import copy
import pickle

import numpy as np
import scipy as stats
import random
from rich import print as rprint
from rich.progress import track

import seaborn as sns
import matplotlib.pyplot as plt

import scipy as sp
from scipy.stats import chi2
from sklearn.covariance import MinCovDet

from sklearn.decomposition import PCA
from sklearn.decomposition import TruncatedSVD

from sklearn.linear_model import LinearRegression
from sklearn.feature_selection import f_regression

'''
Step 1: Filter out the non-protein coding genes from GTEX

Remove genes from the data which do not encode a protein and to log transform the gene expression matrix.

Source:  https://www.genenames.org/download/statistics-and-files/ retrieved January 19 2022 </p>
Selected 19208 protein coding genes as txt file </p>
Input: protein coding gene list with current list of know protein coding genes
       GTEX matrix contain at total of 56200 genes and 17382 samples 

Output: 
A matrix containing on the ENSMBL description and only the genes that were in the intersection of the protein coding genes and GTEx matrix. 
A dictionary from ENSMBL to description.
'''

def filter_non_protein_coding(HGNC_list, GTEx_tpm):
    
    HGNC_list.rename(columns={'symbol':'PCG Symbols'},inplace=True)
    PCG_list = HGNC_list['PCG Symbols'].tolist()
    pcg_symbols = set(HGNC_list['PCG Symbols'].dropna().tolist())

    df = GTEx_tpm[GTEx_tpm['Description'].isin(pcg_symbols)].copy()

    df['ensembel_id']

    matrix_pcg_df = GTEx_tpm[GTEx_tpm['Description'].isin(PCG_list)]
    pcg_dictionary = matrix_pcg_df[['Name','Description']]
    pcg_dictionary = pcg_dictionary.set_index('Name')
    
    matrix_pcg_df =matrix_pcg_df.drop(['Description'], axis=1)
    matrix_pcg_df.set_index('Name',inplace=True)
    
    log_transform(matrix_pcg_df)
    
    return pcg_dictionary,matrix_pcg_df

    

def log_transform(matrix_pcg_df):
    matrix_pcg_df = matrix_pcg_df +1
    matrix_pcg_df = np.log2(matrix_pcg_df)
    return matrix_pcg_df

'''
Step 2: Sample attributes
1. Filter out death types for both samples and attributes 
2. Filter out low RIN

Input: Sample annotations (attributes csv file),
       Subject phenotypes (pheno csv file),
       GTEX matrix
'''

def setup_attributes (pheno_df, attributes_df, matrix_pcg_df):
    # Filter death types
    pheno_a = pheno_df[(pheno_df['DTHHRDY'] > 0) & (pheno_df['DTHHRDY'] < 3)]
    important_attributes = attributes_df[['SMRIN','SMTSISCH','SMTSD','SMGEBTCH']]
    important_attributes = important_attributes.assign(SUBJID = '')
    sampids = important_attributes.index.tolist()
    for i in range(important_attributes.shape[0]):
        x = sampids[i].split('-')
        xx = ''
        xx = xx+x[0]+'-'+x[1]
        important_attributes.iloc[i,4] = xx
    important_attributes.dropna(inplace = True)
    
    # Filter RIN
    important_attributes = important_attributes.drop(important_attributes[important_attributes['SMRIN'] <5.7].index, axis = 0)

    new_df  =pd.merge(important_attributes.reset_index(), pheno_a, on=['SUBJID'], how='inner')
    new_df =new_df.set_index('SAMPID')
    pheno_aa = new_df.dropna(axis=0,how='any')
    pheno_aa['AGE'] = pheno_aa['AGE'].str.split('-', n=1).str[0].astype(int)

    aa_samples = pheno_aa.index.tolist()
    matrix_pcg_aa = matrix_pcg_df.loc[:, matrix_pcg_df.columns.isin(aa_samples)]
    aa_matrix_samples = matrix_pcg_aa.columns.to_list()
    sample_attributes = new_df.loc[new_df.index.isin(aa_matrix_samples),:]
    return matrix_pcg_aa, sample_attributes

'''
Step 3: Create data for tissues
'''
def create_data_for_tissue (sample_attributes, matrix_pcg_aa, number_of_tissues = 13, length_threshold = 100):

    samples = sample_attributes
    sample_tissues = samples['SMTSD'].value_counts()

    tissues = sample_tissues[sample_tissues>100].index.to_list()
    pcg = matrix_pcg_aa



    for i in range(number_of_tissues):

        current_tissue = tissues[i]
        current_tissue = current_tissue[:length_threshold] #@TODO: Check why to limit the length of the tissue name

        if not os.path.exists(f'old_scripts/new_outputs/{current_tissue}'):
            os.makedirs(f'old_scripts/new_outputs/{current_tissue}')
        tissue_type = current_tissue[:length_threshold] #@TODO: Check why to limit the length of the tissue name again

        tissue_samples = samples[samples['SMTSD'].str.startswith(tissue_type)]
        tissue_ids = tissue_samples.index.tolist()

        # print(pheno.shape)
        matrix = pcg[tissue_ids]

        tissue_samples.to_csv(f'old_scripts/new_outputs/{current_tissue}/tissue_sample_DT1DT2.csv')
        matrix.to_csv(f'old_scripts/new_outputs/{current_tissue}/tissue_matrix_DT1DT2.csv')

    
'''
Step 4: Process tissue data

1. Filter genes with low expression and low veriabilty 
2. Remove outliers
3. Quantile normalization
4. Make modifications to phenotype attributes
'''

def create_exclusion_list(genes, values, limit_val = 0.02):
    exclude_genes_list = []
    for i in range(len(genes)):
        if (values[i] < limit_val):
            exclude_genes_list.append(genes[i])
    return exclude_genes_list 


# filter Genes with zero variance  
def filter_genes (tissue_matrix, sample_attributes, variability_threshold):
    sample_attributes['AGE'] = sample_attributes['AGE'].str.split('-', n=1).str[0].astype(int)

    print(str(len(tissue_matrix[(tissue_matrix.T <np.log2(0.1+1) ).sum()>0.2*tissue_matrix.shape[1]])) + " genes filtered out")
    df1 = tissue_matrix[(tissue_matrix.T <np.log2(0.1+1) ).sum()<0.2*tissue_matrix.shape[1]] 
    df1 = df1.T
    variability_df = df1.var()
    print(variability_df)
    gene_names = variability_df.index.tolist()
    gene_vars = variability_df.tolist()

    excluded_list = create_exclusion_list(gene_names,gene_vars,variability_threshold)
    print('length of excluded list:' + str(len(excluded_list)))

    df1 = df1.drop(excluded_list, axis=1)
    sample1 = sample_attributes.copy(deep=True)
    return df1, sample1
  

#remove outlier samples using mahanalobis on specific pca
def sample_outliers_df(matrix_val, sample_val):
    print('sample outliers')
    print(matrix_val.shape) # rows are genes col are sampels
    print(matrix_val.head())
    print(sample_val.head())
    print('isna')
    print(matrix_val.isna().sum())
    matrix_val.isna().sum()
    matrix_val.fillna(value=0.0000001,inplace=True)
    print(matrix_val.isna().sum())
    
    sample_cum_pca = []
    gene_cum_pca = []
    current_sum = 0
    print("removing outliers")
    print(matrix_val.shape)
    #pca_gene = PCA(n_components=gene_pca_components)
    pca_gene = TruncatedSVD(n_components=20, random_state=1001)
    pca_gene.fit(matrix_val)
    components =  pca_gene.components_
    components = components.T
    print(type(components))
    print('components shape')
    print(components.shape)
    comp_df = pd.DataFrame(data=components)
    
    samples_before_outlier_removal = matrix_val.columns.tolist()
    print(samples_before_outlier_removal)
    samples_to_remove = []
    
     
    # Covariance matrix
    covariance  = np.cov(components , rowvar=False)
    # Covariance matrix power of -1
    covariance_pm1 = np.linalg.matrix_power(covariance, -1)
    # Center point
    centerpoint = np.mean(components , axis=0)
    
    # Distances between center point and 
    distances = []
    for i, val in enumerate(components):
        p1 = val
        
        p2 = centerpoint
       
        distance = (p1-p2).T.dot(covariance_pm1).dot(p1-p2)
        distances.append(distance)
    distances = np.array(distances)

    
    # Cutoff (threshold) value from Chi-Sqaure Distribution for detecting outliers 
    cutoff = chi2.ppf(0.99, components.shape[1]*2)
    print('distances')
    print(distances)
    print('cutoff')
    print(cutoff)
    print(matrix_val.head())
    print(matrix_val.shape)
    print(sample_val.head())
    print(sample_val.shape)

    # Index of outliers
    outlierIndexes = np.where(distances > cutoff )
    outlierIndexes = list(outlierIndexes[0])

    print('--- Index of Outliers ----')
    print(outlierIndexes)
    print(type(outlierIndexes))
    print(len(outlierIndexes))
    # for q2 in range(len(outlierIndexes)):
    #     print(outlierIndexes[q2])
    #     samples_to_remove.append(samples_before_outlier_removal[q2])
    for q2 in outlierIndexes:
        print(q2)
        samples_to_remove.append(samples_before_outlier_removal[q2])
    print(samples_to_remove)
    
    df = pd.DataFrame(data=components)
    print(df.head())

    print('dropping outliers')
    print('matrix shape before drop')
    print(matrix_val.shape)
    matrix_val = matrix_val.drop(labels = samples_to_remove, axis=1)
    print('matrix shape after drop')
    print(matrix_val.shape)
    print('sample shape before drop')
    print(sample_val.shape)
    sample_val = sample_val.drop(labels = samples_to_remove, axis=0)
    print('sample shape after drop')
    print(sample_val.shape)

    return matrix_val, sample_val


def quantile_normalize(matrix_val):
    df = matrix_val
    df_image1 = df.iloc[:,0:3]
    print('df image before quantile normalization: '+str(df_image1.shape))
    sns.set(rc={'figure.figsize':(11.7,8.27)})
    ax = sns.boxplot(data=df_image1, linewidth=2.5).set(title='Before quantile')
    plt.savefig('before_quantile.png')
    # plt.show()
    print(matrix_val.shape)
    print(df_image1.head(2))

    df_sorted = pd.DataFrame(np.sort(df.values,
                                     axis=0), 
                             index=df.index, 
                             columns=df.columns)
    df_mean = df_sorted.mean(axis=1)
    df_mean.index = np.arange(1, len(df_mean) + 1)
    df_qn =df.rank(method="min").stack().astype(int).map(df_mean).unstack()
    
    
    df_image2 = df_qn.iloc[:,0:3]
    print('df image after quantile normalization: '+str(df_image2.shape))
    sns.set(rc={'figure.figsize':(11.7,8.27)})
    ax = sns.boxplot(data=df_image2, linewidth=2.5).set(title='After quantile')
    plt.savefig('after_quantile.png')
    # plt.show()
    print(df_image2.head(2))
    
    return df_qn



def process_tissue (tissue_matrix , sample_attributes, current_tissue, variability_threshold):

    # 1.filter genes with low expression and low veriabilty
    print("---1.filtering---")
    ###remove genes that have value less than 0.1 transcripts per million (TPM) in more than 80% of the samples
    df1, sample1 = filter_genes(tissue_matrix, sample_attributes, variability_threshold)


    # 2. remove outliers
    print("---2.remove outliers---")
    df2, sample2 = sample_outliers_df(df1.T,sample1)

    # 3. quantile normalization
    print("---3.quantile normalization---")

    df3 = quantile_normalize(df2)

    # 4. make modifications to phenotype attributes

    sample3 = sample2.copy(deep=True)
    # agegroup = sample3['AGE'].tolist()
    # sample3['AGE'] = [elem[:2] for elem in agegroup]
    SMGEBTCH= sample3['SMGEBTCH'].value_counts(normalize=False, sort=True)
    SMGEBTCH_names = SMGEBTCH.index.tolist()
    SMGEBTCH_value =  SMGEBTCH.tolist()
    SMGEBTCH_new_list  = []
    for i in range (len(SMGEBTCH_names)):
        if SMGEBTCH_value[i] > 1:
            SMGEBTCH_new_list.append(SMGEBTCH_names[i])

    batches = []
    for i in range(sample3.shape[0]):

        if sample3['SMGEBTCH'].iloc[i] in SMGEBTCH_new_list:
            batches.append(sample3['SMGEBTCH'].iloc[i])
        else:
            batches.append('ASINGLETON_SMGEBTCH')

    sample3['SMGEBTCH'] = batches
    sample4 = sample3.copy(deep=True)
    SMGEBTCH = pd.get_dummies(sample4['SMGEBTCH'],drop_first=True)
    sample4.drop(['SMTSD','SMGEBTCH','SUBJID'], axis=1, inplace=True)

    new_sample4 = [sample3, SMGEBTCH]
    sample5 = pd.concat(new_sample4 , join='inner', axis=1)


    df3.to_csv(f'old_scripts/new_outputs/{current_tissue}/tissue_matrix_preReg.csv')
    sample5.to_csv(f'old_scripts/new_outputs/{current_tissue}/tissue_sample_preReg.csv')
    return df3, sample5

'''
Step 5: Confouding factors
'''

def regress_confounding (matrix , sample, current_tissue):
    matrix = matrix.T # sampels X genes 
    sk_resid_mat=matrix.copy(deep=True) # sampels X genes 
    for col in sk_resid_mat.columns: # loop on genes
        sk_resid_mat[col].values[:] = 0

    age =  sample['AGE']
    death_type =  sample['DTHHRDY']

    x = sample.drop(['SUBJID','AGE', 'SMRIN', 'DTHHRDY','SMTSD', 'SMGEBTCH'], axis =1) #@TODO: The non numeric columns didn't droped in the original script

    for col in x.columns: #@TODO: The bool columns were not converted to int [1,0] in the original script
        if x[col].dtype == 'bool':
            x[col] = x[col].astype(int)
    # Add process bar
    rprint(f"[bold cyan]Regressing confounding factors for tissue: {current_tissue}[/bold cyan]")
    for i in track(range(matrix.shape[1]), description="Processing genes..."): # loop on cols = genes
        reg = LinearRegression()
        y = matrix.iloc[:,i] # all sampels, one gene
        reg.fit(x,y)
        prediction = reg.predict(x)
        coefs = reg.coef_.tolist()
        sk_resid_mat.iloc[:,i] = y - prediction # replace gene i with residual

    # for i in range(matrix.shape[1]): # loop on cols = genes
    #     reg = LinearRegression()
    #     y = matrix.iloc[:,i] # all sampels, one gene
    #     reg.fit(x,y)
    #     prediction = reg.predict(x)
    #     coefs = reg.coef_.tolist()
    #     sk_resid_mat.iloc[:,i] = y - prediction # replace gene i with residual
              
    feature_names = sample.columns.tolist()

    sk_resid_mat_norm = quantile_normalize(sk_resid_mat.T)
    sk_resid_mat_norm.to_csv(f'old_scripts/new_outputs/{current_tissue}/sk_residual_matrix_w_const.csv')
    return sk_resid_mat_norm

    
'''
Step 6: split age groups
'''

def split_age_group(matrix, sample_processed, current_tissue):

    # matrix = matrix.T #TODO : why transpose here?
    # pool_of_samples = sample_processed.index.tolist() 
    # Matrix is genes x samples, coming from regress_confounding
    
    if sample_processed['AGE'].dtype == object:
        sample_processed = sample_processed.copy()
        sample_processed['AGE'] = pd.to_numeric(sample_processed['AGE'], errors='coerce')
        if sample_processed['AGE'].isnull().any():
            rprint(f"[bold cyan][split_age_group] Warning: Non-numeric values found in 'AGE' column for tissue {current_tissue}. Converting to numeric with NaN for non-convertible values.[/bold cyan]")


    pool_of_samples = sample_processed.index.tolist()
    overlap_samples = matrix.columns.intersection(pool_of_samples)
    
    if overlap_samples.empty:
        rprint(f"[bold yellow][split_age_group] No overlap between matrix and sample for tissue {current_tissue}.[/bold yellow]")
        rprint(f"[bold yellow]Matrix columns:[/bold yellow] {matrix.columns.tolist()}")
        rprint(f"[bold yellow]Sample index:[/bold yellow] {sample_processed.index.tolist()}")


    young_sample =  sample_processed[sample_processed['AGE'] < 60]
    old_sample =    sample_processed[sample_processed['AGE'] >= 60]  
    
    young_list =  young_sample.index.tolist()
    old_list =  old_sample.index.tolist()
    
    young_matrix = matrix[matrix.columns.intersection(young_list)].T

    old_matrix = matrix[matrix.columns.intersection(old_list)].T
    # old_matrix =  old_matrix.T
    
    young_sample.to_csv(f'old_scripts/new_outputs/{current_tissue}/young_sample.csv')
    young_matrix.to_csv(f'old_scripts/new_outputs/{current_tissue}/young_matrix.csv')
    old_sample.to_csv(f'old_scripts/new_outputs/{current_tissue}/old_sample.csv')
    old_matrix.to_csv(f'old_scripts/new_outputs/{current_tissue}/old_matrix.csv')

    return young_sample, young_matrix, old_sample, old_matrix
 


if __name__ == '__main__':
    
    # with open("expression_df.pkl", "rb") as f:
    #     GTEx_tpm = pickle.load(f)
    GTEx_tpm = pd.read_csv('data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct', sep='\t',skiprows=2)
    HGNC_list = pd.read_csv('data/hgnc_complete_set.tsv', sep='\t')
    pheno_df = pd.read_csv('data/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.tsv', sep='\t', index_col=0)
    attributes_df = pd.read_csv('data/GTEx_Analysis_v8_Annotations_SampleAttributesDS.tsv', sep='\t', index_col=0)

    number_of_tissues = 3
    length_threshold = 100
    variability_threshold = 0.02
    sample_pca_components = 3
    gene_pca_components = 50

    tissues = ['Muscle - Skeletal',
    #  'Whole Blood',
    #  'Skin - Sun Exposed (Lower leg)',
    #  'Skin - Not Sun Exposed (Suprapubic)',
     'Adipose - Subcutaneous',
    #  'Thyroid',
    #  'Artery - Tibial',
    #  'Nerve - Tibial',
    #  'Lung',
    #  'Brain - Cerebellum',
    #  'Heart - Atrial Appendage',
     'Brain - Cortex']
    #  'Adipose - Visceral (Omentum)']


    pcg_dictionary,matrix_pcg_df = filter_non_protein_coding(HGNC_list, GTEx_tpm)
    # matrix_pcg_df = log_transform(matrix_pcg_df) # @TODO : Check why log transformation is applied twice in the original script

    matrix_pcg_aa, sample_attributes = setup_attributes(pheno_df, attributes_df, matrix_pcg_df)
    create_data_for_tissue(sample_attributes, matrix_pcg_aa, number_of_tissues, length_threshold )


    for i in range(number_of_tissues):
        current_tissue = tissues[i] 
        # rprint("**********"+current_tissue+"************")
        rprint(f"[bold green]Processing tissue: {current_tissue}[/bold green]")
        current_tissue = current_tissue[:length_threshold]
        tissue_matrix = pd.read_csv(f'old_scripts/new_outputs/{current_tissue}/tissue_matrix_DT1DT2.csv', index_col=0)
        sample_attributes = pd.read_csv(f'old_scripts/new_outputs/{current_tissue}/tissue_sample_DT1DT2.csv', index_col=0)
        rprint(f"[bold blue]Loaded tissue matrix and sample attributes for {current_tissue}[/bold blue]")
        tissue_matrix_preReg, tissue_sample_preReg = process_tissue (tissue_matrix , sample_attributes, current_tissue, variability_threshold)
        rprint(f"[bold magenta]Processed tissue data for {current_tissue}[/bold magenta]")
        rprint(f"[bold cyan]Regressing confounding factors for {current_tissue}[/bold cyan]")
        sk_resid_mat_norm = regress_confounding (tissue_matrix_preReg , tissue_sample_preReg, current_tissue)
        rprint(f"[bold yellow]Splitting age groups for {current_tissue}[/bold yellow]")
        # split_age_group(sk_resid_mat_norm, sample_attributes, current_tissue)
        split_age_group(sk_resid_mat_norm, tissue_sample_preReg, current_tissue)



