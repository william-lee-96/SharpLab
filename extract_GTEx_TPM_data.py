#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 10 14:26:09 2023

@author: williamlee
"""

import os
import pandas as pd
import numpy as np

os.chdir('/Users/williamlee/SharpLab')

gene_df = pd.read_excel('Data/GTEx/ensembl_ID_to_gene_name.xlsx', header=None)
gene_df.columns = ['GENE_NAME', 'GENE_ID']
gene_name_list = gene_df['GENE_NAME'].to_list()

tpm_df = pd.DataFrame()
for chunk in pd.read_csv('Data/GTEx/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz', sep='\t', skiprows=2, chunksize=1000):
    
    chunk = chunk[chunk['Description'].isin(gene_name_list)]
    
    tpm_df = pd.concat([tpm_df, chunk], ignore_index=True)
    
tpm_df = tpm_df.T

tpm_df.rename(columns=tpm_df.iloc[1], inplace=True)

tpm_df = tpm_df.iloc[2:,:]

tpm_df.reset_index(inplace=True)

tpm_df.rename(columns={'index':'sample_ID'}, inplace=True)

tpm_df.to_excel('out/GTEx_TPM_25_genes.xlsx', index=False)

###############################################################################

import os
import pandas as pd
import numpy as np

os.chdir('/Users/williamlee/SharpLab')

tpm_df = pd.read_excel('out/GTEx_TPM_25_genes.xlsx')

tpm_df['subject_ID'] = tpm_df['sample_ID'].apply(lambda x: '-'.join(x.split('-', 2)[:2]))

meta_df = pd.read_csv('Data/GTEx/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt', sep='\t', usecols=['SAMPID', 'SMTSD'])
meta_df.columns = ['sample_ID', 'tissue']

tpm_df = tpm_df.merge(meta_df, on='sample_ID', how='left')

tpm_df = tpm_df[['subject_ID', 'sample_ID', 'tissue', 'ADARB1', 'C11orf80', 'CSNK1E', 'FRA10AC1', 'LINGO3', 'PCMTD2', 'ZNF713']]

case_dict = {'ADARB1': ['GTEX-1CB4H'],
             'C11orf80': ['GTEX-1F48J'],
             'CSNK1E': ['GTEX-X62O'],
             'FRA10AC1': ['GTEX-15EO6', 'GTEX-16AAH', 'GTEX-14JG1', 'GTEX-ZP4G'],
             'LINGO3': ['GTEX-1NHNU', 'GTEX-14PQA'],
             'PCMTD2': ['GTEX-13D11'],
             'ZNF713': ['GTEX-18D9U', 'GTEX-1EN7A', 'GTEX-1F75B', 'GTEX-1JK1U', 'GTEX-13NZA', 'GTEX-WH7G']}

for col in tpm_df.columns[3:]:
        
    case_list = case_dict[col]
    
    count_df = []
    for tissue in tpm_df['tissue'].unique():
        count_df.append([tissue, sum(tpm_df[tpm_df['tissue'] == tissue]['subject_ID'].isin(case_list))])
        
    count_df = pd.DataFrame(count_df, columns =['tissue', 'case_count'])
        
    mean_df = tpm_df[['tissue', col]].groupby(by='tissue', sort=False, as_index=False).mean()
    mean_df.columns = ['tissue', 'mean_TPM']
    
    sort_df = count_df.merge(mean_df, on='tissue')
    
    sort_df.sort_values(by=['case_count', 'mean_TPM'], ascending=[False, False], inplace=True)
    
    sort_df.to_excel(f'out/GTEx_TPM_metadata_tables/{col}_TPM_metadata.xlsx', index=False)
        
    out_df = tpm_df[['subject_ID', 'sample_ID', 'tissue', col]]
    out_df.rename(columns={col:'TPM'}, inplace=True)
    out_df = out_df[out_df['tissue'] == sort_df.iloc[0,0]]
    
    out_df['status'] = np.nan
    out_df.loc[out_df.subject_ID.isin(case_list), 'status'] = 'case'
    out_df.loc[~out_df.subject_ID.isin(case_list), 'status'] = 'control'
    out_df.sort_values(by='TPM', ascending=False, inplace=True)
    out_df = out_df[['subject_ID', 'sample_ID', 'tissue', 'status', 'TPM']]

    out_df.to_excel(f'out/GTEx_TPM_epimutation_tables/{col}_TPM_epimutation.xlsx', index=False)
    
###############################################################################

import os
import pandas as pd
import numpy as np

os.chdir('/Users/williamlee/SharpLab')

tpm_df = pd.read_excel('out/GTEx_TPM_25_genes.xlsx')

tpm_df['subject_ID'] = tpm_df['sample_ID'].apply(lambda x: '-'.join(x.split('-', 2)[:2]))

meta_df = pd.read_csv('Data/GTEx/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt', sep='\t', usecols=['SAMPID', 'SMTSD'])
meta_df.columns = ['sample_ID', 'tissue']

tpm_df = tpm_df.merge(meta_df, on='sample_ID', how='left')

tpm_df = tpm_df[['subject_ID', 'sample_ID', 'tissue', 'ADARB1', 'C11orf80', 'CSNK1E', 'FRA10AC1', 'LINGO3', 'PCMTD2', 'ZNF713']]

case_dict = {'ADARB1': ['GTEX-1CB4H'],
             'C11orf80': ['GTEX-1F48J'],
             'CSNK1E': ['GTEX-X62O'],
             'FRA10AC1': ['GTEX-15EO6', 'GTEX-16AAH', 'GTEX-14JG1', 'GTEX-ZP4G'],
             'LINGO3': ['GTEX-1NHNU', 'GTEX-14PQA'],
             'PCMTD2': ['GTEX-13D11'],
             'ZNF713': ['GTEX-18D9U', 'GTEX-1EN7A', 'GTEX-1F75B', 'GTEX-1JK1U', 'GTEX-13NZA', 'GTEX-WH7G']}

summary_df = []
for col in tpm_df.columns[3:]:
        
    case_list = case_dict[col]
    
    for tissue in tpm_df['tissue'].unique():
        
        if sum(tpm_df[tpm_df['tissue'] == tissue]['subject_ID'].isin(case_list)) > 0:
        
            tissue_df = tpm_df[tpm_df['tissue'] == tissue]
            tissue_df = tissue_df[['subject_ID', 'sample_ID', 'tissue', col]]
            
            tissue_df['status'] = np.nan
            tissue_df.loc[tissue_df.subject_ID.isin(case_list), 'status'] = 'case'
            
            tissue_df.sort_values(by=col, ascending=False, inplace=True)
            
            tissue_df['rank'] = range(1, len(tissue_df)+1)
            
            case_df = tissue_df[tissue_df['status'] == 'case']
                        
            case_df['N'] = len(tissue_df)
            
            case_df['percentile'] = 1-(case_df['rank']/case_df['N'])
            case_df['percentile'] = case_df['percentile'].round(3)
            
            case_df['gene'] = col
            
            case_df = case_df[['subject_ID', 'tissue', 'sample_ID', 'gene', 'N', 'rank', 'percentile']]
            
            summary_df.append(case_df)
            
summary_df = pd.concat(summary_df, ignore_index=True)

summary_df.sort_values(by=['percentile', 'N'], ascending=[True, False], inplace=True)
            
summary_df.to_excel('out/GTEx_TPM_epimutation_summary_table.xlsx', index=False)