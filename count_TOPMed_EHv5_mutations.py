#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 12 09:10:19 2022

@author: williamlee
"""

import pandas as pd
import numpy as np
import os

wd = '/sc/arion/projects/sharpa01a/William/'
os.chdir(wd)

os.mkdir('REGENIE_EHv5_pheno')

qc_status_list = ['Phewas_GWAS_FinalFiles', 'GWAS_Will_Unfiltered']
threshold_list = ['VeryLow', 'Low', 'Medium', 'High', 'VeryHigh']
cohort_list = ['ARIC', 'BioMe_Baylor', 'BioMe_MGI', 'CHS', 'COPD_Broad', 'COPD_UW', 'FHS', 'HyperGen', 'JHS', 'MESA', 'MGH', 'VTE', 'VUAF', 'WHI']

for qc_status in qc_status_list:
    
    os.mkdir(f'REGENIE_EHv5_pheno/{qc_status}')
    os.mkdir(f'REGENIE_EHv5_pheno/{qc_status}/mutation_count')

    for threshold in threshold_list:
        
        complete_df = pd.DataFrame()
        
        for cohort in cohort_list:
            
            cohort_df = pd.DataFrame()
            
            for file in os.scandir(f'/sc/arion/projects/sharpa01a/TRE_DataRepository/WGS/hg38/EHv5_GenomewidePolymorphic/{cohort}/{qc_status}/{threshold}'):
                if not file.name.startswith('chrY'):
                    
                    contig_df = pd.read_csv(file.path, sep='\t', usecols=['SampleId', 'Status'])
                    
                    contig_df['Status'].replace('.', 0, inplace=True)
    
                    contig_df['Status'] = pd.to_numeric(contig_df['Status'])
                    
                    contig_df = contig_df.groupby('SampleId', as_index=False).sum()
                    
                    cohort_df = pd.concat([cohort_df, contig_df], ignore_index=True)
                    
            complete_df = pd.concat([complete_df, cohort_df], ignore_index=True)
            
        complete_df = complete_df.groupby('SampleId', as_index=False).sum()
    
        complete_df.rename(columns={'SampleId': 'FID', 'Status': 'pheno'}, inplace=True)
    
        complete_df['IID'] = complete_df['FID']
    
        complete_df = complete_df[['FID', 'IID', 'pheno']]
        
        complete_df.to_csv(f'REGENIE_EHv5_pheno/{qc_status}/mutation_count/REGENIE_pheno_{threshold}.tsv', sep='\t', index=False)
        