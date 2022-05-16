#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 13 15:59:46 2022

@author: williamlee
"""

import pandas as pd
import numpy as np
import os

wd = '/sc/arion/projects/sharpa01a/William/'
os.chdir(wd)

os.mkdir('REGENIE_EHv5_covar')

qc_status_list = ['Phewas_GWAS_FinalFiles', 'GWAS_Will_Unfiltered']
cohort_list = ['ARIC', 'BioMe_Baylor', 'BioMe_MGI', 'CHS', 'COPD_Broad', 'COPD_UW', 'FHS', 'HyperGen', 'JHS', 'MESA', 'MGH', 'VTE', 'VUAF', 'WHI']

for qc_status in qc_status_list:
    
    os.mkdir(f'REGENIE_EHv5_covar/{qc_status}')
        
    complete_df = pd.DataFrame()
    
    for cohort in cohort_list:
                
        for file in os.scandir(f'/sc/arion/projects/sharpa01a/TRE_DataRepository/WGS/hg38/EHv5_GenomewidePolymorphic/{cohort}/{qc_status}/VeryLow'):
            if file.name.startswith('chr21'):
                
                contig_df = pd.read_csv(file.path, sep='\t', usecols=['SampleId', 'Cohort', 'InsertSize', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'gender'])
                
                contig_df.drop_duplicates(inplace=True)
                                                             
                complete_df = pd.concat([complete_df, contig_df], ignore_index=True)
        
    complete_df['gender'].replace('male', 1, inplace=True)
    complete_df['gender'].replace('female', 2, inplace=True)
    
    complete_df['FID'] = complete_df['SampleId']
    complete_df.rename(columns={'SampleId': 'IID'}, inplace=True)

    covar_df = complete_df[['FID', 'IID', 'InsertSize', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'gender']]
    
    cohort_df = complete_df[['FID', 'IID', 'Cohort']]
    
    covar_df.to_csv(f'REGENIE_EHv5_covar/{qc_status}/REGENIE_covar.tsv', sep='\t', index=False)
    
    cohort_df.to_csv(f'REGENIE_EHv5_covar/{qc_status}/REGENIE_cohort.tsv', sep='\t', index=False)  