#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  4 11:19:03 2022

@author: williamlee
"""

import pandas as pd
import numpy as np
import os

wd = '/sc/arion/projects/sharpa01a/William/'
os.chdir(wd)

eve_complete_df = pd.read_csv('EVE_DVA_call_set.tsv', sep='\t')

cadd_complete_df = pd.DataFrame(columns=['PHRED','merge_col'])
for directory in os.scandir('PLINK_to_VCF_out/'):
    for file in os.scandir(f'PLINK_to_VCF_out/{directory.name}/'):
        if file.name.endswith('.tsv.gz'):
            
            cadd_df = pd.read_csv(file.path, sep='\t', skiprows=1)
            cadd_df['merge_col'] = cadd_df['#Chrom'].astype('string') + ':' + cadd_df['Pos'].astype('string') + '_' +  cadd_df['Alt'] + '_' + cadd_df['Ref']
            cadd_df = cadd_df[['PHRED','merge_col']]
            
            cadd_complete_df = pd.concat([cadd_complete_df, cadd_df], ignore_index=True)

cadd_complete_df.drop_duplicates(inplace=True, ignore_index=True)

revel_complete_df = pd.DataFrame(columns=['REVEL','merge_col'])
for file in os.scandir('REVEL_by_chr/'):

    revel_df = pd.read_csv(file.path, sep='\t')
    revel_complete_df = pd.concat([revel_complete_df, revel_df], ignore_index=True)
    
revel_complete_df.drop_duplicates(inplace=True, ignore_index=True)

skato_var_df = pd.read_excel('SKATO_top5_for_scoring.xlsx')

skato_var_df = skato_var_df.merge(eve_complete_df, on='merge_col', how='left')
skato_var_df = skato_var_df.merge(cadd_complete_df, on='merge_col', how='left')
skato_var_df = skato_var_df.merge(revel_complete_df, on='merge_col', how='left')

skato_var_df.to_excel('SKATO_variants_scored.xlsx', index=False)