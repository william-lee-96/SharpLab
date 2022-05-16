#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  3 10:29:45 2022

@author: williamlee
"""

import pandas as pd
import numpy as np
import os

wd = '/sc/arion/projects/sharpa01a/William/'
os.chdir(wd)

os.system('mkdir DVA_PLINKs_REGENIE')

eve_df = pd.read_csv('EVE_DVA_call_set.tsv', sep='\t')

for directory in os.scandir('DVA_PLINKs_cleaned/'):
        os.system(f'mkdir DVA_PLINKs_REGENIE/{directory.name}')
        for file in os.scandir(f'DVA_PLINKs_cleaned/{directory.name}/'):
            if file.name.endswith('.bim'):
                
                prefix = file.name.replace('.bim','')
                contig = prefix.split('_')[0]
                
                bim_df = pd.read_csv(file.path, sep='\t', header=None, usecols=[1])
                bim_df.columns = ['merge_col']
                
                cadd_df = pd.read_csv(f'PLINK_to_VCF_out/{directory.name}/{prefix}.tsv.gz', sep='\t', skiprows=1)
                cadd_df['merge_col'] = cadd_df['#Chrom'].astype('string') + ':' + cadd_df['Pos'].astype('string') + '_' +  cadd_df['Alt'] + '_' + cadd_df['Ref']
                cadd_df = cadd_df[['PHRED','merge_col']]
                
                revel_df = pd.read_csv(f'REVEL_by_chr/{contig}_REVEL.tsv', sep='\t')
                
                merge_df = bim_df.merge(eve_df, on='merge_col', how='left')
                merge_df = merge_df.merge(cadd_df, on='merge_col', how='left')
                merge_df = merge_df.merge(revel_df, on='merge_col', how='left')
                
                merge_df = merge_df[merge_df['EVE'] != 'Benign']
    
                merge_df = merge_df[(merge_df['EVE'] == 'Pathogenic') |
                                    ((merge_df['PHRED'] > 25) |
                                     (merge_df['REVEL'] > 0.5))]
                
                merge_df['merge_col'].to_csv('DVA_SNV_call_set.tsv', sep='\t', header=False)
                
                os.system(f'plink --bfile DVA_PLINKs_cleaned/{directory.name}/{prefix} --extract DVA_SNV_call_set.tsv --make-bed --out DVA_PLINKs_REGENIE/{directory.name}/{prefix}')