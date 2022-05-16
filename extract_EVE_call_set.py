#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  2 17:32:38 2022

@author: williamlee
"""

import pandas as pd
import numpy as np
import os

wd = '/sc/arion/projects/sharpa01a/William/'
os.chdir(wd)

###############################################################################

eve_combined_df = pd.DataFrame(columns=['EVE','merge_col'])

for file in os.scandir('EVE_70_vcfs/'):
    if file.name.endswith('.vcf'):

        eve_df = pd.read_csv(file.path, sep='\t', skiprows=20)
        
        pos1 = eve_df['REF'].str[0] != eve_df['ALT'].str[0]
        pos2 = eve_df['REF'].str[1] != eve_df['ALT'].str[1]
        pos3 = eve_df['REF'].str[2] != eve_df['ALT'].str[2]
        
        snv_idx = pd.concat([pos1,pos2,pos3], axis=1).sum(axis=1)
        snv_idx = snv_idx[snv_idx == 1].index
        
        eve_df = eve_df.iloc[snv_idx,:]
        
        eve_df.loc[eve_df['REF'].str[0] != eve_df['ALT'].str[0], 'bias'] = 0
        eve_df.loc[eve_df['REF'].str[1] != eve_df['ALT'].str[1], 'bias'] = 1
        eve_df.loc[eve_df['REF'].str[2] != eve_df['ALT'].str[2], 'bias'] = 2
        
        eve_df.loc[eve_df['bias'] == 0, 'REF'] = eve_df['REF'].str[0]
        eve_df.loc[eve_df['bias'] == 0, 'ALT'] = eve_df['ALT'].str[0]
        eve_df.loc[eve_df['bias'] == 1, 'REF'] = eve_df['REF'].str[1]
        eve_df.loc[eve_df['bias'] == 1, 'ALT'] = eve_df['ALT'].str[1]
        eve_df.loc[eve_df['bias'] == 2, 'REF'] = eve_df['REF'].str[2]
        eve_df.loc[eve_df['bias'] == 2, 'ALT'] = eve_df['ALT'].str[2]
        
        eve_df['POS'] = eve_df['POS'] + eve_df['bias']
        eve_df['POS'] = eve_df['POS'].astype(int)
        
        eve_df['EVE'] = eve_df['INFO'].str.split(pat=';').apply(lambda x: x[12].split('=')[1])
        
        eve_df = eve_df[eve_df['EVE'] != 'Uncertain']
        
        eve_df['merge_col'] = eve_df['#CHROM'].astype('string') + ':' + eve_df['POS'].astype('string') + '_' +  eve_df['ALT'] + '_' + eve_df['REF']
        
        eve_df = eve_df[['EVE','merge_col']]
        
        eve_combined_df = pd.concat([eve_combined_df, eve_df], ignore_index=True)
        
eve_combined_df.to_csv('EVE_DVA_call_set.tsv', sep='\t', index=False)