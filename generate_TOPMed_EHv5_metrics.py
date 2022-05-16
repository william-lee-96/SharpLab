#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 30 14:08:22 2022

@author: williamlee
"""

import pandas as pd
import numpy as np
import os
import sys

wd = '/sc/arion/projects/sharpa01a/William/'
os.chdir(wd)

contig = sys.argv[1]

contig_df = pd.DataFrame()

for directory in os.scandir('/sc/arion/scratch/jadhab01/ForParas/'):
    for file in os.scandir(f'/sc/arion/scratch/jadhab01/ForParas/{directory.name}/'):
        if file.name.split('.')[0] == contig:
            
            cohort_df = pd.read_csv(file.path, sep='\t')
            
            contig_df = pd.concat([contig_df, cohort_df], ignore_index=True)
            
contig_df['tandem_repeat'] = contig_df['VARID'] + '_' + contig_df['RefUnit']

contig_df['AvgCopy'] = (contig_df['A1'] + contig_df['A2']) / 2

copy_dict = {}
for repeat in contig_df['tandem_repeat'].unique():
    
    repeat_df = contig_df[contig_df['tandem_repeat'] == repeat]
    
    copy_dict[repeat] = repeat_df['AvgCopy'].median()

copy_df = pd.DataFrame.from_dict(copy_dict, orient='index')

copy_df.reset_index(inplace=True)

copy_df.columns = ['tandem_repeat', 'PopCopy']

contig_df = contig_df.merge(copy_df, on='tandem_repeat', how='left')

contig_df['diff'] = contig_df['AvgCopy'] - contig_df['PopCopy']

contig_df = contig_df[contig_df['diff'].notna()]

sample_data_list = []
for sample in contig_df['SampleId'].unique():

    sample_df = contig_df[contig_df['SampleId'] == sample]
    
    diff_list = sample_df['diff'].tolist()
        
    exp_list = [x for x in diff_list if x > 0]
        
    con_list = [x for x in diff_list if x < 0]
        
    sample_data_list.append([sample, sum(diff_list), sum(exp_list), sum(con_list), len(diff_list), len(exp_list), len(con_list)])
    
sample_data_df = pd.DataFrame(sample_data_list, columns=['sample_ID', 'sum_total', 'sum_expansions', 'sum_contractions', 'num_TRs_total', 'num_expansions', 'num_contractions'])

sample_data_df.to_excel(f'TOPMed_EHv5_summary_data/{contig}_TOPMed_EHv5.xlsx', index=False)