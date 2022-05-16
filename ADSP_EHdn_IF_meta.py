#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 24 05:39:21 2021

@author: williamlee
"""

import pandas as pd
import numpy as np
import os

wd = '/sc/arion/projects/sharpa01a/William/'
os.chdir(wd)

###############################################################################

case_batches = ['hiseqx_batch_4_EHdn_IF.xlsx', 'hiseqx_batch_7_EHdn_IF.xlsx']

complete_df = pd.DataFrame()
i = 1
for file in os.scandir('ADSP_EHdn_IF'):
    
    if file.name not in case_batches:
    
        batch_df = pd.read_excel(file.path)
        
        batch_df['batch'] = i
        
        complete_df = pd.concat([complete_df, batch_df])
        
        i += 1
        
complete_df.drop(columns=['gene','region','avg_norm_count','avg_outlier_norm_count'], inplace=True) 
    
complete_df.reset_index(drop=True, inplace=True)

###############################################################################

def detect_overlap(row, df):
    
    df['overlap'] = np.nan
    for i in range(len(df)):
        if (row['start'] > df.iloc[i]['end'] or row['end'] < df.iloc[i]['start']):
            df.iloc[i,7] = 0
        else:
            df.iloc[i,7] = 1   
            
    return df

###############################################################################
    
results_list = []
while not complete_df.empty:
    
    row = complete_df.iloc[0]
    
    contig = row['chr']
    motif = row['motif']
    
    temp_df = complete_df[(complete_df['chr']==contig) & (complete_df['motif']==motif)]
    
    temp_df = detect_overlap(row, temp_df)
    
    temp_df = temp_df[temp_df['overlap']==1]
    
    if len(temp_df) == 1:
        
        row.drop(labels='batch', inplace=True)
        row['N_batches'] = 1
                
        results_list.append(row.to_list())
        
        complete_df.drop(index=row.name, inplace=True)

    else:
        
        temp_df['N_outliers'] = temp_df['N_case_outliers']/temp_df['prop_cases']
        
        start = temp_df['start'].min()
        end = temp_df['end'].max()
        num_case_out = temp_df['N_case_outliers'].sum()                                          
        prop_case = num_case_out/temp_df['N_outliers'].sum()  
        num_batch = len(temp_df['batch'].unique())
                        
        results_list.append([contig, start, end, motif, num_case_out, prop_case, num_batch])
        
        complete_df.drop(index=temp_df.index, inplace=True)

results_df = pd.DataFrame(results_list, columns=['chr', 'start', 'end', 'motif', 'N_case_outliers', 'prop_cases', 'N_batches'])
    
results_df.sort_values(by=['prop_cases','N_case_outliers'], ascending=[False,False], inplace=True, ignore_index=True)

results_df.insert(3, 'ref', pd.Series(['0']*len(results_df)))
results_df.insert(4, 'alt', pd.Series(['0']*len(results_df)))

results_df.iloc[:,0:5].to_csv('out/ADSP_EHdn_unannotated.txt', sep='\t', index=False)

results_df.drop(columns=['ref','alt'], inplace=True) 

os.system('annovar/annotate_variation.pl -out out/ADSP_EHdn_annotated -build hg38 out/ADSP_EHdn_unannotated.txt humandb/')

annotations = pd.read_csv('out/ADSP_EHdn_annotated.variant_function', sep='\t', header=None, usecols=[0,1], names=['region','gene'])
    
results_df.insert(3, 'gene', annotations['gene'])
results_df.insert(4, 'region', annotations['region'])

results_df.to_excel('ADSP_EHdn_IF_meta.xlsx', index=False)