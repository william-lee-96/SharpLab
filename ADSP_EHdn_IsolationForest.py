#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 23 12:28:48 2021

@author: williamlee
"""

import pandas as pd
import numpy as np
import json
import os

from sklearn.ensemble import IsolationForest

wd = '/sc/arion/projects/sharpa01a/William/'
os.chdir(wd)

os.system('mkdir ADSP_EHdn_IF')

human_chr = ('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22')

batch_list = ['hiseqx_batch_1', 'hiseqx_batch_2', 'hiseqx_batch_3', 'hiseqx_batch_4', 'hiseqx_batch_5', 'hiseqx_batch_6', 'hiseqx_batch_7', 'hiseqx_batch_8', 'novaseq']

def split_locus(locus):
    str_list_1 = locus.split(':')
    contig = str_list_1[0]
    start_end = str_list_1[1]
    str_list_2 = start_end.split('-')
    start = str_list_2[0]
    end = str_list_2[1]
    return [contig, start, end]

for batch in batch_list:
    
    with open('{}.multisample_profile.json'.format(batch)) as f:
        EHdn_json_dict = json.load(f)
    
    mstr_counts_df = pd.DataFrame.from_dict(EHdn_json_dict['Counts'], orient='index')
    mstr_counts_df.reset_index(inplace=True)
    mstr_counts_df.columns.values[0] = 'motif'
    mstr_counts_df = mstr_counts_df[mstr_counts_df['RegionsWithIrrAnchors'].notna()]
    mstr_counts_df.reset_index(drop=True, inplace=True)
    
    mstr_parameters_df = pd.DataFrame.from_dict(EHdn_json_dict['Parameters'])
    mstr_parameters_df.index.names = ['sample_id']
    mstr_parameters_df['norm_factor'] = 40/(mstr_parameters_df['Depths'])
    
    label_df = pd.read_csv('{}_manifest.txt'.format(batch), sep='\t', header=None, usecols=[0,1])
    label_df.columns = ['sample_id', 'status']
    label_df.set_index('sample_id', inplace=True)
    
    total_case = len(label_df[label_df['status']=='case'])
    total_control = len(label_df[label_df['status']=='control'])
    
    results_list = []
    for row in range(len(mstr_counts_df)):
        
        motif = mstr_counts_df.loc[row,'motif']
        region_dict = mstr_counts_df.loc[row,'RegionsWithIrrAnchors']
        
        for region, sample_dict in region_dict.items():
            
            if len(sample_dict) < 5:
                continue # skip to next loop if no. samples < 5
        
            locus = split_locus(region)
            
            if locus[0] not in human_chr: 
                continue # skip to next loop if chr is not chr1-22
                                            
            sample_counts = pd.DataFrame(zip(sample_dict.keys(), sample_dict.values()), columns = ['sample_id', 'counts'])
            sample_counts.set_index('sample_id', inplace=True)
            
            sample_counts_merge = sample_counts.join(label_df, how='left')
            
            if sum(sample_counts_merge['status']=='case') == 0:
                continue # skip to next loop if there are no cases
            
            sample_counts_merge = sample_counts_merge.join(mstr_parameters_df, how='left')
            sample_counts_merge['norm_counts'] = sample_counts_merge['counts']*sample_counts_merge['norm_factor']
            
            sample_counts_merge.reset_index(inplace=True)
            
            avg_norm_count = sample_counts_merge['norm_counts'].mean()
    
            iso = IsolationForest(random_state=0)
            
            X = sample_counts_merge['norm_counts'].to_numpy().reshape(-1, 1)
            
            sample_counts_merge['iso_lab'] = iso.fit_predict(X)
            
            outlier_df = sample_counts_merge[(sample_counts_merge['iso_lab']==-1) & (sample_counts_merge['norm_counts']>avg_norm_count)]
            
            if sum(outlier_df['status']=='case') == 0:
                continue # skip to next loop if there are no cases
            
            avg_out_norm_count = outlier_df['norm_counts'].mean()
            
            num_case_out = sum(outlier_df['status']=='case')
                                            
            prop_case = num_case_out/len(outlier_df)
                            
            results_list.append([locus[0], locus[1], locus[2], motif, num_case_out, prop_case, avg_norm_count, avg_out_norm_count])
            
    results_df = pd.DataFrame(results_list, columns=['chr', 'start', 'end', 'motif', 'N_case_outliers', 'prop_cases', 'avg_norm_count', 'avg_outlier_norm_count'])
    
    results_df.sort_values(by=['prop_cases','N_case_outliers'], ascending=[False,False], inplace=True, ignore_index=True)
        
    results_df.insert(3, 'ref', pd.Series(['0']*len(results_df)))
    results_df.insert(4, 'alt', pd.Series(['0']*len(results_df)))
    
    results_df.iloc[:,0:5].to_csv('out/ADSP_EHdn_unannotated.txt', sep='\t', index=False)
    
    results_df.drop(columns=['ref','alt'], inplace=True) 

    os.system('annovar/annotate_variation.pl -out out/ADSP_EHdn_annotated -build hg38 out/ADSP_EHdn_unannotated.txt humandb/')
    
    annotations = pd.read_csv('out/ADSP_EHdn_annotated.variant_function', sep='\t', header=None, usecols=[0,1], names=['region','gene'])
        
    results_df.insert(3, 'gene', annotations['gene'])
    results_df.insert(4, 'region', annotations['region'])
    
    results_df.to_excel('ADSP_EHdn_IF/{}_EHdn_IF.xlsx'.format(batch), index=False)