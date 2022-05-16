#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 10 16:47:59 2021

@author: williamlee
"""

import pandas as pd 
import json
import os

wd = '/sc/arion/projects/sharpa01a/William/'
os.chdir(wd)

human_chr = ('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12',
             'chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22')

def split_locus(locus):
    str_list_1 = locus.split(':')
    contig = str_list_1[0]
    start_end = str_list_1[1]
    str_list_2 = start_end.split('-')
    start = str_list_2[0]
    end = str_list_2[1]
    return [contig, start, end]

mstr_prefix = ['ADSP_hiseqx', 'ADSP_novaseq']

for prefix in mstr_prefix:

    with open('{}.multisample_profile.json'.format(prefix)) as f:
        EHdn_json_dict = json.load(f)
    
    mstr_counts_df = pd.DataFrame.from_dict(EHdn_json_dict['Counts'], orient='index')
    mstr_counts_df.reset_index(inplace=True)
    mstr_counts_df.columns.values[0] = 'motif'
    mstr_counts_df = mstr_counts_df[mstr_counts_df['RegionsWithIrrAnchors'].notna()]
    mstr_counts_df.reset_index(drop=True, inplace=True)
    
    mstr_parameters_df = pd.DataFrame.from_dict(EHdn_json_dict['Parameters'])
    mstr_parameters_df.index.names = ['sample_id']
    mstr_parameters_df['norm_factor'] = 40/(mstr_parameters_df['Depths'])
    
    label_df = pd.read_csv('{}_manifest.txt'.format(prefix), sep='\t', header=None, usecols=[0,1])
    label_df.columns = ['sample_id', 'status']
    label_df.set_index('sample_id', inplace=True)
        
    supplementary_df = pd.DataFrame(index=label_df.index, columns=['num_uniq_exp','sum_norm_IRR','motif:num_uniq_exp','motif:sum_norm_IRR'])
    supplementary_df['num_uniq_exp'] = 0
    supplementary_df['sum_norm_IRR'] = 0
    supplementary_df['motif:num_uniq_exp'] = supplementary_df['motif:num_uniq_exp'].apply(lambda x: {})
    supplementary_df['motif:sum_norm_IRR'] = supplementary_df['motif:sum_norm_IRR'].apply(lambda x: {})
    
    motif_list = []
    for row in range(len(mstr_counts_df)):
        
        motif = mstr_counts_df.loc[row,'motif']
        region_dict = mstr_counts_df.loc[row,'RegionsWithIrrAnchors']
        
        for region, sample_dict in region_dict.items():
            
            locus = split_locus(region)
            
            if locus[0] not in human_chr: 
                continue # skip to next loop if contig is not chr1-22
            
            else:
                
                motif_list.append(motif) #*#
                                
                sample_counts = pd.DataFrame(zip(sample_dict.keys(), sample_dict.values()), columns = ['sample_id', 'counts'])
                sample_counts.set_index('sample_id', inplace=True)
                
                sample_counts_merge = sample_counts.join(label_df, how='left')
                sample_counts_merge = sample_counts_merge.join(mstr_parameters_df, how='left')
                sample_counts_merge['norm_counts'] = sample_counts_merge['counts']*sample_counts_merge['norm_factor']
                            
                for idx in sample_counts_merge.index:
                    
                    supplementary_df.loc[idx,'num_uniq_exp'] += 1
                    supplementary_df.loc[idx,'sum_norm_IRR'] += sample_counts_merge.loc[idx,'norm_counts']
                    
                    if motif in supplementary_df.loc[idx,'motif:num_uniq_exp']:
                        supplementary_df.loc[idx,'motif:num_uniq_exp'][motif] += 1
                    else:
                        supplementary_df.loc[idx,'motif:num_uniq_exp'][motif] = 1
                        
                    if motif in supplementary_df.loc[idx,'motif:sum_norm_IRR']:
                        supplementary_df.loc[idx,'motif:sum_norm_IRR'][motif] += sample_counts_merge.loc[idx,'norm_counts']
                    else:
                        supplementary_df.loc[idx,'motif:sum_norm_IRR'][motif] = sample_counts_merge.loc[idx,'norm_counts']
                        
    supplementary_df = label_df.join(supplementary_df, how='inner')
    
    supplementary_df.sort_values(by='num_uniq_exp', ascending=False, inplace=True)
    
    supplementary_df.reset_index(inplace=True) 
    
    ###############################################################################
    
    motif_uniq = sorted(set(motif_list), key=len)
    
    ###############################################################################
    
    motif_exp_df = pd.DataFrame(index=motif_uniq, columns=supplementary_df['sample_id'].tolist())
    
    j = 0
    for samp_dict in supplementary_df['motif:num_uniq_exp']:
        
        samp_id = supplementary_df['sample_id'][j]
        for key in samp_dict:
            motif_exp_df.loc[key, samp_id] = samp_dict[key]
        
        j += 1
      
    motif_exp_df.reset_index(inplace=True)
    motif_exp_df.columns.values[0] = 'motif'
    
    ###############################################################################
    
    motif_IRR_df = pd.DataFrame(index=motif_uniq, columns=supplementary_df['sample_id'].tolist())
    
    j = 0
    for samp_dict in supplementary_df['motif:sum_norm_IRR']:
        
        samp_id = supplementary_df['sample_id'][j]
        for key in samp_dict:
            motif_IRR_df.loc[key, samp_id] = samp_dict[key]
        
        j += 1
    
    motif_IRR_df.reset_index(inplace=True)
    motif_IRR_df.columns.values[0] = 'motif'
    
    ###############################################################################
    
    supplementary_df.drop(columns=['motif:num_uniq_exp','motif:sum_norm_IRR'], inplace=True)
    
    ###############################################################################
    
    supplementary_df.to_csv('out/{}_EHdn_supp_main.tsv'.format(prefix), sep='\t', index=False)
    motif_exp_df.to_csv('out/{}_EHdn_motif_locus.tsv'.format(prefix), sep='\t', index=False)
    motif_IRR_df.to_csv('out/{}_EHdn_motif_IRR.tsv'.format(prefix), sep='\t', index=False)