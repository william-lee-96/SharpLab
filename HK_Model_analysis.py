#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  8 15:44:38 2021

@author: williamlee
"""

import os
import numpy as np
import pandas as pd

import seaborn as sns
import matplotlib.pyplot as plt

wd = '/Users/williamlee/SharpLab/'
os.chdir(wd)

###############################################################################

me_90_10_loci = pd.read_excel('Data/HK_model/regions_for_Will.xlsx', usecols=[1,2,3,8])
me_90_10_loci.columns = ['chr', 'start', 'end', 'region']

me_90_10_list = []
for region in me_90_10_loci['region'].unique():
    region_df = me_90_10_loci[me_90_10_loci['region'] == region]
    
    if '_unmeth' in region:
        me_90_10_list.append([region_df['chr'].iloc[0], min(region_df['start']), max(region_df['end']), 0.1])
    else:
        me_90_10_list.append([region_df['chr'].iloc[0], min(region_df['start']), max(region_df['end']), 0.9])

###############################################################################

me_50_50_loci = pd.read_csv('Data/HK_model/50-50_methylated_loci.bed', sep='\t', header=None)

me_50_50_loci['pred_me'] = 0.5

me_50_50_list = me_50_50_loci.to_numpy().tolist()

###############################################################################

me_imprinted_loci = pd.read_excel('Data/HK_model/ImprintedDMRs_MethylSeqhg19_PMID27911167.xlsx', usecols=[1,2,3])

me_imprinted_loci['Chr'] = 'chr' + me_imprinted_loci['Chr'].astype('string')

me_imprinted_loci['pred_me'] = 0.5

me_imprinted_list = me_imprinted_loci.to_numpy().tolist()

###############################################################################

me_complete_list = me_90_10_list + me_50_50_list + me_imprinted_list
[l.append([]) for l in me_complete_list]

###############################################################################
###############################################################################
###############################################################################

hk_sample_out = pd.read_csv('out/HK_fin_out.txt', sep='\t')

df1, df2, df3 = np.array_split(hk_sample_out[hk_sample_out['#rname']=='chr1'], 3)
with pd.ExcelWriter('out/HK_out_chr1.xlsx') as writer:  
    df1.to_excel(writer, sheet_name='chr1 - pt.1', index=False)
    df2.to_excel(writer, sheet_name='chr1 - pt.2', index=False)
    df3.to_excel(writer, sheet_name='chr1 - pt.3', index=False)
                        
hk_sample_out['cpg_posn'] = hk_sample_out['#rname'] + '_' + hk_sample_out['start'].astype('string')

hk_merge_df = hk_sample_out[['#rname','start','cpg_posn']]
hk_merge_df.drop_duplicates(subset='cpg_posn', inplace=True, ignore_index=True)
hk_merge_df.set_index('cpg_posn', inplace=True)

hk_temp_df = hk_sample_out[['depthW','depthC','y_pred','cpg_posn']]
hk_temp_df = hk_temp_df.groupby(['cpg_posn']).mean()

hk_complete_df = hk_merge_df.join(hk_temp_df, how='inner')

hk_complete_df['prediction'] = hk_complete_df['y_pred'].apply(lambda x: 1 if x>=0.5 else 0)

hk_complete_list = hk_complete_df.to_numpy().tolist()

for l in me_complete_list:
    [l[4].append(hk[5]) for hk in hk_complete_list if
     hk[0] == l[0] and hk[1] >= l[1] and hk[1] <= l[2]]
    
###############################################################################
###############################################################################
###############################################################################
    
hk_summ_df = pd.DataFrame(me_complete_list, columns=['chr','start','end','expected','observed'])

hk_summ_df['num_calls'] = hk_summ_df['observed'].apply(lambda x: len(x))

hk_summ_df['observed'] = hk_summ_df['observed'].apply(lambda x: np.nan if len(x)==0 else sum(x)/len(x))

hk_summ_df.to_excel('out/HK_out_summary.xlsx', index=False)