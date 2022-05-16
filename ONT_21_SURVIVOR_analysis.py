#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  1 13:00:32 2021

@author: williamlee
"""

import pandas as pd
import numpy as np
import os

wd = '/Users/williamlee/SharpLab/'
os.chdir(wd)

survivor_df = pd.read_csv('Data/ONT/ONT_21_out_merged.vcf', sep='\t', skiprows=3399)

survivor_df.columns.values[9:30] = [x.split('.')[0] for x in survivor_df.columns.values[9:30].tolist()]

human_chr = ('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22')

survivor_df = survivor_df[survivor_df['#CHROM'].isin(human_chr)]

survivor_df['temp'] = survivor_df['INFO'].str.split(pat=';')

survivor_df['end_chr'] = survivor_df['temp'].apply(lambda x: x[5].split('=')[1])

survivor_df = survivor_df[survivor_df['end_chr'].isin(human_chr)]

survivor_df.reset_index(drop=True, inplace=True)

survivor_df['N_supp_samp'] = survivor_df['temp'].apply(lambda x: x[0].split('=')[1])
survivor_df['SV_len'] = survivor_df['temp'].apply(lambda x: x[2].split('=')[1])
survivor_df['SV_type'] = survivor_df['temp'].apply(lambda x: x[3].split('=')[1])
survivor_df['end'] = survivor_df['temp'].apply(lambda x: x[6].split('=')[1])
survivor_df['CI_pos'] = survivor_df['temp'].apply(lambda x: x[7].split('=')[1])
survivor_df['CI_end'] = survivor_df['temp'].apply(lambda x: x[8].split('=')[1])
survivor_df['strands'] = survivor_df['temp'].apply(lambda x: x[9].split('=')[1])

supp_samp = []
for row in range(len(survivor_df)):
    
    samp_sv = survivor_df.iloc[row][9:30]
    samp_sv = samp_sv[samp_sv != './.:NaN:0:0,0:--:NaN:NaN:NaN:NAN:NAN:NAN']
    
    supp_samp.append(samp_sv.index.tolist())
   
survivor_df['supp_samp'] = supp_samp

survivor_df = survivor_df[['#CHROM', 'POS', 'REF', 'ALT', 'end_chr', 'N_supp_samp', 'SV_len', 'SV_type', 'end', 'CI_pos', 'CI_end', 'strands', 'supp_samp']]

survivor_df.columns.values[0:4] = ['chr', 'pos', 'ref_bases', 'alt_bases']

survivor_df = survivor_df[['chr', 'pos', 'end', 'end_chr', 'SV_type', 'SV_len', 'ref_bases', 'alt_bases', 'CI_pos', 'CI_end', 'strands', 'supp_samp', 'N_supp_samp']]

survivor_df = survivor_df.apply(pd.to_numeric, errors='ignore')

survivor_df.sort_values(by=['N_supp_samp','SV_len'], ascending=[False,False], inplace=True, ignore_index=True)

insertion_df = survivor_df[survivor_df['SV_type']=='INS']
deletion_df = survivor_df[survivor_df['SV_type']=='DEL']
translocation_df = survivor_df[survivor_df['SV_type']=='TRA']
duplication_df = survivor_df[survivor_df['SV_type']=='DUP']
inversion_df = survivor_df[survivor_df['SV_type']=='INV']

with pd.ExcelWriter('out/ONT_21_LR_SV_merge.xlsx') as writer:  
    insertion_df.to_excel(writer, sheet_name='insertion', index=False)
    deletion_df.to_excel(writer, sheet_name='deletion', index=False)
    translocation_df.to_excel(writer, sheet_name='translocation', index=False)
    duplication_df.to_excel(writer, sheet_name='duplication', index=False)
    inversion_df.to_excel(writer, sheet_name='inversion', index=False)