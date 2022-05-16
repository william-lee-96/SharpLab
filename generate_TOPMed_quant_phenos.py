#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 15:26:34 2022

@author: williamlee
"""

import pandas as pd
import numpy as np
import os

wd = '/sc/arion/projects/sharpa01a/William/'
os.chdir(wd)

combined_df = pd.DataFrame()
for file in os.scandir('TOPMed_EHv5_summary_data/'):
    if not file.name.startswith('.'):
                
        contig_df = pd.read_excel(file.path)
        
        combined_df = pd.concat([combined_df, contig_df], ignore_index=True)
        
sample_data_list = []
for sample in combined_df['sample_ID'].unique():

    sample_df = combined_df[combined_df['sample_ID'] == sample]
    
    all_size = sum(sample_df['sum_total']) / sum(sample_df['num_TRs_total'])
    exp_size = sum(sample_df['sum_expansions']) / sum(sample_df['num_expansions'])
    con_size = sum(sample_df['sum_contractions']) / sum(sample_df['num_contractions'])
    
    exp_prop = sum(sample_df['num_expansions']) / sum(sample_df['num_TRs_total'])
    con_prop = sum(sample_df['num_contractions']) / sum(sample_df['num_TRs_total'])
    
    sample_data_list.append([sample, all_size, exp_size, con_size, exp_prop, con_prop, len(sample_df)])

sample_data_df = pd.DataFrame(sample_data_list, columns=['sample_ID', 'mutation_size', 'expansion_size', 'contraction_size', 'expansion_propensity', 'contraction_propensity', 'num_chr'])

sample_data_df.to_excel('TOPMed_quant_phenotypes_TEST.xlsx', index=False)