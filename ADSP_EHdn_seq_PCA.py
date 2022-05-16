#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 11 11:03:33 2021

@author: williamlee
"""

import os
import pandas as pd
import numpy as np

from sklearn.decomposition import PCA

import seaborn as sns
import matplotlib.pyplot as plt

wd = '/sc/arion/projects/sharpa01a/William/'
os.chdir(wd)

pc_list = ['PC' + str(num) for num in list(range(1, 11))]

num_list = [(0,1),(2,3),(4,5),(6,7),(8,9)]

os.system('mkdir out/ADSP_seq_PCA_plots')
os.system('mkdir out/ADSP_seq_PCA_plots/locus')
os.system('mkdir out/ADSP_seq_PCA_plots/IRR')

mstr_prefix = ['ADSP_hiseqx', 'ADSP_novaseq']

for prefix in mstr_prefix:

    batch_df = pd.read_csv('{}_label.txt'.format(prefix), sep='\t')
    batch_df.set_index('sample', inplace=True)
        
    ###############################################################################
    
    os.system('mkdir out/ADSP_seq_PCA_plots/locus/{}'.format(prefix))
    
    locus_matrix = pd.read_csv('out/{}_EHdn_motif_locus.tsv'.format(prefix), sep='\t')
        
    locus_matrix.replace(np.nan, 0, inplace=True)
    locus_matrix.set_index('motif', inplace=True)
    
    locus_pca = PCA().fit(locus_matrix)
    
    loadings = locus_pca.components_
    
    exp_var_ratio = list(locus_pca.explained_variance_ratio_[0:10])
    exp_var_ratio = [round(x*100, 2) for x in exp_var_ratio]
    
    loadings_df = pd.DataFrame.from_dict(dict(zip(pc_list, loadings)))
    loadings_df['sample'] = locus_matrix.columns.values
    loadings_df.set_index('sample', inplace=True)
    
    for col in batch_df:
        
        if col != 'Sequencer':
        
            os.system('mkdir out/ADSP_seq_PCA_plots/locus/{}/{}'.format(prefix, col))
            
            temp_df = loadings_df.join(batch_df[col], how='inner')
        
            for j,k in num_list:
                plt.figure()
                ax = sns.scatterplot(x=temp_df['{}'.format(pc_list[j])], y=temp_df['{}'.format(pc_list[k])], hue=temp_df[col])
                ax.set_xlabel('{} ({}%)'.format(pc_list[j], exp_var_ratio[j]))
                ax.set_ylabel('{} ({}%)'.format(pc_list[k], exp_var_ratio[k]))
                ax.set_title('{}, {} vs. {}'.format('ADSP locus case-control', pc_list[j], pc_list[k]))
                ax.figure.savefig('out/ADSP_seq_PCA_plots/locus/{}/{}/{}_{}_locus.pdf'.format(prefix, col, pc_list[j], pc_list[k]))
        
    ###############################################################################
    
    os.system('mkdir out/ADSP_seq_PCA_plots/IRR/{}'.format(prefix))
    
    irr_matrix = pd.read_csv('out/{}_EHdn_motif_IRR.tsv'.format(prefix), sep='\t')
        
    irr_matrix.replace(np.nan, 0, inplace=True)
    irr_matrix.set_index('motif', inplace=True)
    
    irr_pca = PCA().fit(irr_matrix)
    
    loadings = irr_pca.components_
    
    exp_var_ratio = list(irr_pca.explained_variance_ratio_[0:10])
    exp_var_ratio = [round(x*100, 2) for x in exp_var_ratio]
    
    loadings_df = pd.DataFrame.from_dict(dict(zip(pc_list, loadings)))
    loadings_df['sample'] = irr_matrix.columns.values
    loadings_df.set_index('sample', inplace=True)
    
    for col in batch_df:
        
        if col != 'Sequencer':
        
            os.system('mkdir out/ADSP_seq_PCA_plots/IRR/{}/{}'.format(prefix, col))
            
            temp_df = loadings_df.join(batch_df[col], how='inner')
        
            for j,k in num_list:
                plt.figure()
                ax = sns.scatterplot(x=temp_df['{}'.format(pc_list[j])], y=temp_df['{}'.format(pc_list[k])], hue=temp_df[col])
                ax.set_xlabel('{} ({}%)'.format(pc_list[j], exp_var_ratio[j]))
                ax.set_ylabel('{} ({}%)'.format(pc_list[k], exp_var_ratio[k]))
                ax.set_title('{}, {} vs. {}'.format('ADSP IRR case-control', pc_list[j], pc_list[k]))
                ax.figure.savefig('out/ADSP_seq_PCA_plots/IRR/{}/{}/{}_{}_irr.pdf'.format(prefix, col, pc_list[j], pc_list[k]))