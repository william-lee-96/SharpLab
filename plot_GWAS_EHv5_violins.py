#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 28 12:14:14 2022

@author: williamlee
"""

import pandas as pd
import numpy as np
import os

import seaborn as sns
import matplotlib.pyplot as plt

wd = '/sc/arion/projects/sharpa01a/William/'
os.chdir(wd)

os.mkdir('GWAS_EHv5_violin_plots')

###############################################################################

pheno_df = pd.read_csv('REGENIE_EHv5_pheno/Phewas_GWAS_FinalFiles/mutation_count/REGENIE_pheno_Medium.tsv', sep='\t', usecols=['IID', 'pheno'])

cohort_df = pd.read_csv('REGENIE_EHv5_covar/Phewas_GWAS_FinalFiles/REGENIE_cohort.tsv', sep='\t', usecols=['IID', 'Cohort'])

merge_df = pheno_df.merge(cohort_df, on='IID')

###############################################################################

meta_df = pd.read_csv('GWAS_EHv5_Medium_meta1.txt', sep='\t')

meta_df.sort_values('P-value', inplace=True)

meta_df = meta_df.iloc[:50,:]

for snp in meta_df['MarkerName']:
    
    pval = meta_df[meta_df['MarkerName'] == snp].iloc[0,5]
    
    snp_fn = snp.replace(':', '_')
    
    contig = snp.split(':')[0]
    
    genotype_df = pd.DataFrame(columns=['IID', 'genotype'])
    
    ###############################################################################
    
    if os.system(f'plink --bfile GWAS_PLINKs/EUR/chr{contig}_merged --list --snp {snp} --out plink_EUR') == 0:
    
        with open('plink_EUR.list') as file:
            data = [line.rstrip() for line in file]
        
        data = data[:3]
        
        eur_df = pd.DataFrame(columns=['IID', 'genotype'])
        for line in data:
            
            line_list = line.split(' ')
            
            temp_df = pd.DataFrame(list(set(line_list[3:])), columns=['IID'])
            
            temp_df['genotype'] = line_list[2]
            
            eur_df = pd.concat([eur_df, temp_df], ignore_index=True)
            
        genotype_df = pd.concat([genotype_df, eur_df], ignore_index=True)
                    
    ###############################################################################
    
    if os.system(f'plink --bfile GWAS_PLINKs/AFR/chr{contig}_merged --list --snp {snp} --out plink_AFR') == 0:

        with open('plink_AFR.list') as file:
            data = [line.rstrip() for line in file]
        
        data = data[:3]
        
        afr_df = pd.DataFrame(columns=['IID', 'genotype'])
        for line in data:
            
            line_list = line.split(' ')
            
            temp_df = pd.DataFrame(list(set(line_list[3:])), columns=['IID'])
            
            temp_df['genotype'] = line_list[2]
            
            afr_df = pd.concat([afr_df, temp_df], ignore_index=True)
            
        genotype_df = pd.concat([genotype_df, afr_df], ignore_index=True)
        
    ###############################################################################
    
    if os.system(f'plink --bfile GWAS_PLINKs/AMR/chr{contig}_merged --list --snp {snp} --out plink_AMR') == 0:
    
        with open('plink_AMR.list') as file:
            data = [line.rstrip() for line in file]
        
        data = data[:3]
        
        amr_df = pd.DataFrame(columns=['IID', 'genotype'])
        for line in data:
            
            line_list = line.split(' ')
            
            temp_df = pd.DataFrame(list(set(line_list[3:])), columns=['IID'])
            
            temp_df['genotype'] = line_list[2]
            
            amr_df = pd.concat([amr_df, temp_df], ignore_index=True)
            
        genotype_df = pd.concat([genotype_df, amr_df], ignore_index=True)
        
    ###############################################################################
    
    plot_df = genotype_df.merge(merge_df, on='IID', how='left')
    
    plot_df.dropna(inplace=True)
    
    sample_size = len(plot_df)
    
    ###############################################################################
    
    fig, ax = plt.subplots()
    ax = sns.violinplot(data=plot_df, x='genotype', y='pheno', inner='quartile', color='.8')
    ax = sns.stripplot(data=plot_df, x='genotype', y='pheno', hue='ancestry', size=1.25)
    plt.legend('upper left', bbox_to_anchor=(1, 1))
    ax.legend(fontsize=5)
    ax.set_title(f'{snp} | N = {sample_size} | p-value = {pval}')
    ax.set_xlabel('genotype')
    ax.set_ylabel('mutation count')
    ax.figure.savefig(f'GWAS_EHv5_violin_plots/{snp_fn}_violin.pdf')