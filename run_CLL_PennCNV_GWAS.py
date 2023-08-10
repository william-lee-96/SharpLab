#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  5 16:08:28 2022

@author: williamlee
"""

import pandas as pd
import numpy as np
import os

wd = '/Users/williamlee/SharpLab/'
os.chdir(wd)

###############################################################################
###############################################################################
###############################################################################

CLL_df = pd.read_csv('out/rerun/out/CLL_report.rawcnv', delim_whitespace=True, header=None, usecols=[0,3,4])

CLL_df['sample'] = CLL_df[4].str.split('/').str[3]
CLL_df.drop(columns=4, inplace=True)

CLL_df['chr'] = CLL_df[0].str.split(':').str[0]
CLL_df['start'] = pd.to_numeric(CLL_df[0].str.split(':').str[1].str.split('-').str[0])
CLL_df['end'] = pd.to_numeric(CLL_df[0].str.split(':').str[1].str.split('-').str[1])
CLL_df.drop(columns=0, inplace=True)

CLL_df['cn'] = pd.to_numeric(CLL_df[3].str.split('=').str[1])
CLL_df.drop(columns=3, inplace=True)

###############################################################################

control_df = pd.read_csv('out/rerun/out/CONTROL_report.rawcnv', delim_whitespace=True, header=None, usecols=[0,3,4])

control_df['sample'] = control_df[4].str.split('/').str[3]
control_df.drop(columns=4, inplace=True)

control_df['chr'] = control_df[0].str.split(':').str[0]
control_df['start'] = pd.to_numeric(control_df[0].str.split(':').str[1].str.split('-').str[0])
control_df['end'] = pd.to_numeric(control_df[0].str.split(':').str[1].str.split('-').str[1])
control_df.drop(columns=0, inplace=True)

control_df['cn'] = pd.to_numeric(control_df[3].str.split('=').str[1])
control_df.drop(columns=3, inplace=True)

###############################################################################

cn_df = pd.concat([CLL_df, control_df], ignore_index=True)
cn_df['chr'] = pd.to_numeric(cn_df['chr'].str.lstrip('chr'))

del CLL_df; del control_df

###############################################################################
###############################################################################
###############################################################################

CLL_exclude_list = pd.read_csv('out/rerun/cases_to_exclude.tsv', sep='\t', header=None)[0].to_list()
CLL_fam_df = pd.read_csv('Data/CLL/ALL_NHL_OmniEx_CLL-cg3.fam', delim_whitespace=True, usecols=[1,4], names=['sample','sex'])
CLL_fam_df = CLL_fam_df[~CLL_fam_df['sample'].isin(CLL_exclude_list)]
CLL_fam_df['status'] = 1

control_exclude_list = pd.read_csv('out/rerun/controls_to_exclude.tsv', sep='\t', header=None)[0].to_list()
control_fam_df = pd.read_csv('Data/CLL/ALL_NHL_OmniEx_CONTROL-cg3.fam', delim_whitespace=True, usecols=[1,4], names=['sample','sex'])
control_fam_df = control_fam_df[~control_fam_df['sample'].isin(control_exclude_list)]
control_fam_df['status'] = 0

pheno_df = pd.concat([CLL_fam_df, control_fam_df], ignore_index=True)
pheno_df.set_index('sample', inplace=True)

del CLL_exclude_list; del control_exclude_list; del CLL_fam_df; del control_fam_df

###############################################################################
###############################################################################
###############################################################################

import statsmodels.api as sma

autosomes = (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22)

snp_df = pd.read_csv('out/rerun/OmniExpress_hg38.pfb', sep='\t', usecols=['Name', 'Chr', 'Position'])

y = pheno_df['status']
pheno_df.drop(columns='status', inplace=True)

cn_pval_list = []
del_pval_list = []
dup_pval_list = []
for autosome in autosomes:
    
    snp_temp_df = snp_df[snp_df['Chr'] == autosome]
    snp_temp_df.reset_index(drop=True, inplace=True)
    
    cn_temp_df = cn_df[cn_df['chr'] == autosome]
    cn_temp_df.reset_index(drop=True, inplace=True)
    
    for i in range(len(snp_temp_df)):
                        
        cn_temp_df['pos'] = snp_temp_df.loc[i, 'Position']
        
        cn_temp_df['pass'] = (cn_temp_df['start'] <= cn_temp_df['pos']) & (cn_temp_df['pos'] <= cn_temp_df['end'])
        
        cn_merge = cn_temp_df[cn_temp_df['pass'] == True]
        
        cn_merge = cn_merge[['sample', 'cn']]
        
        cn_merge.set_index('sample', inplace=True)
        
        del_merge = cn_merge[cn_merge['cn'] < 2]
        del_merge['cn'] = 1
        del_merge.rename(columns={'cn':'del'}, inplace=True)
        
        dup_merge = cn_merge[cn_merge['cn'] > 2]
        dup_merge['cn'] = 1
        dup_merge.rename(columns={'cn':'dup'}, inplace=True)
                
        pheno_df = pheno_df.join(cn_merge, how='left')
        pheno_df = pheno_df.join(del_merge, how='left')
        pheno_df = pheno_df.join(dup_merge, how='left')
        
        pheno_df['cn'].fillna(2, inplace=True)
        pheno_df['del'].fillna(0, inplace=True)
        pheno_df['dup'].fillna(0, inplace=True)
        
        X_cn = pheno_df[['cn', 'sex']]
        X_del = pheno_df[['del', 'sex']]
        X_dup = pheno_df[['dup', 'sex']]

        X_cn = sma.add_constant(X_cn)
        X_del = sma.add_constant(X_del)
        X_dup = sma.add_constant(X_dup)
        
        if len(X_cn['cn'].unique()) > 1:
            cn_pval_list.append(sma.Logit(y, X_cn).fit_regularized().pvalues[1])
        else:
            cn_pval_list.append(np.nan)
        
        if len(X_del['del'].unique()) > 1:
            del_pval_list.append(sma.Logit(y, X_del).fit_regularized().pvalues[1])
        else:
            del_pval_list.append(np.nan)
            
        if len(X_dup['dup'].unique()) > 1:
            dup_pval_list.append(sma.Logit(y, X_dup).fit_regularized().pvalues[1])
        else:
            dup_pval_list.append(np.nan)
        
        pheno_df.drop(columns=['cn', 'del', 'dup'], inplace=True)
        
###############################################################################

snp_df['P_cnv'] = cn_pval_list
snp_df['P_del'] = del_pval_list
snp_df['P_dup'] = dup_pval_list

snp_df.to_excel('out/CLL_SNP_logreg_pval.xlsx', index=False)

snp_df.to_csv('out/CLL_SNP_logreg_pval.tsv', sep='\t', index=False)