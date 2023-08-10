#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 27 16:10:11 2022

@author: williamlee
"""

import pandas as pd
import numpy as np
import os

wd = '/Users/williamlee/SharpLab/'
os.chdir(wd)

autosomes = ('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22')

sex_chr_list = [23, 24, 25]

###############################################################################
###############################################################################
###############################################################################

CLL_imiss_df = pd.read_csv('out/QC/CLL/CLL.imiss', delim_whitespace=True, usecols=['IID','F_MISS'])
CLL_imiss_df = CLL_imiss_df[CLL_imiss_df['F_MISS'] > 0.07]

CLL_het_df = pd.read_csv('out/QC/CLL/CLL.het', delim_whitespace=True, usecols=['IID','O(HOM)','N(NM)'])
CLL_het_df['mean_het'] = 1 - (CLL_het_df['O(HOM)']/CLL_het_df['N(NM)'])
CLL_het_df = CLL_het_df[(CLL_het_df['mean_het'] < 0.25) | (CLL_het_df['mean_het'] > 0.33)]

CLL_sex_df = pd.read_csv('out/QC/CLL/CLL.sexcheck', delim_whitespace=True, usecols=['IID','STATUS'])
CLL_sex_df = CLL_sex_df[CLL_sex_df['STATUS'] == 'PROBLEM']

CLL_genome_df = pd.read_csv('out/QC/CLL/CLL.genome', delim_whitespace=True, usecols=['IID1','IID2','PI_HAT'])
CLL_genome_df = CLL_genome_df[(CLL_genome_df['PI_HAT'] > 0.40)]

CLL_exclude_list = CLL_imiss_df['IID'].to_list() + CLL_het_df['IID'].to_list() + CLL_sex_df['IID'].to_list() + CLL_genome_df['IID1'].to_list() + CLL_genome_df['IID2'].to_list() 

CLL_exclude_list = list(set(CLL_exclude_list))

pd.Series(CLL_exclude_list).to_csv('out/rerun/cases_to_exclude.tsv', sep='\t', index=False, header=False)

del CLL_imiss_df; del CLL_het_df; del CLL_sex_df; del CLL_genome_df; del CLL_exclude_list

###############################################################################

control_imiss_df = pd.read_csv('out/QC/control/control.imiss', delim_whitespace=True, usecols=['IID','F_MISS'])
control_imiss_df = control_imiss_df[control_imiss_df['F_MISS'] > 0.07]

control_het_df = pd.read_csv('out/QC/control/control.het', delim_whitespace=True, usecols=['IID','O(HOM)','N(NM)'])
control_het_df['mean_het'] = 1 - (control_het_df['O(HOM)']/control_het_df['N(NM)'])
control_het_df = control_het_df[(control_het_df['mean_het'] < 0.25) | (control_het_df['mean_het'] > 0.33)]

control_sex_df = pd.read_csv('out/QC/control/control.sexcheck', delim_whitespace=True, usecols=['IID','STATUS'])
control_sex_df = control_sex_df[control_sex_df['STATUS'] == 'PROBLEM']

control_genome_df = pd.read_csv('out/QC/control/control.genome', delim_whitespace=True, usecols=['IID1','IID2','PI_HAT'])
control_genome_df = control_genome_df[(control_genome_df['PI_HAT'] > 0.40)]

control_exclude_list = control_imiss_df['IID'].to_list() + control_het_df['IID'].to_list() + control_sex_df['IID'].to_list() + control_genome_df['IID1'].to_list() + control_genome_df['IID2'].to_list() 

control_exclude_list = list(set(control_exclude_list))

pd.Series(control_exclude_list).to_csv('out/rerun/controls_to_exclude.tsv', sep='\t', index=False, header=False)

del control_imiss_df; del control_het_df; del control_sex_df; del control_genome_df; del control_exclude_list

###############################################################################
###############################################################################
###############################################################################

CLL_bim_df = pd.read_csv('Data/CLL/ALL_NHL_OmniEx_CLL-cg3.bim', delim_whitespace=True, header=None, usecols=[0,1,3])
control_bim_df = pd.read_csv('Data/CLL/ALL_NHL_OmniEx_CONTROL-cg3.bim', delim_whitespace=True, header=None, usecols=[0,1,3])

combined_bim_df = pd.concat([CLL_bim_df, control_bim_df], ignore_index=True)
combined_bim_df.drop_duplicates(inplace=True)
combined_bim_df.columns = ['Chr', 'Name', 'Position']

combined_bim_df = combined_bim_df[~combined_bim_df['Chr'].isin(sex_chr_list)]

del CLL_bim_df; del control_bim_df

###############################################################################
###############################################################################
###############################################################################

combined_hwe_df = pd.read_csv('out/CLL_control_merged.hwe', delim_whitespace=True, usecols=['CHR', 'SNP', 'P'])

combined_hwe_df = combined_hwe_df[~combined_hwe_df['CHR'].isin(sex_chr_list)]

combined_hwe_df = combined_hwe_df[combined_hwe_df['P'] > 0.001]

hwe_keep_list = combined_hwe_df['SNP'].to_list()

del combined_hwe_df

###############################################################################

combined_lmiss_df = pd.read_csv('out/CLL_control_merged.lmiss', delim_whitespace=True, usecols=['CHR', 'SNP', 'F_MISS'])

combined_lmiss_df = combined_lmiss_df[~combined_lmiss_df['CHR'].isin(sex_chr_list)]

combined_lmiss_df = combined_lmiss_df[combined_lmiss_df['F_MISS'] < 0.05]

lmiss_keep_list = combined_lmiss_df['SNP'].to_list()

del combined_lmiss_df

###############################################################################

combined_bim_df = combined_bim_df[combined_bim_df['Name'].isin(hwe_keep_list)]
combined_bim_df = combined_bim_df[combined_bim_df['Name'].isin(lmiss_keep_list)]

combined_bim_df['Name'].to_csv('out/rsIDs_to_convert.tsv', sep='\t', header=False, index=False)

###############################################################################
###############################################################################
###############################################################################

hg38_dup_df = pd.read_excel('../Downloads/GRCh38GenomicSuperDup.xlsx', header=None, usecols=[0,1,2])
hg38_dup_df.columns = ['chr', 'start', 'end']
hg38_dup_df = hg38_dup_df[hg38_dup_df['chr'].isin(autosomes)]

hg38_df = pd.read_csv('out/CLL_dbSNP155_hg38.tsv', sep='\t', header=None)
hg38_df.columns = ['chr', 'pos', 'rsID']
hg38_df = hg38_df[hg38_df['chr'].isin(autosomes)]

temp_df = pd.DataFrame()
for autosome in autosomes:
    
    hg38_temp_df = hg38_df[hg38_df['chr'] == autosome]
    hg38_temp_df.reset_index(drop=True, inplace=True)
    
    hg38_dup_temp_df = hg38_dup_df[hg38_dup_df['chr'] == autosome]
    hg38_dup_temp_df.reset_index(drop=True, inplace=True)
    
    dup_list = []
    for i in range(len(hg38_temp_df)):
        
        snv_pos = hg38_temp_df.loc[i, 'pos']
        
        hg38_dup_temp_df['pos'] = snv_pos
        
        hg38_dup_temp_df['pass'] = (hg38_dup_temp_df['start'] <= hg38_dup_temp_df['pos']) & (hg38_dup_temp_df['pos'] <= hg38_dup_temp_df['end'])
        
        dup_list.append(hg38_dup_temp_df['pass'].any())
        
    hg38_temp_df['dup'] = dup_list
    
    temp_df = pd.concat([temp_df, hg38_temp_df], ignore_index=True)

temp_df = temp_df[temp_df['dup'] == False]

temp_df.drop(columns='dup', inplace=True)

del hg38_dup_df; del hg38_dup_temp_df; del hg38_temp_df; del dup_list; del i; del snv_pos

###############################################################################
###############################################################################
###############################################################################

pfb_df = pd.read_csv('../Downloads/OmniExpress_hg19.pfb', sep='\t', usecols=['Name', 'PFB'])

temp_df.columns = ['Chr', 'Position', 'Name']

pfb_df = temp_df.merge(pfb_df, on='Name', how='left')

pfb_df['Chr'] = pd.to_numeric(pfb_df['Chr'].str.lstrip('chr'))

pfb_df = pfb_df[['Name', 'Chr', 'Position', 'PFB']]

pfb_df = pfb_df[pfb_df['PFB'].notna()]

pfb_df.to_csv('out/rerun/OmniExpress_hg38.pfb', sep='\t', index=False)

###############################################################################
###############################################################################
###############################################################################