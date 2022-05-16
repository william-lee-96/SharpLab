#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  3 15:00:52 2022

@author: williamlee
"""

import pandas as pd
import numpy as np
import os

from scipy import stats

wd = '/Users/williamlee/SharpLab/'
os.chdir(wd)

###############################################################################
###############################################################################
###############################################################################

human_chr = ('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12',
             'chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX')

exon_coord_df = pd.read_csv('Data/nsSNV_TRE/Refseq_exon_hg38.txt.gz', sep='\t', header=None)
exon_coord_df.columns = ['chr','start','end','gene','strand','region']
exon_coord_df.drop(columns=['strand','region'], inplace=True)

exon_coord_df['gene'] = exon_coord_df['gene'].str.split('_').str[0]

mda_hgnc = pd.read_csv('Data/nsSNV_TRE/MDA_to_HGNC.csv', sep=',', skiprows=1)

mda_hgnc = mda_hgnc[mda_hgnc['Match type'].isin(['Approved symbol', 'Previous symbol'])] # (219, 6)

mda_merge = mda_hgnc['Approved symbol']
mda_merge.name = 'gene'

dva_exon_df = exon_coord_df.merge(mda_merge, on='gene')

dva_exon_df = dva_exon_df[dva_exon_df['chr'].isin(human_chr)]

dva_exon_df.drop_duplicates(inplace=True, ignore_index=True)

dva_exon_df['start'] = dva_exon_df['start'] - 10
dva_exon_df['end'] = dva_exon_df['end'] + 10

dva_exon_df['chr'] = dva_exon_df['chr'].str.replace('chr', '')

###############################################################################
###############################################################################
###############################################################################

top_genes = ['SLX1A', 'PNKP', 'XAB2', 'ALKBH2', 'RNF168']

eur_cutoff = 27912*0.85
afr_cutoff = 10623*0.85
amr_cutoff = 4880*0.85

expansion_df = pd.read_excel('Data/nsSNV_TRE/DVA_case_TREs.xlsx')
expansion_df.columns = ['sample_id', 'TRE']

###############################################################################
###############################################################################
###############################################################################

eur_df = pd.read_csv('out/REGENIE_step2_EUR_DVA_minMAC1_pheno.regenie', delim_whitespace=True)

eur_df = eur_df[eur_df['N'] > eur_cutoff]
eur_df.reset_index(drop=True, inplace=True)

eur_df['P-value'] = 10**(-eur_df['LOG10P'])

for row in range(len(eur_df)):
    if np.isnan(eur_df.loc[row,'P-value']):
        eur_df.loc[row,'P-value'] = stats.norm.sf(abs(eur_df.loc[row,'BETA']/eur_df.loc[row,'SE']))*2

eur_df = eur_df[['CHROM', 'GENPOS', 'ID', 'N', 'A1FREQ', 'P-value']]

eur_df['CHROM'] = eur_df['CHROM'].astype('string')

eur_df['CHROM'].replace('23', 'X', inplace=True)

eur_df.columns = ['chr', 'pos', 'ID', 'N', 'MAF', 'P-value']

eur_df['gene'] = np.nan
for row in range(len(eur_df)):
    eur_df.loc[row, 'gene'] = dva_exon_df[(dva_exon_df['chr'] == eur_df.loc[row,'chr']) &
                                          (dva_exon_df['start'] <= eur_df.loc[row,'pos']) &
                                          (dva_exon_df['end'] >= eur_df.loc[row,'pos'])].iloc[0,3]

###############################################################################

afr_df = pd.read_csv('out/REGENIE_step2_AFR_DVA_minMAC1_pheno.regenie', delim_whitespace=True)

afr_df = afr_df[afr_df['N'] > afr_cutoff]
afr_df.reset_index(drop=True, inplace=True)

afr_df['P-value'] = 10**(-afr_df['LOG10P'])

for row in range(len(afr_df)):
    if np.isnan(afr_df.loc[row,'P-value']):
        afr_df.loc[row,'P-value'] = stats.norm.sf(abs(afr_df.loc[row,'BETA']/afr_df.loc[row,'SE']))*2

afr_df = afr_df[['CHROM', 'GENPOS', 'ID', 'N', 'A1FREQ', 'P-value']]

afr_df['CHROM'] = afr_df['CHROM'].astype('string')

afr_df['CHROM'].replace('23', 'X', inplace=True)

afr_df.columns = ['chr', 'pos', 'ID', 'N', 'MAF', 'P-value']

afr_df['gene'] = np.nan
for row in range(len(afr_df)):
    afr_df.loc[row, 'gene'] = dva_exon_df[(dva_exon_df['chr'] == afr_df.loc[row,'chr']) &
                                          (dva_exon_df['start'] <= afr_df.loc[row,'pos']) &
                                          (dva_exon_df['end'] >= afr_df.loc[row,'pos'])].iloc[0,3]
    
###############################################################################

amr_df = pd.read_csv('out/REGENIE_step2_AMR_DVA_minMAC1_pheno.regenie', delim_whitespace=True)

amr_df = amr_df[amr_df['N'] > amr_cutoff]
amr_df.reset_index(drop=True, inplace=True)

amr_df['P-value'] = 10**(-amr_df['LOG10P'])

for row in range(len(amr_df)):
    if np.isnan(amr_df.loc[row,'P-value']):
        amr_df.loc[row,'P-value'] = stats.norm.sf(abs(amr_df.loc[row,'BETA']/amr_df.loc[row,'SE']))*2

amr_df = amr_df[['CHROM', 'GENPOS', 'ID', 'N', 'A1FREQ', 'P-value']]

amr_df['CHROM'] = amr_df['CHROM'].astype('string')

amr_df['CHROM'].replace('23', 'X', inplace=True)

amr_df.columns = ['chr', 'pos', 'ID', 'N', 'MAF', 'P-value']

amr_df['gene'] = np.nan
for row in range(len(amr_df)):
    amr_df.loc[row, 'gene'] = dva_exon_df[(dva_exon_df['chr'] == amr_df.loc[row,'chr']) &
                                          (dva_exon_df['start'] <= amr_df.loc[row,'pos']) &
                                          (dva_exon_df['end'] >= amr_df.loc[row,'pos'])].iloc[0,3]

###############################################################################
###############################################################################
###############################################################################

eur_ped_df = pd.read_csv('Data/nsSNV_TRE/DVA_PED_files/DVA_EUR_REGENIE_input.ped', sep='\t', header=None)

eur_map_df = pd.read_csv('Data/nsSNV_TRE/DVA_PED_files/DVA_EUR_REGENIE_input.map', sep='\t', header=None)

rename_dict = dict(zip(list(eur_ped_df.iloc[:,6:].columns), list(eur_map_df[1])))

eur_ped_df.rename(columns=rename_dict, inplace=True)

eur_ped_df.drop(columns=[1,2,3,4], inplace=True)

eur_ped_df.rename(columns={0:'sample_id', 5:'phenotype'}, inplace=True)

###############################################################################

afr_ped_df = pd.read_csv('Data/nsSNV_TRE/DVA_PED_files/DVA_AFR_REGENIE_input.ped', sep='\t', header=None)

afr_map_df = pd.read_csv('Data/nsSNV_TRE/DVA_PED_files/DVA_AFR_REGENIE_input.map', sep='\t', header=None)

rename_dict = dict(zip(list(afr_ped_df.iloc[:,6:].columns), list(afr_map_df[1])))

afr_ped_df.rename(columns=rename_dict, inplace=True)

afr_ped_df.drop(columns=[1,2,3,4], inplace=True)

afr_ped_df.rename(columns={0:'sample_id', 5:'phenotype'}, inplace=True)

###############################################################################

amr_ped_df = pd.read_csv('Data/nsSNV_TRE/DVA_PED_files/DVA_AMR_REGENIE_input.ped', sep='\t', header=None)

amr_map_df = pd.read_csv('Data/nsSNV_TRE/DVA_PED_files/DVA_AMR_REGENIE_input.map', sep='\t', header=None)

rename_dict = dict(zip(list(amr_ped_df.iloc[:,6:].columns), list(amr_map_df[1])))

amr_ped_df.rename(columns=rename_dict, inplace=True)

amr_ped_df.drop(columns=[1,2,3,4], inplace=True)

amr_ped_df.rename(columns={0:'sample_id', 5:'phenotype'}, inplace=True)

###############################################################################
###############################################################################
###############################################################################

skato_var_list = []
for gene in top_genes:

    ###########################################################################

    eur_gene_df = eur_df[eur_df['gene']==gene]

    case_counts = []
    ctrl_counts = []
    tre_list = []
    for snv in eur_gene_df['ID']:
        
        alt_allele = snv.split('_')[1]
        
        snv_df = eur_ped_df[['sample_id', 'phenotype', snv]]
        
        carrier_df = snv_df[snv_df[snv].str.contains(alt_allele, regex=False)]
        
        case_counts.append(sum(carrier_df['phenotype']==2))
        ctrl_counts.append(sum(carrier_df['phenotype']==1))
        
        case_df = carrier_df[carrier_df['phenotype']==2]
            
        if not case_df.empty:
            
            tre_df = case_df.merge(expansion_df, on='sample_id')
            
            tre_list.append(tre_df['TRE'].value_counts().to_dict())
            
        else:
            
            tre_list.append({})
        
    eur_gene_df['N_case_carriers'] = case_counts
    eur_gene_df['N_control_carriers'] = ctrl_counts
    eur_gene_df['TREs_in_cases'] = tre_list

    eur_gene_df.drop(columns='gene', inplace=True)
    eur_gene_df.sort_values(by='P-value', inplace=True)

    ###########################################################################

    afr_gene_df = afr_df[afr_df['gene']==gene]

    case_counts = []
    ctrl_counts = []
    tre_list = []
    for snv in afr_gene_df['ID']:
        
        alt_allele = snv.split('_')[1]
        
        snv_df = afr_ped_df[['sample_id', 'phenotype', snv]]
        
        carrier_df = snv_df[snv_df[snv].str.contains(alt_allele, regex=False)]
        
        case_counts.append(sum(carrier_df['phenotype']==2))
        ctrl_counts.append(sum(carrier_df['phenotype']==1))
        
        case_df = carrier_df[carrier_df['phenotype']==2]
            
        if not case_df.empty:
            
            tre_df = case_df.merge(expansion_df, on='sample_id')
            
            tre_list.append(tre_df['TRE'].value_counts().to_dict())
            
        else:
            
            tre_list.append({})
        
    afr_gene_df['N_case_carriers'] = case_counts
    afr_gene_df['N_control_carriers'] = ctrl_counts
    afr_gene_df['TREs_in_cases'] = tre_list

    afr_gene_df.drop(columns='gene', inplace=True)
    afr_gene_df.sort_values(by='P-value', inplace=True)

    ###########################################################################

    amr_gene_df = amr_df[amr_df['gene']==gene]

    case_counts = []
    ctrl_counts = []
    tre_list = []
    for snv in amr_gene_df['ID']:
        
        alt_allele = snv.split('_')[1]
        
        snv_df = amr_ped_df[['sample_id', 'phenotype', snv]]
        
        carrier_df = snv_df[snv_df[snv].str.contains(alt_allele, regex=False)]
        
        case_counts.append(sum(carrier_df['phenotype']==2))
        ctrl_counts.append(sum(carrier_df['phenotype']==1))
        
        case_df = carrier_df[carrier_df['phenotype']==2]
            
        if not case_df.empty:
            
            tre_df = case_df.merge(expansion_df, on='sample_id')
            
            tre_list.append(tre_df['TRE'].value_counts().to_dict())
            
        else:
            
            tre_list.append({})
        
    amr_gene_df['N_case_carriers'] = case_counts
    amr_gene_df['N_control_carriers'] = ctrl_counts
    amr_gene_df['TREs_in_cases'] = tre_list

    amr_gene_df.drop(columns='gene', inplace=True)
    amr_gene_df.sort_values(by='P-value', inplace=True)

    ###############################################################################
    
    skato_var_list.extend(eur_gene_df['ID'].tolist())
    skato_var_list.extend(afr_gene_df['ID'].tolist())
    skato_var_list.extend(amr_gene_df['ID'].tolist())
        
    ###############################################################################

    with pd.ExcelWriter(f'out/SKATO_top_hits/{gene}_SKATO_input.xlsx') as writer:  
        eur_gene_df.to_excel(writer, sheet_name='EUR', index=False)
        afr_gene_df.to_excel(writer, sheet_name='AFR', index=False)
        amr_gene_df.to_excel(writer, sheet_name='AMR', index=False)

###############################################################################
###############################################################################
###############################################################################

skato_var_df = pd.DataFrame(list(set(skato_var_list)), columns=['merge_col'])

skato_var_df.to_excel('out/SKATO_top5_for_scoring.xlsx', index=False)

skato_var_df['chr'] = pd.to_numeric(skato_var_df['merge_col'].str.split(':').apply(lambda x: x[0]))

skato_var_df['start'] = skato_var_df['merge_col'].str.split(':').apply(lambda x: x[1])
skato_var_df['start'] = pd.to_numeric(skato_var_df['start'].str.split('_').apply(lambda x: x[0]))

skato_var_df['end'] = skato_var_df['start']

skato_var_df['ref'] = skato_var_df['merge_col'].str.split('_').apply(lambda x: x[2])
skato_var_df['alt'] = skato_var_df['merge_col'].str.split('_').apply(lambda x: x[1])

skato_var_df = skato_var_df[['chr','start','end','ref','alt','merge_col']]

skato_var_df.to_csv('out/SKATO_top5_for_annotation.txt', sep='\t', index=False)