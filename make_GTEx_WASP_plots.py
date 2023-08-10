#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 10 11:07:03 2023

@author: williamlee
"""

import os
import pandas as pd
import numpy as np

import seaborn as sns
import matplotlib.pyplot as plt

os.chdir('/Users/williamlee/SharpLab')

wasp_df = pd.read_excel('Data/GTEx/GTEx_WASP_25_genes.xlsx')

wasp_df = wasp_df[wasp_df['TOTAL_COUNT'] >= 10]

###############################################################################

ADARB1_df = wasp_df[wasp_df['GENE_NAME'] == 'ADARB1']

ADARB1_df = ADARB1_df[ADARB1_df['TISSUE_ID'] == 'ARTTBL']

ADARB1_df['STATUS'] = np.nan

ADARB1_df.loc[ADARB1_df.SUBJECT_ID == 'GTEX-1CB4H', 'STATUS'] = 'case'
ADARB1_df.loc[ADARB1_df.SUBJECT_ID != 'GTEX-1CB4H', 'STATUS'] = 'control'

var_list = list(ADARB1_df[ADARB1_df['STATUS']=='case']['VARIANT_ID'].unique())

ADARB1_df = ADARB1_df[ADARB1_df['VARIANT_ID'].isin(var_list)]

ADARB1_df.sort_values(by=['POS', 'STATUS'], inplace=True)

# ADARB1_df[ADARB1_df['STATUS']=='case']['SAMPLE_ID'].unique()
# array(['GTEX-1CB4H-0826-SM-7SB7Z'], dtype=object)

# generate and save scatterplot
fig, ax = plt.subplots()
sns.scatterplot(data=ADARB1_df[ADARB1_df['STATUS']=='control'], x='REF_COUNT', y='ALT_COUNT', fc='none', ec='#a2bffe', legend=False, ax=ax)
sns.scatterplot(data=ADARB1_df[ADARB1_df['STATUS']=='case'], x='REF_COUNT', y='ALT_COUNT', alpha=0.75, ec='black', legend=False, ax=ax)
ax.set_xlabel('reference count')
ax.set_ylabel('alternate count')
ax.set_title('ADARB1: reference count vs. alternate count\nepimutation carrier(s): GTEX-1CB4H-0826-SM-7SB7Z')
ax.figure.savefig('out/GTEx_WASP_scatterplots/ADARB1_GTEx_WASP_scatterplot.pdf')

del ADARB1_df

###############################################################################

CSNK1E_df = wasp_df[wasp_df['GENE_NAME'] == 'CSNK1E']

CSNK1E_df = CSNK1E_df[CSNK1E_df['TISSUE_ID'] == 'ESPMCS']

CSNK1E_df['STATUS'] = np.nan

CSNK1E_df.loc[CSNK1E_df.SUBJECT_ID == 'GTEX-X62O', 'STATUS'] = 'case'
CSNK1E_df.loc[CSNK1E_df.SUBJECT_ID != 'GTEX-X62O', 'STATUS'] = 'control'

var_list = list(CSNK1E_df[CSNK1E_df['STATUS']=='case']['VARIANT_ID'].unique())

CSNK1E_df = CSNK1E_df[CSNK1E_df['VARIANT_ID'].isin(var_list)]

CSNK1E_df.sort_values(by=['POS', 'STATUS'], inplace=True)

# CSNK1E_df[CSNK1E_df['STATUS']=='case']['SAMPLE_ID'].unique()
# array(['GTEX-X62O-2226-SM-46MW3'], dtype=object)

# generate and save scatterplot
fig, ax = plt.subplots()
sns.scatterplot(data=CSNK1E_df[CSNK1E_df['STATUS']=='control'], x='REF_COUNT', y='ALT_COUNT', fc='none', ec='#a2bffe', legend=False, ax=ax)
sns.scatterplot(data=CSNK1E_df[CSNK1E_df['STATUS']=='case'], x='REF_COUNT', y='ALT_COUNT', alpha=0.75, ec='black', legend=False, ax=ax)
ax.set_xlabel('reference count')
ax.set_ylabel('alternate count')
ax.set_title('CSNK1E: reference count vs. alternate count\nepimutation carrier(s): GTEX-X62O-2226-SM-46MW3')
ax.figure.savefig('out/GTEx_WASP_scatterplots/CSNK1E_GTEx_WASP_scatterplot.pdf')

del CSNK1E_df

###############################################################################

FRA10AC1_df = wasp_df[wasp_df['GENE_NAME'] == 'FRA10AC1']

FRA10AC1_df = FRA10AC1_df[FRA10AC1_df['TISSUE_ID'] == 'NERVET']

FRA10AC1_df['STATUS'] = np.nan

subject_list = ['GTEX-15EO6', 'GTEX-16AAH', 'GTEX-14JG1', 'GTEX-ZP4G']

FRA10AC1_df.loc[FRA10AC1_df.SUBJECT_ID.isin(subject_list), 'STATUS'] = 'case'
FRA10AC1_df.loc[~FRA10AC1_df.SUBJECT_ID.isin(subject_list), 'STATUS'] = 'control'

var_list = list(FRA10AC1_df[FRA10AC1_df['STATUS']=='case']['VARIANT_ID'].unique())

FRA10AC1_df = FRA10AC1_df[FRA10AC1_df['VARIANT_ID'].isin(var_list)]

FRA10AC1_df.sort_values(by=['POS', 'STATUS'], inplace=True)

FRA10AC1_df[FRA10AC1_df['STATUS']=='case']['SAMPLE_ID'].unique()
# array(['GTEX-15EO6-2626-SM-6PALK', 'GTEX-16AAH-2126-SM-7LG4D', 'GTEX-ZP4G-2226-SM-57WFB'], dtype=object)

# generate and save scatterplot
fig, ax = plt.subplots()
sns.scatterplot(data=FRA10AC1_df[FRA10AC1_df['STATUS']=='control'], x='REF_COUNT', y='ALT_COUNT', fc='none', ec='#a2bffe', legend=False, ax=ax)
sns.scatterplot(data=FRA10AC1_df[FRA10AC1_df['STATUS']=='case'], x='REF_COUNT', y='ALT_COUNT', alpha=0.75, ec='black', legend=False, ax=ax)
ax.set_xlabel('reference count')
ax.set_ylabel('alternate count')
ax.set_title('FRA10AC1: reference count vs. alternate count\nepimutation carrier(s): GTEX-15EO6-2626-SM-6PALK,\nGTEX-16AAH-2126-SM-7LG4D, GTEX-ZP4G-2226-SM-57WFB')
ax.figure.savefig('out/GTEx_WASP_scatterplots/FRA10AC1_GTEx_WASP_scatterplot.pdf', bbox_inches='tight')

del FRA10AC1_df

###############################################################################

PCMTD2_df = wasp_df[wasp_df['GENE_NAME'] == 'PCMTD2']

PCMTD2_df = PCMTD2_df[PCMTD2_df['TISSUE_ID'] == 'ADPVSC']

PCMTD2_df['STATUS'] = np.nan

PCMTD2_df.loc[PCMTD2_df.SUBJECT_ID == 'GTEX-13D11', 'STATUS'] = 'case'
PCMTD2_df.loc[PCMTD2_df.SUBJECT_ID != 'GTEX-13D11', 'STATUS'] = 'control'

var_list = list(PCMTD2_df[PCMTD2_df['STATUS']=='case']['VARIANT_ID'].unique())

PCMTD2_df = PCMTD2_df[PCMTD2_df['VARIANT_ID'].isin(var_list)]

PCMTD2_df.sort_values(by=['POS', 'STATUS'], inplace=True)

PCMTD2_df[PCMTD2_df['STATUS']=='case']['SAMPLE_ID'].unique()
# array(['GTEX-13D11-0526-SM-5LZYM'], dtype=object)

# generate and save scatterplot
fig, ax = plt.subplots()
sns.scatterplot(data=PCMTD2_df[PCMTD2_df['STATUS']=='control'], x='REF_COUNT', y='ALT_COUNT', fc='none', ec='#a2bffe', legend=False, ax=ax)
sns.scatterplot(data=PCMTD2_df[PCMTD2_df['STATUS']=='case'], x='REF_COUNT', y='ALT_COUNT', alpha=0.75, ec='black', legend=False, ax=ax)
ax.set_xlabel('reference count')
ax.set_ylabel('alternate count')
ax.set_title('PCMTD2: reference count vs. alternate count\nepimutation carrier(s): GTEX-13D11-0526-SM-5LZYM')
ax.figure.savefig('out/GTEx_WASP_scatterplots/PCMTD2_GTEx_WASP_scatterplot.pdf')

del PCMTD2_df

###############################################################################

ZNF713_df = wasp_df[wasp_df['GENE_NAME'] == 'ZNF713']

ZNF713_df = ZNF713_df[ZNF713_df['TISSUE_ID'] == 'PTTARY']

ZNF713_df['STATUS'] = np.nan

subject_list = ['GTEX-18D9U', 'GTEX-1EN7A', 'GTEX-1F75B', 'GTEX-1JK1U', 'GTEX-13NZA', 'GTEX-WH7G']

ZNF713_df.loc[ZNF713_df.SUBJECT_ID.isin(subject_list), 'STATUS'] = 'case'
ZNF713_df.loc[~ZNF713_df.SUBJECT_ID.isin(subject_list), 'STATUS'] = 'control'

var_list = list(ZNF713_df[ZNF713_df['STATUS']=='case']['VARIANT_ID'].unique())

ZNF713_df = ZNF713_df[ZNF713_df['VARIANT_ID'].isin(var_list)]

ZNF713_df.sort_values(by=['POS', 'STATUS'], inplace=True)

ZNF713_df[ZNF713_df['STATUS']=='case']['SAMPLE_ID'].unique()
# array(['GTEX-1F75B-3126-SM-9JGGE'], dtype=object)

# generate and save scatterplot
fig, ax = plt.subplots()
sns.scatterplot(data=ZNF713_df[ZNF713_df['STATUS']=='control'], x='REF_COUNT', y='ALT_COUNT', fc='none', ec='#a2bffe', legend=False, ax=ax)
sns.scatterplot(data=ZNF713_df[ZNF713_df['STATUS']=='case'], x='REF_COUNT', y='ALT_COUNT', alpha=0.75, ec='black', legend=False, ax=ax)
ax.set_xlabel('reference count')
ax.set_ylabel('alternate count')
ax.set_title('ZNF713: reference count vs. alternate count\nepimutation carrier(s): GTEX-1F75B-3126-SM-9JGGE')
ax.figure.savefig('out/GTEx_WASP_scatterplots/ZNF713_GTEx_WASP_scatterplot.pdf')

del ZNF713_df