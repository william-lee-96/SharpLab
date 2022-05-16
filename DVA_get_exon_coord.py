#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 26 12:05:19 2022

@author: williamlee
"""

import pandas as pd
import numpy as np
import os

wd = '/Users/williamlee/SharpLab/'
os.chdir(wd)

human_chr = ('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12',
             'chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX')

###############################################################################

exon_coord_df = pd.read_csv('Data/nsSNV_TRE/Refseq_exon_hg38.txt.gz', sep='\t', header=None)
exon_coord_df.columns = ['chr','start','end','gene','strand','region']
exon_coord_df.drop(columns=['strand','region'], inplace=True)

exon_coord_df['gene'] = exon_coord_df['gene'].str.split('_').str[0]

###############################################################################

mda_hgnc = pd.read_csv('Data/nsSNV_TRE/MDA_to_HGNC.csv', sep=',', skiprows=1)

mda_hgnc = mda_hgnc[mda_hgnc['Match type'].isin(['Approved symbol', 'Previous symbol'])] # (219, 6)

mda_merge = mda_hgnc['Approved symbol']
mda_merge.name = 'gene'

###############################################################################

dva_exon_df = exon_coord_df.merge(mda_merge, on='gene')

dva_exon_df = dva_exon_df[dva_exon_df['chr'].isin(human_chr)]

dva_exon_df.drop_duplicates(inplace=True, ignore_index=True)

dva_exon_df['start'] = dva_exon_df['start'] - 10
dva_exon_df['end'] = dva_exon_df['end'] + 10

dva_exon_df.drop(columns='gene', inplace=True)

dva_exon_df['chr'] = dva_exon_df['chr'].str.replace('chr', '')

label_col = []
for idx in range(len(dva_exon_df)):
    label_col.append('R' + str(idx))

dva_exon_df['label'] = label_col

dva_exon_df.to_csv('out/DVA_exon_ranges.tsv', sep='\t', header=False, index=False)