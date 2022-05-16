#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  7 12:05:06 2021

@author: williamlee
"""

import pandas as pd
import numpy as np
import os
import sys

import seaborn as sns

wd = '/sc/arion/projects/sharpa01a/William/'
os.chdir(wd)

tre_df = pd.read_csv('TOPMed_ListOf14PathogenicLoci_hg38.tsv', sep='\t')

contig = sys.argv[1]

for file in os.scandir('REGENIE_meta/'):
    if file.name.endswith('.txt') and file.name.split('_')[1] == contig:
        
        chr_df = pd.read_csv(file.path, sep='\t')
        
        chr_df['pos'] = chr_df['MarkerName'].str.split(':').apply(lambda x: x[1])
        chr_df['pos'] = pd.to_numeric(chr_df['pos'].str.split('_').apply(lambda x: x[0]))

        chr_df['legend'] = 'normal'

        chr_df.loc[chr_df['P-value'] < 1E-6, 'legend'] = 'P < 1E-6'
        chr_df.loc[chr_df['P-value'] < 5E-8, 'legend'] = 'P < 5E-8'

        temp_df = tre_df[tre_df['chrom'] == contig]
        if not temp_df.empty:
            for row in range(len(temp_df)):
                
                start = temp_df.iloc[row, 1] - 1000000
                end = temp_df.iloc[row, 2] + 1000000
                
                chr_df.loc[(chr_df['pos'] > start) & (chr_df['pos'] < end), 'legend'] = 'within TRE +/- 1Mb'
                
        custom_palette = {}
        for key in chr_df['legend'].unique():
            if key == 'within TRE +/- 1Mb':
                custom_palette[key] = 'green'
            elif key == 'P < 5E-8':
                custom_palette[key] = 'red'
            elif key == 'P < 1E-6':
                custom_palette[key] = 'orange'
            else:
                custom_palette[key] = 'gray'

        ax = sns.scatterplot(x=chr_df['pos'], y=-np.log10(chr_df['P-value']), hue=chr_df['legend'], palette=custom_palette)
        ax.set_xlabel('position')
        ax.set_ylabel('-log10(P)')
        ax.set_title(f'{contig}')
        ax.figure.savefig(f'REGENIE_METAL_plots/{contig}_GWAS.png', dpi=300)