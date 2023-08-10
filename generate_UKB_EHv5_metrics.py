#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 23 09:58:31 2022

@author: williamlee
"""

import pandas as pd
import numpy as np
import os
import sys

wd = '/sc/arion/projects/sharpa01a/William/'
os.chdir(wd)

cohort = sys.argv[1]
contig = sys.argv[2]

contig_df = pd.read_csv(f'/sc/arion/projects/sharpa01a/TRE_DataRepository/WGS/hg38/EHv5_GenomewidePolymorphic/UKBB/PheWAS_GWAS_Tables_Draft2/{cohort}/{contig}.{cohort}_UKBB_EHv5_GenomewidePolymorphic_MediumTable_ForPheWAS_GWAS.tsv.gz', sep='\t', usecols=['VARID','RefUnit','SampleId','A1','A2'])

contig_df['tandem_repeat'] = contig_df['VARID'] + '_' + contig_df['RefUnit']
contig_df.drop(columns=['VARID', 'RefUnit'], inplace=True)

contig_df['AvgCopy'] = (contig_df['A1'] + contig_df['A2']) / 2
contig_df.drop(columns=['A1', 'A2'], inplace=True)

copy_df = contig_df[['tandem_repeat', 'AvgCopy']].groupby(by='tandem_repeat', as_index=False, sort=False).median()
copy_df.rename(columns={'AvgCopy':'PopCopy'}, inplace=True)

contig_df = contig_df.merge(copy_df, on='tandem_repeat', how='left')

del copy_df
contig_df.drop(columns=['tandem_repeat'], inplace=True)

contig_df['AvgCopy'] = contig_df['AvgCopy'].round()
contig_df['PopCopy'] = contig_df['PopCopy'].round()

contig_df['diff'] = contig_df['AvgCopy'] - contig_df['PopCopy']

contig_df.drop(columns=['AvgCopy', 'PopCopy'], inplace=True)

contig_df = contig_df[contig_df['diff'].notna()]

contig_df['pos_diff'] = contig_df['diff'].where(contig_df['diff'] > 0, np.nan)
contig_df['neg_diff'] = contig_df['diff'].where(contig_df['diff'] < 0, np.nan)

contig_count = contig_df.groupby(by='SampleId', as_index=False, sort=False).count()
contig_count.columns = ['sample_ID', 'num_TRs_total', 'num_expansions', 'num_contractions']

contig_sum = contig_df.groupby(by='SampleId', as_index=False, sort=False).sum()
contig_sum.columns = ['sample_ID', 'sum_total', 'sum_expansions', 'sum_contractions']

contig_df = contig_count.merge(contig_sum, on='sample_ID')

contig_df.to_excel(f'UKB_{cohort}_EHv5_metrics/{contig}_metrics.xlsx', index=False)