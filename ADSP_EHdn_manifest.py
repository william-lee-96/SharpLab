#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  5 14:15:01 2021

@author: williamlee
"""

import os
import pandas as pd

# set wd
wd = '/Users/williamlee/SharpLab/'
os.chdir(wd)

###############################################################################
###############################################################################
###############################################################################

adsp_case_ctrl_df = pd.read_excel('Data/ADSP/ADSP_14239WGS_FiltCaseControlSamplesWithAnsPCsSequencingInfo_4August21.xlsx')

adsp_case_ctrl_df.dropna(subset=['Status'], inplace=True)
adsp_case_ctrl_df = adsp_case_ctrl_df[adsp_case_ctrl_df['ReadLength']==150]
adsp_case_ctrl_df = adsp_case_ctrl_df[~adsp_case_ctrl_df['Ancestry'].isin(['EAS', 'SAS'])]

adsp_ehdn_out_df = pd.read_csv('Data/ADSP/ADSP_EHdn_out.txt', sep='\t', header=None, names=['file'])
adsp_ehdn_out_df['SampleId'] = adsp_ehdn_out_df['file'].str.split('_vcpa1').apply(lambda x: x[0])

adsp_cc_df = pd.merge(adsp_case_ctrl_df, adsp_ehdn_out_df, on='SampleId')

adsp_cc_df['Status'].replace('Case', 'case', inplace=True)
adsp_cc_df['Status'].replace('Control', 'control', inplace=True)

adsp_cc_df['file'] = '/sc/arion/projects/sharpa01a/William/ADSP/EHdn/' + adsp_cc_df['file']

###############################################################################

adsp_cc_manifest = adsp_cc_df[['SampleId', 'Status', 'file']]

adsp_cc_manifest.to_csv('Data/ADSP/ADSP_150bp_cc_manifest.txt', sep='\t', header=False, index=False)

###############################################################################

adsp_cc_label = adsp_cc_df[['SampleId', 'Ancestry', 'SeqCenter', 'Sequencer', 'Study', 'BodySite', 'Sex']]
adsp_cc_label.columns.values[0] = 'sample'

adsp_cc_label.to_csv('Data/ADSP/ADSP_150bp_cc_label.txt', sep='\t', index=False)

###############################################################################
###############################################################################
###############################################################################

adsp_case_ctrl_df = pd.read_excel('Data/ADSP/ADSP_14239WGS_FiltCaseControlSamplesWithAnsPCsSequencingInfo_4August21.xlsx')

adsp_case_ctrl_df.dropna(subset=['Status'], inplace=True)
adsp_case_ctrl_df = adsp_case_ctrl_df[adsp_case_ctrl_df['ReadLength']==150]
adsp_case_ctrl_df = adsp_case_ctrl_df[~adsp_case_ctrl_df['Ancestry'].isin(['EAS', 'SAS'])]

rm_outliers = pd.read_csv('Data/ADSP/s22_to_rm.tsv', sep='\t')['sample'].to_list()
adsp_case_ctrl_df = adsp_case_ctrl_df[~adsp_case_ctrl_df['SampleId'].isin(rm_outliers)]

adsp_hiseqx_df = adsp_case_ctrl_df[adsp_case_ctrl_df['Sequencer']=='HiSeqX']
adsp_novaseq_df = adsp_case_ctrl_df[adsp_case_ctrl_df['Sequencer']=='Novaseq']

adsp_ehdn_out_df = pd.read_csv('Data/ADSP/ADSP_EHdn_out.txt', sep='\t', header=None, names=['file'])
adsp_ehdn_out_df['SampleId'] = adsp_ehdn_out_df['file'].str.split('_vcpa1').apply(lambda x: x[0])

adsp_hiseqx_df = pd.merge(adsp_hiseqx_df, adsp_ehdn_out_df, on='SampleId')
adsp_hiseqx_df['Status'].replace(['Case','Control'], ['case','control'], inplace=True)
adsp_hiseqx_df['file'] = '/sc/arion/projects/sharpa01a/William/ADSP/EHdn/' + adsp_hiseqx_df['file']

adsp_novaseq_df = pd.merge(adsp_novaseq_df, adsp_ehdn_out_df, on='SampleId')
adsp_novaseq_df['Status'].replace(['Case','Control'], ['case','control'], inplace=True)
adsp_novaseq_df['file'] = '/sc/arion/projects/sharpa01a/William/ADSP/EHdn/' + adsp_novaseq_df['file']

###############################################################################

adsp_hiseqx_manifest = adsp_hiseqx_df[['SampleId', 'Status', 'file']]
adsp_hiseqx_manifest.to_csv('Data/ADSP/ADSP_hiseqx_manifest.txt', sep='\t', header=False, index=False)

adsp_novaseq_manifest = adsp_novaseq_df[['SampleId', 'Status', 'file']]
adsp_novaseq_manifest.to_csv('Data/ADSP/ADSP_novaseq_manifest.txt', sep='\t', header=False, index=False)

###############################################################################

adsp_hiseqx_label = adsp_hiseqx_df[['SampleId', 'Ancestry', 'SeqCenter', 'Sequencer', 'Study', 'BodySite', 'Sex']]
adsp_hiseqx_label.columns.values[0] = 'sample'
adsp_hiseqx_label.to_csv('Data/ADSP/ADSP_hiseqx_label.txt', sep='\t', index=False)

adsp_novaseq_label = adsp_novaseq_df[['SampleId', 'Ancestry', 'SeqCenter', 'Sequencer', 'Study', 'BodySite', 'Sex']]
adsp_novaseq_label.columns.values[0] = 'sample'
adsp_novaseq_label.to_csv('Data/ADSP/ADSP_novaseq_label.txt', sep='\t', index=False)