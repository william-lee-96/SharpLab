#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 24 16:10:55 2023

@author: williamlee
"""

import os
import pandas as pd
import numpy as np

os.chdir('/Users/williamlee/SharpLab')

###############################################################################
###############################################################################

amppd_df = pd.read_csv('../Downloads/AMP_PD_24k_TR_List_24April.tsv', sep='\t')

amppd_df['TRGT_ID'] = amppd_df['chrom'] + '_' + amppd_df['start'].astype('string') + '_' + amppd_df['end'].astype('string') + '_' + amppd_df['RefUnit']

amppd_df['info'] = 'ID=' + amppd_df['TRGT_ID'] + ';'
amppd_df['info'] = amppd_df['info'] + 'MOTIFS=' + amppd_df['RefUnit'] + ';'
amppd_df['info'] = amppd_df['info'] + 'STRUC=(' + amppd_df['RefUnit'] + ')n'

amppd_df = amppd_df[['chrom', 'start', 'end', 'info']]

###############################################################################

pcgc_df = pd.read_json('../Downloads/Combined_TRListFinalList.n6566.json')

pcgc_df[['chrom', 'start', 'end']] = pcgc_df['VariantId'].str.split('_').tolist()

pcgc_df['motif'] = pcgc_df['LocusStructure'].str.replace('(', '', regex=False)
pcgc_df['motif'] = pcgc_df['motif'].str.replace(')*', '', regex=False)

pcgc_df['LocusStructure'] = pcgc_df['LocusStructure'].str.replace('*', 'n', regex=False)

pcgc_df['info'] = 'ID=' + pcgc_df['VariantId'] + '_' + pcgc_df['motif'] + ';'
pcgc_df['info'] = pcgc_df['info'] + 'MOTIFS=' + pcgc_df['motif'] + ';'
pcgc_df['info'] = pcgc_df['info'] + 'STRUC=' + pcgc_df['LocusStructure']

pcgc_df = pcgc_df[['chrom', 'start', 'end', 'info']]

###############################################################################

gmkf_df = pd.read_csv('../Downloads/All_FINAL_FinalList.tsv', sep='\t')

gmkf_df['TRGT_ID'] = gmkf_df['chrom'] + '_' + gmkf_df['start'].astype('string') + '_' + gmkf_df['end'].astype('string') + '_' + gmkf_df['motif']

gmkf_df['info'] = 'ID=' + gmkf_df['TRGT_ID'] + ';'
gmkf_df['info'] = gmkf_df['info'] + 'MOTIFS=' + gmkf_df['motif'] + ';'
gmkf_df['info'] = gmkf_df['info'] + 'STRUC=(' + gmkf_df['motif'] + ')n'

gmkf_df = gmkf_df[['chrom', 'start', 'end', 'info']]

###############################################################################

ehv5_df = pd.read_csv('../Downloads/EHv5_Feb23_TRLIst_ForWill.bed', sep='\t')

ehv5_df['TRGT_ID'] = ehv5_df['chrom'] + '_' + ehv5_df['start'].astype('string') + '_' + ehv5_df['end'].astype('string') + '_' + ehv5_df['RefUnit']

ehv5_df['info'] = 'ID=' + ehv5_df['TRGT_ID'] + ';'
ehv5_df['info'] = ehv5_df['info'] + 'MOTIFS=' + ehv5_df['RefUnit'] + ';'
ehv5_df['info'] = ehv5_df['info'] + 'STRUC=(' + ehv5_df['RefUnit'] + ')n'

ehv5_df = ehv5_df[['chrom', 'start', 'end', 'info']]

###############################################################################
###############################################################################

sharp_df = pd.concat([amppd_df, pcgc_df, gmkf_df, ehv5_df], ignore_index=True)
sharp_df = sharp_df.astype({'chrom':'string', 'start':'int32', 'end':'int32', 'info':'string'})
sharp_df.drop_duplicates(inplace=True, ignore_index=True)

# sharp_df.to_csv('out/SharpLab_TR_catalog.bed', sep='\t', index=False, header=False)

trgt_df = pd.read_csv('out/TRGT_TR_catalog_no_overlap.bed', sep='\t', header=None)
trgt_df.columns = ['chrom', 'start', 'end', 'info']
trgt_df = trgt_df.astype({'chrom':'string', 'start':'int32', 'end':'int32', 'info':'string'})

###############################################################################
###############################################################################

out_df = pd.concat([trgt_df, sharp_df], ignore_index=True)

out_df.to_csv('out/TRGT_complete_catalog.hg38.bed', sep='\t', index=False, header=False)