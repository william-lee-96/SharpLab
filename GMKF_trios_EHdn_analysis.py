#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 10 17:23:05 2021

@author: williamlee
"""

import pandas as pd 
import numpy as np
import json
import os

# set to a whole number between 0 and 5 (0 = European_OFC_trios, 5 = DSD_trios)
i = 5

mstr_json = ['European_OFC_trios', 'Latin_American_OFC_trios', 'African_Asian_OFC_trios', 'CCDD_trios', 'CDH_trios', 'DSD_trios']
label_tsv = ['European_OFC_labels', 'Latin_American_OFC_labels', 'African_Asian_OFC_labels', 'CCDD_labels', 'CDH_labels', 'DSD_labels']

wd = "/Users/williamlee/SharpLab"
os.chdir(wd)

with open('Data/EHdn-output/{}.multisample_profile.json'.format(mstr_json[i])) as f:
    EHdn_json_dict = json.load(f)
    
mstr_counts_df = pd.DataFrame.from_dict(EHdn_json_dict['Counts'], orient='index')
mstr_counts_df.reset_index(inplace=True)
mstr_counts_df.columns.values[0] = 'motif'
mstr_counts_df = mstr_counts_df[mstr_counts_df['RegionsWithIrrAnchors'].notna()]
mstr_counts_df.reset_index(drop=True, inplace=True)

mstr_parameters_df = pd.DataFrame.from_dict(EHdn_json_dict['Parameters'])
mstr_parameters_df.index.names = ['sample_id']
mstr_parameters_df['norm_factor'] = 40/(mstr_parameters_df['Depths'])

label_df = pd.read_csv('Data/EHdn-output/{}.tsv'.format(label_tsv[i]), sep='\t', header=None)
label_df.columns = ['sample_id', 'status', 'family_id']
label_df.set_index('sample_id', inplace=True)

human_chr = ('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12',
             'chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22')

def split_locus(locus):
    str_list_1 = locus.split(':')
    contig = str_list_1[0]
    start_end = str_list_1[1]
    str_list_2 = start_end.split('-')
    start = str_list_2[0]
    end = str_list_2[1]
    return [contig, start, end]

total_case = len(label_df[label_df['status']=='case'])
total_control = len(label_df[label_df['status']=='control'])
results_list = []

for row in range(len(mstr_counts_df)):
    
    motif = mstr_counts_df.loc[row,'motif']
    region_dict = mstr_counts_df.loc[row,'RegionsWithIrrAnchors']
    
    for region, sample_dict in region_dict.items():
        
        locus = split_locus(region)
        
        if locus[0] not in human_chr: 
            continue # skip to next loop if contig is not chr1-22
        
        else:
                            
            sample_counts = pd.DataFrame(zip(sample_dict.keys(), sample_dict.values()), columns = ['sample_id', 'counts'])
            sample_counts.set_index('sample_id', inplace=True)
            
            sample_counts_merge = sample_counts.join(label_df, how='left')
            sample_counts_merge = sample_counts_merge.join(mstr_parameters_df, how='left')
            sample_counts_merge['norm_counts'] = round(sample_counts_merge['counts']*sample_counts_merge['norm_factor'], 2)
            
            sample_counts_merge_case = sample_counts_merge[sample_counts_merge['status']=='case']
            num_case_pos = len(sample_counts_merge_case)
            
            if num_case_pos == 0:
                continue # skip to next loop if there are 0 case counts
            
            else:
                
                repeat_region = region + ':' + motif
                
                sample_counts_merge_control = sample_counts_merge[sample_counts_merge['status']=='control']
                num_control_pos = len(sample_counts_merge_control)
                
                if len(sample_counts_merge_control) == 0:
                                        
                    trios_list = [sample_counts_merge_case['norm_counts'].tolist(), [0]*num_case_pos, [0]*num_case_pos,
                                  sample_counts_merge_case['family_id'].tolist(), [num_case_pos]*num_case_pos, [1]*num_case_pos]
                    
                    trios_list = list(map(list, zip(*trios_list)))
                    
                else:
                    
                    num_flagged_trios = 0
                    trios_list = []
                    for fam_id in sample_counts_merge_case['family_id']:
                        
                        case_norm_count = sample_counts_merge_case[sample_counts_merge_case['family_id']==fam_id]['norm_counts'][0]
                        
                        control_fam_sub = sample_counts_merge_control[sample_counts_merge_control['family_id']==fam_id]           
                        
                        if len(control_fam_sub) == 0:
                            num_flagged_trios += 1
                            temp_list = [case_norm_count, 0, 0, fam_id]
                            trios_list.append(temp_list)
                            
                        elif len(control_fam_sub) == 1:
                            if control_fam_sub['norm_counts'][0] < case_norm_count:
                                num_flagged_trios += 1
                            temp_list = [case_norm_count, control_fam_sub['norm_counts'][0], 0, fam_id]
                            trios_list.append(temp_list)
                            
                        else:
                            if (control_fam_sub['norm_counts'][0] < case_norm_count) and (control_fam_sub['norm_counts'][1] < case_norm_count):
                                num_flagged_trios += 1
                            temp_list = [case_norm_count, control_fam_sub['norm_counts'][0], control_fam_sub['norm_counts'][1], fam_id]
                            trios_list.append(temp_list)
                    
                    nft_proportion = num_flagged_trios/len(sample_counts_merge['family_id'].unique())
                    [l.extend([num_flagged_trios, nft_proportion]) for l in trios_list]
                            
                current_rows = [l1 + l2 for l1, l2 in zip([[locus[0], locus[1], locus[2], motif, repeat_region, num_case_pos, num_control_pos]]*len(trios_list), trios_list)]
                                
                results_list.extend(current_rows)
                
results_df = pd.DataFrame(results_list) # 486597, 1130624, 627203, 406109, NA, 175050

results_df.columns = ['contig', 'start', 'end', 'motif', 'repeat_region', 'num_case_w_count', 'num_ctrl_w_count', 'child_norm_count', 'parent1_norm_count', 'parent2_norm_count', 'family_id', 'num_flagged_trios', 'prop_flagged_trios']

results_df.sort_values(by=['prop_flagged_trios','num_flagged_trios'], ascending=[False,False], inplace=True)

ref = pd.Series(['0']*len(results_df))
results_df.insert(3, 'ref', ref)

alt = pd.Series(['0']*len(results_df))
results_df.insert(4, 'alt', alt)

results_df.iloc[:,0:5].to_csv('out/to_be_annotated.txt', sep='\t', index=False)

### call ANNOVAR on Minvera; scp '{mstr_json[i]}.variant_function' to ~/SharpLab/out ###
# $ scp ~/SharpLab/out/to_be_annotated.txt leew13@minerva.hpc.mssm.edu:/sc/arion/projects/sharpa01a/William
# $ annovar/annotate_variation.pl -out out/{mstr_json[i]} -build hg38 to_be_annotated.txt humandb/
# $ scp leew13@minerva.hpc.mssm.edu:/sc/arion/projects/sharpa01a/William/out/{mstr_json[i]}.variant_function ~/SharpLab/out

annotations = pd.read_csv('out/{}.variant_function'.format(mstr_json[i]), sep='\t', header=None).iloc[:,0:2] 
annotations.columns =  ['region', 'gene']

results_df.drop(columns = ['ref','alt'], inplace=True) 
results_df.reset_index(inplace=True, drop=True)

gene_annotations = annotations['gene']
region_annotations = annotations['region']

results_df.insert(3, 'gene', gene_annotations)
results_df.insert(4, 'region', region_annotations)

results_df.to_excel('out/{}_EHdn.xlsx'.format(mstr_json[i]), index=False)

###############################################################################
###############################################################################
###############################################################################

# for LA, sheet is too large and must be trimmed - save full result in TSV:
results_df.to_csv('out/{}_EHdn.tsv'.format(mstr_json[i]), sep='\t', index=False)

# and output trimmed sheet here:
latin_american_df = results_df.loc[0:1048575]
latin_american_df.to_excel('out/{}_EHdn.xlsx'.format(mstr_json[i]), index=False)