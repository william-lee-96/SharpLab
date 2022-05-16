#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  3 04:11:36 2021

@author: williamlee
"""

import os
import numpy as np
import pandas as pd

# set wd
wd = "/Users/williamlee/SharpLab"
os.chdir(wd)

###############################################################################
###############################################################################
###############################################################################

def true_trio(dict_val):
    
    if len(dict_val) != 3:
        return False
    
    else:
        i = 0
        j = 0
        for l in dict_val:
            
            if [True, 'Parent'] == l[1:3]:
                i += 1
            
            elif [False, 'Child'] == l[1:3]:
                j += 1
        
        if i == 1 and j == 2:
            return True
        else:
            return False

###############################################################################
###############################################################################
###############################################################################

CDH_family_df = pd.read_csv("Data/Cavatica/CDH_trios/sd_46sk55a3_20210504_0.tsv", sep='\t')

CDH_relationship_df = pd.read_csv("Data/Cavatica/CDH_trios/sd_46sk55a3_20210504_1.tsv", sep='\t')
CDH_relationship_df['familyrelationship_relationship_1'].replace(['Father','Mother'], ['Parent','Parent'], inplace=True)
CDH_relationship_df.drop(columns=['participant_affectedstatus_1', 'familyrelationship'], inplace=True)
CDH_relationship_df.drop_duplicates(subset='participant', inplace=True, ignore_index=True)

CDH_merge_df = pd.merge(CDH_family_df, CDH_relationship_df, on='participant')

CDH_file_df = pd.read_csv("Data/Cavatica/CDH_trios/sd_46sk55a3_20210504_2.tsv", sep='\t')
CDH_file_df.drop(columns = 'participant_affectedstatus_1', inplace=True)
CDH_file_df.drop_duplicates(subset='participant', inplace=True, ignore_index=True)

CDH_complete_df = pd.merge(CDH_merge_df, CDH_file_df, on='participant')
CDH_complete_df.set_index('family', inplace=True)

CDH_dict = {} 
for fam_id in CDH_complete_df.index.unique():
    CDH_dict[fam_id] = CDH_complete_df.loc[fam_id,:].values.tolist()
    
# 731 pre-removal
CDH_dict = {k:v for k,v in CDH_dict.items() if true_trio(CDH_dict[k])}
# 485 post-removal

###############################################################################

CDH_file_list = []
{CDH_file_list.extend([l[3] for l in v]) for k,v in CDH_dict.items()}

CDH_file_series = pd.Series(CDH_file_list, name='file')
CDH_complete_df.reset_index(inplace=True)

CDH_trios_df = pd.merge(CDH_complete_df, CDH_file_series, on='file')

CDH_filenames = CDH_trios_df['file'].str.split(pat='.')
sample_id = pd.Series([filename[0] for filename in CDH_filenames], name='sample_id')
STR_names = sample_id + '.str_profile.json'
STR_names.rename('STR_profile', inplace=True)

CDH_trios_df.reset_index(drop=True, inplace=True)

sample_id = pd.Series([filename[0] for filename in CDH_trios_df['file'].str.split(pat='.')], name='sample_id')

case_control_status = pd.Series(['case' if status==True else 'control' for status in CDH_trios_df['participant_affectedstatus_1']], name='case_control_status')

absolute_path = '/sc/arion/projects/sharpa01a/William/CDH_STR_profiles/' + sample_id + '.str_profile.json'
absolute_path.rename('absolute_path', inplace=True)

CDH_manifest = pd.concat([sample_id, case_control_status, absolute_path], axis=1)

CDH_manifest.to_csv('CDH_manifest.tsv', sep='\t', index=False, header=False)

CDH_labels = pd.concat([sample_id, case_control_status, CDH_trios_df['family']], axis=1)
CDH_labels.to_csv('Data/EHdn-output/CDH_labels.tsv', sep='\t', index=False, header=False)

###############################################################################
###############################################################################
###############################################################################

CCDD_family_df = pd.read_csv("Data/Cavatica/CCDD_trios/sd_dztb5hrr_20210505_0.tsv", sep='\t')

CCDD_diagnosis_df = pd.read_csv("Data/Cavatica/CCDD_diagnosis/sd_dztb5hrr_20210506_0.tsv", sep='\t')
CCDD_diagnosis_df.drop(columns='diagnosis', inplace=True)
CCDD_diagnosis_df['diagnosis_diagnosiscategory_1'] = True

CCDD_merge_df = pd.merge(CCDD_family_df, CCDD_diagnosis_df, on='participant', how='left')
CCDD_merge_df['diagnosis_diagnosiscategory_1'].replace(np.nan, False, inplace=True)

CCDD_relationship_df = pd.read_csv("Data/Cavatica/CCDD_trios/sd_dztb5hrr_20210505_1.tsv", sep='\t')
CCDD_relationship_df['familyrelationship_relationship_1'].replace(['Father','Mother','Son','Daughter'], ['Parent','Parent','Child','Child'], inplace=True)
CCDD_relationship_df.drop(columns='familyrelationship', inplace=True)
CCDD_relationship_df.drop_duplicates(subset='participant', inplace=True, ignore_index=True)

CCDD_merge_df = pd.merge(CCDD_merge_df, CCDD_relationship_df, on='participant')

CCDD_file_df = pd.read_csv("Data/Cavatica/CCDD_trios/sd_dztb5hrr_20210505_2.tsv", sep='\t')

CCDD_complete_df = pd.merge(CCDD_merge_df, CCDD_file_df, on='participant')
CCDD_complete_df.set_index('family', inplace=True)

CCDD_dict = {} 
for fam_id in CCDD_complete_df.index.unique():
    CCDD_dict[fam_id] = CCDD_complete_df.loc[fam_id,:].values.tolist()
    
# 248 pre-removal
CCDD_dict = {k:v for k,v in CCDD_dict.items() if true_trio(CCDD_dict[k])}
# 196 post-removal

###############################################################################

CCDD_file_list = []
for k,v in CCDD_dict.items():
    for row in v:
        CCDD_file_list.append(row[3])
        
CCDD_file_series = pd.Series(CCDD_file_list, name='file')
CCDD_complete_df.reset_index(inplace=True)

CCDD_trios_df = pd.merge(CCDD_complete_df, CCDD_file_series, on='file')

CCDD_filenames = CCDD_trios_df['file'].str.split(pat='.')
sample_id = pd.Series([filename[0] for filename in CCDD_filenames], name='sample_id')
STR_names = sample_id + '.str_profile.json'
STR_names.rename('STR_profile', inplace=True)

CCDD_trios_df.reset_index(drop=True, inplace=True)

sample_id = pd.Series([filename[0] for filename in CCDD_trios_df['file'].str.split(pat='.')], name='sample_id')

case_control_status = pd.Series(['case' if status==True else 'control' for status in CCDD_trios_df['diagnosis_diagnosiscategory_1']], name='case_control_status')

absolute_path = '/sc/arion/projects/sharpa01a/William/CCDD_STR_profiles/' + sample_id + '.str_profile.json'
absolute_path.rename('absolute_path', inplace=True)

CCDD_manifest = pd.concat([sample_id, case_control_status, absolute_path], axis=1)

CCDD_manifest.to_csv('CCDD_manifest.tsv', sep='\t', index=False, header=False)

CCDD_labels = pd.concat([sample_id, case_control_status, CCDD_trios_df['family']], axis=1)
CCDD_labels.to_csv('Data/EHdn-output/CCDD_labels.tsv', sep='\t', index=False, header=False)

###############################################################################
###############################################################################
###############################################################################

DSD_family_df = pd.read_csv("Data/Cavatica/DSD_trios/sd_6fpyjqbr_20210505_0.tsv", sep='\t')

DSD_file_df = pd.read_csv("Data/Cavatica/DSD_trios/sd_6fpyjqbr_20210505_1.tsv", sep='\t')
DSD_file_df.drop(columns='participant_affectedstatus_1', inplace=True)

DSD_merge_df = pd.merge(DSD_family_df, DSD_file_df, on='participant')

DSD_relationship_df = pd.read_csv("Data/Cavatica/DSD_trios/sd_6fpyjqbr_20210505_2.tsv", sep='\t')
DSD_relationship_df['familyrelationship_relationship_1'].replace('Father', 'Parent', inplace=True)
DSD_relationship_df['familyrelationship_relationship_1'].replace('Mother', 'Parent', inplace=True)
DSD_relationship_df.drop(columns=['participant_affectedstatus_1', 'familyrelationship'], inplace=True)
DSD_relationship_df.drop_duplicates(subset='participant', inplace=True, ignore_index=True)

DSD_merge_df = pd.merge(DSD_family_df, DSD_relationship_df, on='participant')

DSD_file_df = pd.read_csv("Data/Cavatica/DSD_trios/sd_6fpyjqbr_20210505_1.tsv", sep='\t')
DSD_file_df.drop(columns = 'participant_affectedstatus_1', inplace=True)

DSD_complete_df = pd.merge(DSD_merge_df, DSD_file_df, on='participant')
DSD_complete_df.set_index('family', inplace=True)

DSD_dict = {} 
for fam_id in DSD_complete_df.index.unique():
    DSD_dict[fam_id] = DSD_complete_df.loc[fam_id,:].values.tolist()
    
# 92 pre-removal
DSD_dict = {k:v for k,v in DSD_dict.items() if true_trio(DSD_dict[k])}
# 79 post-removal

###############################################################################

DSD_file_list = []
for k,v in DSD_dict.items():
    for row in v:
        DSD_file_list.append(row[3])
        
DSD_file_series = pd.Series(DSD_file_list, name='file')
DSD_complete_df.reset_index(inplace=True)

DSD_trios_df = pd.merge(DSD_complete_df, DSD_file_series, on='file')

DSD_trios_df.iloc[:,0:5].to_excel('out/DSD_trios_forGaby.xlsx', index=False) #*# this line added on 10/19/2021

DSD_filenames = DSD_trios_df['file'].str.split(pat='.')
sample_id = pd.Series([filename[0] for filename in DSD_filenames], name='sample_id')
STR_names = sample_id + '.str_profile.json'
STR_names.rename('STR_profile', inplace=True)

DSD_trios_df.reset_index(drop=True, inplace=True)

sample_id = pd.Series([filename[0] for filename in DSD_trios_df['file'].str.split(pat='.')], name='sample_id')

case_control_status = pd.Series(['case' if status==True else 'control' for status in DSD_trios_df['participant_affectedstatus_1']], name='case_control_status')

absolute_path = '/sc/arion/projects/sharpa01a/William/DSD_STR_profiles/' + sample_id + '.str_profile.json'
absolute_path.rename('absolute_path', inplace=True)

DSD_manifest = pd.concat([sample_id, case_control_status, absolute_path], axis=1)

DSD_manifest.to_csv('DSD_manifest.tsv', sep='\t', index=False, header=False)

DSD_labels = pd.concat([sample_id, case_control_status, DSD_trios_df['family']], axis=1)
DSD_labels.to_csv('Data/EHdn-output/DSD_labels.tsv', sep='\t', index=False, header=False)