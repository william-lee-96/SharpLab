#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  2 15:41:40 2021

@author: williamlee
"""

import os
import pandas as pd

# set wd
wd = "/Users/williamlee/SharpLab"
os.chdir(wd)

###############################################################################

european_family_df = pd.read_csv("Data/Cavatica/OFC/European/sd_9pyzahhe_20210225_0.tsv", sep='\t')

european_file_df = pd.read_csv("Data/Cavatica/OFC/European/sd_9pyzahhe_20210225_1.tsv", sep='\t')
european_file_df.drop(columns = 'participant_affectedstatus_1', inplace=True)

european_complete_df = pd.merge(european_family_df, european_file_df, on='participant')
european_complete_df.set_index('family', inplace=True)

european_dict = {} 
for fam_id in european_complete_df.index.unique():
    european_dict[fam_id] = european_complete_df.loc[fam_id,:].values.tolist()

# 470 pre-removal
european_dict = {k:v for k,v in european_dict.items() if len(v) == 3 and sum([row[1] for row in european_dict[k]]) == 1}
# 380 post-removal
# 336 post-removal of trios w/ affected parent(s) (1 parent: 40 | 2 parents: 4)

european_file_list = []
for k,v in european_dict.items():
    for row in v:
        european_file_list.append(row[2])
        
pd.Series(european_file_list).to_csv('out/OFC_euro_TT_crams.txt', header=False, index=False) #*# this line added on 10/25/2021

european_file_series = pd.Series(european_file_list, name='file')
european_complete_df.reset_index(inplace=True)

european_trios_df = pd.merge(european_complete_df, european_file_series, on='file')

european_filenames = european_trios_df['file'].str.split(pat='.')
sample_id = pd.Series([filename[0] for filename in european_filenames], name='sample_id')

case_control_status = pd.Series(['case' if status==True else 'control' for status in european_trios_df['participant_affectedstatus_1']], name='case_control_status')

absolute_path = '/sc/arion/projects/sharpa01a/William/European_STR_profiles/' + sample_id + '.str_profile.json'
absolute_path.rename('absolute_path', inplace=True)

european_manifest = pd.concat([sample_id, case_control_status, absolute_path], axis=1)
european_manifest.to_csv('European_OFC_manifest.tsv', sep='\t', index=False, header=False)

european_labels = pd.concat([sample_id, case_control_status, european_trios_df['family']], axis=1)
european_labels.to_csv('Data/EHdn-output/European_OFC_labels.tsv', sep='\t', index=False, header=False)

###############################################################################

latin_american_family_df = pd.read_csv("Data/Cavatica/OFC/LatinAmerican/sd_r0eprsgs_20210225_0.tsv", sep='\t')

latin_american_file_df = pd.read_csv("Data/Cavatica/OFC/LatinAmerican/sd_r0eprsgs_20210225_1.tsv", sep='\t')
latin_american_file_df.drop(columns = 'participant_affectedstatus_1', inplace=True)

latin_american_complete_df = pd.merge(latin_american_family_df, latin_american_file_df, on='participant')
latin_american_complete_df.set_index('family', inplace=True)

latin_american_dict = {} 
for fam_id in latin_american_complete_df.index.unique():
    latin_american_dict[fam_id] = latin_american_complete_df.loc[fam_id,:].values.tolist()

# 271 pre-removal
latin_american_dict = {k:v for k,v in latin_american_dict.items() if len(v) == 3 and sum([row[1] for row in latin_american_dict[k]]) == 1}
# 262 post-removal
# 262 post-removal of trios w/ affected parent(s) (1 parent: 0 | 2 parents: 0)

latin_american_file_list = []
for k,v in latin_american_dict.items():
    for row in v:
        latin_american_file_list.append(row[2])

latin_american_file_series = pd.Series(latin_american_file_list, name='file')
latin_american_complete_df.reset_index(inplace=True)

latin_american_trios_df = pd.merge(latin_american_complete_df, latin_american_file_series, on='file')

latin_american_filenames = latin_american_trios_df['file'].str.split(pat='.')
sample_id = pd.Series([filename[0] for filename in latin_american_filenames], name='sample_id')

case_control_status = pd.Series(['case' if status==True else 'control' for status in latin_american_trios_df['participant_affectedstatus_1']], name='case_control_status')

absolute_path = '/sc/arion/projects/sharpa01a/William/Latin_American_STR_profiles/' + sample_id + '.str_profile.json'
absolute_path.rename('absolute_path', inplace=True)

latin_american_manifest = pd.concat([sample_id, case_control_status, absolute_path], axis=1)
latin_american_manifest.to_csv('Latin_American_OFC_manifest.tsv', sep='\t', index=False, header=False)

latin_american_labels = pd.concat([sample_id, case_control_status, latin_american_trios_df['family']], axis=1)
latin_american_labels.to_csv('Data/EHdn-output/Latin_American_OFC_labels.tsv', sep='\t', index=False, header=False)

###############################################################################

african_asian_family_df = pd.read_csv("Data/Cavatica/OFC/AfricanAsian/sd_dk0krwk8_20210225_0.tsv", sep='\t')

african_asian_file_df = pd.read_csv("Data/Cavatica/OFC/AfricanAsian/sd_dk0krwk8_20210225_1.tsv", sep='\t')
african_asian_file_df.drop(columns = 'participant_affectedstatus_1', inplace=True)

african_asian_complete_df = pd.merge(african_asian_family_df, african_asian_file_df, on='participant')
african_asian_complete_df.set_index('family', inplace=True)

african_asian_dict = {} 
for fam_id in african_asian_complete_df.index.unique():
    african_asian_dict[fam_id] = african_asian_complete_df.loc[fam_id,:].values.tolist()

# 244 pre-removal
african_asian_dict = {k:v for k,v in african_asian_dict.items() if len(v) == 3 and sum([row[1] for row in african_asian_dict[k]]) == 1}
# 238 post-removal
# 231 post-removal of trios w/ affected parent(s) (1 parent: 7 | 2 parents: 0)

for key in african_asian_dict.keys():
    print(key)

african_asian_file_list = []
for k,v in african_asian_dict.items():
    for row in v:
        african_asian_file_list.append(row[2])

african_asian_file_series = pd.Series(african_asian_file_list, name='file')
african_asian_complete_df.reset_index(inplace=True)

african_asian_trios_df = pd.merge(african_asian_complete_df, african_asian_file_series, on='file')

african_asian_filenames = african_asian_trios_df['file'].str.split(pat='.')
sample_id = pd.Series([filename[0] for filename in african_asian_filenames], name='sample_id')

case_control_status = pd.Series(['case' if status==True else 'control' for status in african_asian_trios_df['participant_affectedstatus_1']], name='case_control_status')

absolute_path = '/sc/arion/projects/sharpa01a/William/African_Asian_STR_profiles/' + sample_id + '.str_profile.json'
absolute_path.rename('absolute_path', inplace=True)

african_asian_manifest = pd.concat([sample_id, case_control_status, absolute_path], axis=1)
african_asian_manifest.to_csv('African_Asian_OFC_manifest.tsv', sep='\t', index=False, header=False)

african_asian_labels = pd.concat([sample_id, case_control_status, african_asian_trios_df['family']], axis=1)
african_asian_labels.to_csv('Data/EHdn-output/African_Asian_OFC_labels.tsv', sep='\t', index=False, header=False)