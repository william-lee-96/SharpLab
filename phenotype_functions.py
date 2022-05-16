#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  6 16:06:53 2021

@author: williamlee
"""

import pandas as pd
import numpy as np

###############################################################################

def apply_quant_filters(filter_df, complete_df):
    
    filter_df.dropna(subset=['Andy recommends to filter'], inplace=True)
    
    for row in range(len(filter_df)):
    
        var_id = filter_df.iloc[row,9]
        filters = filter_df.iloc[row,2].split()
        
        float_list = []
        compare_list = []
        for element in filters:
            if '>' in element or '<' in element:
                
                float_list.append(float(element.replace(',','')[1:]))
                compare_list.append(element[0])
                
        for i in range(len(float_list)):
            
            if compare_list[i] == '>':
                complete_df[var_id] = complete_df[var_id].apply(lambda x: np.nan if x > float_list[i] else x)
                                  
            elif compare_list[i] == '<':
                complete_df[var_id] = complete_df[var_id].apply(lambda x: np.nan if x < float_list[i] else x)
    
    return complete_df

###############################################################################

def Levenshtein_distance(string1, string2):

    np_array = np.zeros((len(string1)+1, len(string2)+1))
    
    for i in range(1, len(string1)+1):
        np_array[i, 0] = i
     
    for j in range(1, len(string2)+1):
        np_array[0, j] = j
     
    for j in range(1, len(string2)+1):
        for i in range(1, len(string1)+1): 
            if string1[i-1] == string2[j-1]:
              sub_cost = 0
            else:
              sub_cost = 1
    
            np_array[i, j] = min(np_array[i-1, j] + 1,                  
                                 np_array[i, j-1] + 1,                   
                                 np_array[i-1, j-1] + sub_cost)
            
    return np_array[len(string1), len(string2)]

###############################################################################

def match_within_cohorts(phenotype_var_df):
    
    phenotype_var_df.sort_values(by=['Variable name'], ascending=[True], inplace=True)
    phenotype_var_df.reset_index(inplace=True, drop=True)
    phenotype_var_df['Variable name'] = phenotype_var_df['Variable name'].str.lower()
    phenotype_var_df['Variable description'] = phenotype_var_df['Variable description'].str.lower()
    phenotype_var_df['Phenotype group'] = np.nan
    i = 0
    for row in range(len(phenotype_var_df)):
        if pd.isna(phenotype_var_df.loc[row,'Phenotype group']):
            
            var_name = phenotype_var_df.loc[row,'Variable name']
            var_description = phenotype_var_df.loc[row,'Variable description']
            dataset_accession = phenotype_var_df.loc[row,'Dataset accession']
            
            phenotype_var_df.loc[row,'Phenotype group'] = i
            
            nan_df = phenotype_var_df[pd.isna(phenotype_var_df['Phenotype group'])]
            for row_n in range(len(nan_df)):
                if dataset_accession != nan_df.iloc[row_n,3]:
                    
                    if (((var_description in nan_df.iloc[row_n,2] or
                          nan_df.iloc[row_n,2] in var_description) and
                          Levenshtein_distance(var_description, nan_df.iloc[row_n,2]) < 7) or
                         (var_name in nan_df.iloc[row_n,1] and
                          Levenshtein_distance(var_name, nan_df.iloc[row_n,1]) < 3)):
                        
                        nan_df.iloc[row_n,9] = i
                    
                    else: continue              
                else: continue
            
            idx_df = nan_df[pd.notna(nan_df['Phenotype group'])]
            idx_list = idx_df.index.tolist()
            
            phenotype_var_df.loc[idx_list,'Phenotype group'] = i
            
            i += 1
            
    phenotype_var_df.sort_values(by=['Phenotype group'], ascending=[True], inplace=True)
    phenotype_var_df.reset_index(inplace=True, drop=True)
    
    return phenotype_var_df

###############################################################################

def average_quant_vars(phenotype_group_df, complete_df, quant_var_list):
    
    num_pgs = int(max(phenotype_group_df['Phenotype group']))
    
    quant_dict = {}
    avg_quant_df = pd.DataFrame(complete_df['WGS IDs'])
    for pg in range(num_pgs+1):
    
        pg_var_ids = phenotype_group_df[phenotype_group_df['Phenotype group']==pg]['Variable unique'].tolist()
        if all(x in quant_var_list for x in pg_var_ids): # check if all vars are quant
        
            if len(pg_var_ids) == 1:
                avg_quant_df = pd.concat([avg_quant_df, complete_df[pg_var_ids]], axis=1)
            
            else:
                cohort_pg_df = complete_df[pg_var_ids]
                temp_list = []
                for row in range(len(cohort_pg_df)):
                    
                    temp = cohort_pg_df.loc[row, :].tolist()
                    temp = [x for x in temp if x == x]
                    
                    if not temp:
                        temp_list.append(np.nan)
                    else:
                        temp_list.append(sum(temp)/len(temp))
                        
                temp_series = pd.Series(temp_list, name='{}_mean'.format(min(pg_var_ids, key=len).split(':')[0]))
                avg_quant_df = pd.concat([avg_quant_df, temp_series], axis=1)
                
                quant_dict[temp_series.name] = pg_var_ids

    return avg_quant_df, quant_dict

###############################################################################

def combine_cat_vars(phenotype_group_df, complete_df, quant_var_list):
    
    num_pgs = int(max(phenotype_group_df['Phenotype group']))
    
    cat_dict = {}
    combined_cat_df = pd.DataFrame(complete_df['WGS IDs'])
    for pg in range(num_pgs+1):
    
        pg_var_ids = phenotype_group_df[phenotype_group_df['Phenotype group']==pg]['Variable unique'].tolist()
        if all(x not in quant_var_list for x in pg_var_ids): # check if all vars are categorical
        
            if len(pg_var_ids) == 1:
                combined_cat_df = pd.concat([combined_cat_df, complete_df[pg_var_ids]], axis=1)
            
            else:
                
                cat_pg_dict = {}
                for var_id in pg_var_ids:
                    cat_pg_dict[var_id] = complete_df[var_id].unique().tolist()
                    
                num_cat = min([len(v) for k,v in cat_pg_dict.items()])
                cat_pg_dict = {k:v for k,v in cat_pg_dict.items() if len(v) == num_cat}
                
                pg_var_ids = list(cat_pg_dict.keys()) # update this list post-cleaning
                
                if len(pg_var_ids) == 1:
                    combined_cat_df = pd.concat([combined_cat_df, complete_df[pg_var_ids]], axis=1)
                
                else:
                
                    cat_uniq_list = []
                    for k in cat_pg_dict:
                        cat_uniq_list = cat_uniq_list + cat_pg_dict[k]
            
                    cat_uniq_list = set(cat_uniq_list)
                    cat_uniq_list = [x for x in cat_uniq_list if x == x]
                                 
                    cat_pg_df = pd.DataFrame(index=cat_uniq_list)
                    for k in cat_pg_dict:
                        cat_pg_df = cat_pg_df.join(complete_df[k].value_counts())
                        
                    cat_pg_df["sum"] = cat_pg_df.sum(axis=1)
                    cat_pg_df.sort_values(by=['sum'], ascending=[True], inplace=True)
                    
                    cohort_pg_df = complete_df[pg_var_ids]
                    
                    temp_list = []
                    for row in range(len(cohort_pg_df)):
                        
                        temp = cohort_pg_df.loc[row, :].tolist()
                        temp = [x for x in temp if x == x]
                        
                        if not temp:
                            temp_list.append(np.nan)
                        else:
                            for cat in cat_pg_df.index.tolist():
                                if cat in temp:
                                    temp_list.append(cat)
                                    break # go in idx order, once match is found exit loop
                    
                    temp_series = pd.Series(temp_list, name='{}_combined'.format(min(pg_var_ids, key=len).split(':')[0]))
                    combined_cat_df = pd.concat([combined_cat_df, temp_series], axis=1)
        
                    cat_dict[temp_series.name] = pg_var_ids
        
    return combined_cat_df, cat_dict

###############################################################################

def calc_summary_stats(quant_var_df, cohort):
    
    temp_var_list = []
    temp_mean_list = []
    temp_stdv_list = []
    for column in quant_var_df.iloc[:,1:]:
        temp_var_list.append(column)
        temp_mean_list.append(quant_var_df[column].mean())
        temp_stdv_list.append(quant_var_df[column].std())
        
    temp_df = pd.DataFrame(list(zip(temp_var_list, temp_mean_list, temp_stdv_list)), columns =['Variable name','mean','stdv'])
    temp_df['cohort'] = cohort
    
    return temp_df

###############################################################################

def summarize_cat_vars(cat_var_df, cohort):
    
    temp_var_list = []
    temp_cat_list = []
    temp_len_list = []
    for column in cat_var_df.iloc[:,1:]:
        
        temp_var_list.append(column)
        
        cat_list_rm_na = cat_var_df[column].unique().tolist()
        cat_list_rm_na = [x for x in cat_list_rm_na if x == x]
        
        temp_cat_list.append(cat_list_rm_na)
        temp_len_list.append(len(cat_list_rm_na))
        
    temp_df = pd.DataFrame(list(zip(temp_var_list, temp_cat_list, temp_len_list)), columns =['Variable name','categories','# categories'])
    temp_df['cohort'] = cohort
    
    return temp_df