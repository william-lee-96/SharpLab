#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 10 10:52:59 2021

@author: williamlee
"""

import pandas as pd
import numpy as np
import os

wd = '/sc/arion/projects/sharpa01a/William/'
os.chdir(wd)

os.system('mkdir METAL_PLINKs')

ancestry_list = ['AFR','AMR','EUR']
hwe_dict = {'AFR':1e-12,'AMR':1e-8,'EUR':1e-50}

for directory in os.scandir('PLINKs_merged/'):
        os.system(f'mkdir METAL_PLINKs/{directory.name}')
        for file in os.scandir(f'PLINKs_merged/{directory.name}/'):
            if file.name.endswith('.bed'):
                
                prefix = file.name.replace('.bed','')
                contig = prefix.split('_')[0]
                
                prefix_f1 = prefix + '.filter_1'
                                
                os.system(f'plink --bfile PLINKs_merged/{directory.name}/{prefix} --hardy --out PLINKs_merged/{directory.name}/{prefix}')
                os.system(f'plink --bfile PLINKs_merged/{directory.name}/{prefix} --hwe {hwe_dict[directory.name]} --make-bed --out METAL_PLINKs/{directory.name}/{prefix_f1}')
                
                prefix_f2 = prefix + '.filter_2'

                os.system(f'plink --bfile METAL_PLINKs/{directory.name}/{prefix_f1} --exclude snv_ids_exclude/{contig}_snplist.txt --make-bed --out METAL_PLINKs/{directory.name}/{prefix_f2}')
                os.system(f'rm METAL_PLINKs/{directory.name}/{prefix_f1}*')
                
                prefix_f3 = prefix + '.filter_3'
                                
                os.system(f'plink --bfile METAL_PLINKs/{directory.name}/{prefix_f2} --exclude snv_0AC_by_chr/{contig}_snp_0AC.txt --make-bed --out METAL_PLINKs/{directory.name}/{prefix_f3}')
                os.system(f'rm METAL_PLINKs/{directory.name}/{prefix_f2}*')
                
                os.system(f'plink --bfile METAL_PLINKs/{directory.name}/{prefix_f3} --remove GWAS_samplesToRemove.tsv --make-pheno GWAS_case_ctrl_info.tsv "*" --make-bed --out METAL_PLINKs/{directory.name}/{prefix}')
                os.system(f'rm METAL_PLINKs/{directory.name}/{prefix_f3}*')