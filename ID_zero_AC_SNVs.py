#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 24 12:24:33 2021

@author: williamlee
"""

import pandas as pd
import numpy as np
import os
import sys
import gzip

wd = '/sc/arion/projects/sharpa01a/William/'
os.chdir(wd)

cohort = sys.argv[1]
contig = sys.argv[2]

for file in os.scandir(f'/sc/arion/projects/sharpa01a/TOPMed/TOPMed_vcfs/{cohort}/'):
    if file.name.endswith('.vcf.gz') and file.name.split('.')[2] == contig:
        
        id_list = []
        
        with gzip.open(file.path, 'rt') as infile:
            for line in infile:
                
                if line.startswith('#'):
                    continue
                else:
                    
                    line_list = line.split('\t')
                    info = line_list[7].split(';')
                    
                    ac = int(info[1].split('=')[1])
                    
                    is_snv = (len(line_list[3]) == 1 and len(line_list[4]) == 1)
        
                    if ac == 0 and is_snv:
                    
                        snv_id = line_list[0].replace('chr','') + ':' + line_list[1] + '_' + line_list[4] + '_' + line_list[3]
                        id_list.append(snv_id)
                    
                    else:
                        continue
    
id_series = pd.Series(id_list)

id_series.to_csv(f'snv_zero_AC/{cohort}_{contig}_snp_0AC.txt', sep='\t', header=False, index=False)