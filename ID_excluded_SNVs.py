#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 22 16:13:41 2021

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
                    
                    avgdp = float(info[0].split('=')[1])
                    abz = float(info[15].split('=')[1])
                    qual = int(line_list[5])
                    
                    is_snv = (len(line_list[3]) == 1 and len(line_list[4]) == 1)
        
                    if (avgdp < 10 or abz < -3 or qual < 20) and is_snv:
                    
                        snv_id = line_list[0].replace('chr','') + ':' + line_list[1] + '_' + line_list[3] + '_' + line_list[4]
                        id_list.append(snv_id)
                    
                    else:
                        continue
    
id_series = pd.Series(id_list)

id_series.to_csv(f'snv_exclude_out/{cohort}_{contig}_snplist.txt', sep='\t', header=False, index=False)