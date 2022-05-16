#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  8 14:53:34 2021

@author: williamlee
"""

import pandas as pd
import numpy as np
import os
import sys
import gzip

from statistics import mean, median, stdev

wd = '/sc/arion/projects/sharpa01a/William/'
os.chdir(wd)

cohort = sys.argv[1]
contig = sys.argv[2]

for file in os.scandir(f'/sc/arion/projects/sharpa01a/TOPMed/TOPMed_vcfs/{cohort}/'):
    if file.name.endswith('.vcf.gz') and file.name.split('.')[2] == contig:
        
        avgdp_list = []
        ac_list = []
        af_list = []
        qual_list = []
        
        with gzip.open(file.path, 'rt') as infile:
            num_calls = 0
            for line in infile:
                if line.startswith('#'):
                    continue
                else:
                    line_list = line.split('\t')
                    
                    info = line_list[7].split(';')
                    
                    avgdp_list.append(info[0].split('=')[1])
                    ac_list.append(info[1].split('=')[1])
                    af_list.append(info[3].split('=')[1])
                    
                    qual_list.append(line_list[5])
                    
                    num_calls += 1
        
        avgdp_list = [float(x) for x in avgdp_list]
        ac_list = [int(x) for x in ac_list]
        af_list = [float(x) for x in af_list]
        qual_list = [int(x) for x in qual_list]
                                            
        avgdp_mean = mean(avgdp_list)
        avgdp_median = median(avgdp_list)
        avgdp_std = stdev(avgdp_list)
        
        ac_mean = mean(ac_list)
        ac_median = median(ac_list)
        ac_std = stdev(ac_list)
        
        af_mean = mean(af_list)
        af_median = median(af_list)
        af_std = stdev(af_list)
        
        qual_mean = mean(qual_list)
        qual_median = median(qual_list)
        qual_std = stdev(qual_list)
                
        cur_list = [[cohort,contig,num_calls,avgdp_mean,avgdp_median,avgdp_std,ac_mean,ac_median,ac_std,af_mean,af_median,af_std,qual_mean,qual_median,qual_std]]

        cur_df = pd.DataFrame(cur_list, columns=['cohort','chr','num_calls','AVGDP_mean','AVGDP_median','AVGDP_stdev','AC_mean','AC_median','AC_stdev','AF_mean','AF_median','AF_stdev','QUAL_mean','QUAL_median','QUAL_stdev'])
        
        cur_df.to_csv(f'vcf_metrics_out/{cohort}_{contig}.tsv', sep='\t', index=False)
                    
    else:
        continue