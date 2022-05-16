#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec  5 14:25:22 2021

@author: williamlee
"""

import pandas as pd
import numpy as np
import os
import sys

wd = '/sc/arion/projects/sharpa01a/William/'
os.chdir(wd)

ancestry = sys.argv[1]
contig = sys.argv[2]

# ARIC removed
cohort_list = ['AmishCG','BioMe','CHS.c1c2','COPD','FHS.c1c2','HyperGen','JHS.c1c2','MESA.c1c2','MGH','VTE','VUAF','WHI.c1c2']

for file in os.scandir(f'GWAS_PLINKs/ARIC/{ancestry}/'):
    if file.name.endswith('.bed') and file.name.split('.')[3] == contig:
        
        prefix = file.name.replace('.bed','')
        first_merge = False
        
        for cohort in cohort_list:
            for directory in os.scandir(f'GWAS_PLINKs/{cohort}/'):
                if directory.name == ancestry:
                    for file_2 in os.scandir(f'GWAS_PLINKs/{cohort}/{ancestry}/'):
                        if file_2.name.endswith('.bed'):
                        
                            prefix_2 = file_2.name.replace('.bed','')
                            
                            if prefix_2.split('.')[3] == contig:
                                
                                if first_merge == False:
                                
                                    os.system(f'plink --bfile GWAS_PLINKs/ARIC/{ancestry}/{prefix} --bmerge GWAS_PLINKs/{cohort}/{ancestry}/{prefix_2}.bed GWAS_PLINKs/{cohort}/{ancestry}/{prefix_2}.bim GWAS_PLINKs/{cohort}/{ancestry}/{prefix_2}.fam --make-bed --out PLINKs_merged/{ancestry}/{contig}_merged')
                                    first_merge = True
                                    
                                else:
                                    
                                    os.system(f'plink --bfile PLINKs_merged/{ancestry}/{contig}_merged --bmerge GWAS_PLINKs/{cohort}/{ancestry}/{prefix_2}.bed GWAS_PLINKs/{cohort}/{ancestry}/{prefix_2}.bim GWAS_PLINKs/{cohort}/{ancestry}/{prefix_2}.fam --make-bed --out PLINKs_merged/{ancestry}/{contig}_merged')