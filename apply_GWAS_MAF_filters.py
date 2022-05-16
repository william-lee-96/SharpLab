#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  9 14:35:47 2021

@author: williamlee
"""

import pandas as pd
import numpy as np
import os
import sys

wd = '/sc/arion/projects/sharpa01a/William/'
os.chdir(wd)

cohort = sys.argv[1]

ancestry_list = ['AFR','AMR','EUR']

os.system(f'mkdir GWAS_PLINKs/{cohort}')
for directory in os.scandir(f'/sc/arion/projects/sharpa01a/TOPMed/TOPMed_plinks/{cohort}/'):
    if directory.name in ancestry_list:
        os.system(f'mkdir GWAS_PLINKs/{cohort}/{directory.name}')
        for file in os.scandir(f'/sc/arion/projects/sharpa01a/TOPMed/TOPMed_plinks/{cohort}/{directory.name}/'):
            if file.name.endswith('.bed'):
                
                prefix = file.name.replace('.bed','')
                contig = prefix.split('.')[3]
                
                if contig != 'chrX':
                
                    os.system(f'plink --bfile /sc/arion/projects/sharpa01a/TOPMed/TOPMed_plinks/{cohort}/{directory.name}/{prefix} --snps-only --maf 0.001 --make-bed --out GWAS_PLINKs/{cohort}/{directory.name}/{prefix}')