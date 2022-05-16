#!/bin/bash

SCHEME STDERR

###############################################################################
###############################################################################
###############################################################################

MARKER ID
ALLELE ALLELE0 ALLELE1
EFFECT BETA
STDERR SE
PROCESS GWAS_for_METAL/Medium/EUR/COPD.UW_merged.txt
PROCESS GWAS_for_METAL/Medium/EUR/BioMe.Baylor_merged.txt
PROCESS GWAS_for_METAL/Medium/EUR/MESA_merged.txt
PROCESS GWAS_for_METAL/Medium/EUR/MGH_merged.txt
PROCESS GWAS_for_METAL/Medium/EUR/COPD.Broad_merged.txt
PROCESS GWAS_for_METAL/Medium/EUR/VUAF_merged.txt
PROCESS GWAS_for_METAL/Medium/EUR/BioMe.MGI_merged.txt
PROCESS GWAS_for_METAL/Medium/EUR/ARIC_merged.txt
PROCESS GWAS_for_METAL/Medium/EUR/WHI_merged.txt
PROCESS GWAS_for_METAL/Medium/EUR/VTE_merged.txt
PROCESS GWAS_for_METAL/Medium/EUR/CHS_merged.txt
PROCESS GWAS_for_METAL/Medium/EUR/FHS_merged.txt

OUTFILE GWAS_EHv5_Medium_EUR .txt
ANALYZE
CLEAR

###############################################################################

MARKER ID
ALLELE ALLELE0 ALLELE1
EFFECT BETA
STDERR SE
PROCESS GWAS_for_METAL/Medium/AFR/COPD.UW_merged.txt
PROCESS GWAS_for_METAL/Medium/AFR/BioMe.Baylor_merged.txt
PROCESS GWAS_for_METAL/Medium/AFR/MESA_merged.txt
PROCESS GWAS_for_METAL/Medium/AFR/JHS_merged.txt
PROCESS GWAS_for_METAL/Medium/AFR/COPD.Broad_merged.txt
PROCESS GWAS_for_METAL/Medium/AFR/BioMe.MGI_merged.txt
PROCESS GWAS_for_METAL/Medium/AFR/ARIC_merged.txt
PROCESS GWAS_for_METAL/Medium/AFR/WHI_merged.txt
PROCESS GWAS_for_METAL/Medium/AFR/CHS_merged.txt
PROCESS GWAS_for_METAL/Medium/AFR/HyperGen_merged.txt

OUTFILE GWAS_EHv5_Medium_AFR .txt
ANALYZE
CLEAR

###############################################################################

MARKER ID
ALLELE ALLELE0 ALLELE1
EFFECT BETA
STDERR SE
PROCESS GWAS_for_METAL/Medium/AMR/BioMe.Baylor_merged.txt
PROCESS GWAS_for_METAL/Medium/AMR/MESA_merged.txt
PROCESS GWAS_for_METAL/Medium/AMR/COPD.Broad_merged.txt
PROCESS GWAS_for_METAL/Medium/AMR/BioMe.MGI_merged.txt
PROCESS GWAS_for_METAL/Medium/AMR/WHI_merged.txt

OUTFILE GWAS_EHv5_Medium_AMR .txt
ANALYZE
CLEAR

###############################################################################
###############################################################################
###############################################################################

MARKER MarkerName
ALLELE Allele1 Allele2
EFFECT Effect
STDERR StdErr
PROCESS GWAS_EHv5_Medium_EUR1.txt
PROCESS GWAS_EHv5_Medium_AFR1.txt
PROCESS GWAS_EHv5_Medium_AMR1.txt

OUTFILE GWAS_EHv5_Medium_meta .txt
ANALYZE
CLEAR