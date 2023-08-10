#!/bin/bash

ml regenie/3.2.2

mkdir UKB_REGENIE_step1_5QP_out
mkdir UKB_REGENIE_step1_5QP_log

for pheno in mutation_size expansion_size contraction_size expansion_propensity contraction_propensity
do
	
	mkdir UKB_REGENIE_step1_5QP_out/$pheno
	mkdir UKB_REGENIE_step1_5QP_log/$pheno

	for cohort in SC deCODE
	do

		mkdir UKB_REGENIE_step1_5QP_out/$pheno/$cohort
		mkdir UKB_REGENIE_step1_5QP_log/$pheno/$cohort

		for file in /sc/arion/projects/sharpa01a/Paras/Projects/sharplab2/UK_biobank/Analysis/SNPs/DNAnexus/WGS_PLINK_filtered_SNPs/merged_files/*.bed
		do

			file_bn=$(basename "$file" .bed)

			bsub -P acc_sharpa01a -q express -n 16 -R span[hosts=1] -W 8:00 -o UKB_REGENIE_step1_5QP_log/$pheno/$cohort/"$file_bn.log" \
			regenie     \
			    --step 1 \
			    --bed /sc/arion/projects/sharpa01a/Paras/Projects/sharplab2/UK_biobank/Analysis/SNPs/DNAnexus/WGS_PLINK_filtered_SNPs/merged_files/$file_bn \
			    --extract REGENIE_input/$cohort/snplists/$file_bn.snplist \
			    --keep REGENIE_input/$cohort/UKB_${cohort}_samples_to_keep.tsv \
			    --bsize 1000 \
			    --threads 16 \
			    --phenoFile REGENIE_input/$cohort/quant_phenos/${pheno}.tsv \
			    --covarFile REGENIE_input/$cohort/UKB_${cohort}_covar.tsv \
			    --apply-rint \
			    --force-step1 \
			    --out UKB_REGENIE_step1_5QP_out/$pheno/$cohort/$file_bn

		done
	done
done