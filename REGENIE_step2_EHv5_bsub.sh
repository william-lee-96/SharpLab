#!/bin/bash

mkdir REGENIE_step2_EHv5_out
mkdir REGENIE_step2_EHv5_log

ml regenie/2.0.2

for threshold in VeryLow Low Medium High VeryHigh
do

	mkdir REGENIE_step2_EHv5_out/$threshold
	mkdir REGENIE_step2_EHv5_log/$threshold

	for dir in AFR AMR EUR
	do

		mkdir REGENIE_step2_EHv5_out/$threshold/$dir
		mkdir REGENIE_step2_EHv5_log/$threshold/$dir

		for cohort in ARIC BioMe.Baylor BioMe.MGI CHS COPD.Broad COPD.UW FHS HyperGen JHS MESA MGH VTE VUAF WHI
		do

			mkdir REGENIE_step2_EHv5_out/$threshold/$dir/$cohort
			mkdir REGENIE_step2_EHv5_log/$threshold/$dir/$cohort

			for file in GWAS_PLINKs/$dir/*.bed
			do

				file_bn=$(basename "$file" .bed)

				bsub -P acc_sharpa01a -q premium -n 48 -R span[hosts=1] -W 24:00 -o REGENIE_step2_EHv5_log/$threshold/$dir/$cohort/"$file_bn.log" \
				regenie     \
				    --step 2 \
				    --bed GWAS_PLINKs/$dir/$file_bn \
				    --extract MAC_filter_out/$dir/$cohort/"$file_bn.snplist" \
				    --keep GWAS_filtered_cohorts/${cohort}_samples_to_keep.tsv \
				    --bsize 100 \
				    --threads 48 \
				    --phenoFile REGENIE_EHv5_pheno/Phewas_GWAS_FinalFiles/mutation_count/REGENIE_pheno_${threshold}.tsv \
				    --covarFile REGENIE_EHv5_covar/Phewas_GWAS_FinalFiles/REGENIE_covar.tsv \
				    --pred REGENIE_step1_EHv5_out/$threshold/$dir/$cohort/"${file_bn}_pred.list" \
				    --out REGENIE_step2_EHv5_out/$threshold/$dir/$cohort/$file_bn

			done
		done
	done
done