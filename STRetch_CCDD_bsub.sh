#!/bin/bash

ml stretch/0.4.0

mkdir STRetch_CCDD_out
mkdir STRetch_CCDD_log

cd STRetch_CCDD_out

for file in /sc/arion/scratch/leew13/CCDD_crams/*.cram
do

	file_bn="$(basename $file .cram)"

	mkdir $file_bn
	cd $file_bn

	bsub -P acc_sharpa01a -q premium -n 48 -R span[hosts=1] -W 72:00 -o ../../STRetch_CCDD_log/"$file_bn.stdout" -eo ../../STRetch_CCDD_log/"$file_bn.stderr" \
	bpipe run \
	-p input_regions=/sc/arion/projects/sharpa01a/Genomes/Human/RepeatExpansionB38AltContigs/TandemRepeatCatalog/TRs_MasterTable.hg38.rmdup.1_6bp.0based.forSTRetch.bed \
	-p STR_BED=/sc/arion/projects/sharpa01a/Genomes/Human/RepeatExpansionB38AltContigs/TandemRepeatCatalog/TRs_MasterTable.hg38.rmdup.1_6bp.0based.forSTRetch.bed \
	-p REF=/sc/arion/projects/sharpa01a/Genomes/Human/RepeatExpansionB38AltContigs/STRechReference/hg38.STRdecoys.sorted.fasta \
	-p CRAM_REF=/sc/arion/projects/sharpa01a/Genomes/Human/RepeatExpansionB38AltContigs/ReferenceGenome/Homo_sapiens_assembly38.fasta \
	-n 48 \
	/hpc/packages/minerva-centos7/stretch/0.4.0/STRetch/pipelines/STRetch_wgs_bam_pipeline.groovy $file

	cd ..

done

for file in /sc/arion/scratch/hadelo01/CCDD_Will_files/*.cram
do

	file_bn="$(basename $file .cram)"

	mkdir $file_bn
	cd $file_bn

	bsub -P acc_sharpa01a -q premium -n 48 -R span[hosts=1] -W 72:00 -o ../../STRetch_CCDD_log/"$file_bn.stdout" -eo ../../STRetch_CCDD_log/"$file_bn.stderr" \
	bpipe run \
	-p input_regions=/sc/arion/projects/sharpa01a/Genomes/Human/RepeatExpansionB38AltContigs/TandemRepeatCatalog/TRs_MasterTable.hg38.rmdup.1_6bp.0based.forSTRetch.bed \
	-p STR_BED=/sc/arion/projects/sharpa01a/Genomes/Human/RepeatExpansionB38AltContigs/TandemRepeatCatalog/TRs_MasterTable.hg38.rmdup.1_6bp.0based.forSTRetch.bed \
	-p REF=/sc/arion/projects/sharpa01a/Genomes/Human/RepeatExpansionB38AltContigs/STRechReference/hg38.STRdecoys.sorted.fasta \
	-p CRAM_REF=/sc/arion/projects/sharpa01a/Genomes/Human/RepeatExpansionB38AltContigs/ReferenceGenome/Homo_sapiens_assembly38.fasta \
	-n 48 \
	/hpc/packages/minerva-centos7/stretch/0.4.0/STRetch/pipelines/STRetch_wgs_bam_pipeline.groovy $file

	cd ..

done