#!/bin/bash

# store ONT fastq directory in "ont_fastq"
ont_fastq=$1

# change directory
wd=/sc/arion/projects/sharpa01a/William/
cd "$wd"

# store reference genome in "hg38"
hg38=/sc/arion/projects/sharpa01a/Genomes/Human/RepeatExpansionB38AltContigs/ReferenceGenome/Homo_sapiens_assembly38.fasta

# store unique ID of ont_fastq in "uniq_id"
uniq_id=$(basename "$ont_fastq")

###################################################################################

mkdir "$uniq_id.minimap2_sam"

module load minimap2

for fastq in $ont_fastq/*
do
	minimap2 -ax map-ont -t 48 -2 --MD $hg38 $fastq > "$uniq_id.minimap2_sam"/"$(basename $fastq .fastq).sam"
done

###################################################################################

module load samtools

mkdir "$uniq_id.minimap2_bam"

for sam in "$uniq_id.minimap2_sam"/*
do
	samtools sort -o "$uniq_id.minimap2_bam"/"$(basename $sam .sam).sorted.bam" -@ 12 $sam
done

rm -r "$uniq_id.minimap2_sam"

samtools merge -@ 12 "$uniq_id.minimap2.bam" "$uniq_id.minimap2_bam"/*

###################################################################################

module load sniffles

sniffles -t 12 -m "$uniq_id.minimap2.bam" -v "$uniq_id.minimap2_sniffles.vcf"