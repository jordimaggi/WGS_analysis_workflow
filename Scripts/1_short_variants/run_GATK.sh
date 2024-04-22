#!/usr/bin/bash


CONDA='/home/analyst/anaconda3'

snakemake -s "./Snakefile_GATK" \
	--cores 20 \
	--resources mem_gb=50 \
	--use-conda \
	--conda-prefix $CONDA
