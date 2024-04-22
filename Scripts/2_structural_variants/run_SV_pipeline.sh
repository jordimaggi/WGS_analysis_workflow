#!/usr/bin/bash

CONDA='/home/analyst/anaconda3'

snakemake -s "./Snakefile_SV_pipeline" \
        --cores 20 \
        --resources mem_gb=50 \
        --use-conda \
        --rerun-triggers  mtime input \
      	--conda-prefix $CONDA
