#!/bin/bash

CONDA='/home/USER/anaconda3'

snakemake -s "./Snakefile_Nirvana" \
        --cores 35 \
        --resources mem_gb=80 \
        --use-conda \
      	--conda-prefix $CONDA
