#!/bin/bash
## the only expected argument is <dataset>

snakemake -s data_preprocessing.snake -p -j 16 --config dataset=${1} --configfile config.yaml --use-conda --use-singularity --nolock --rerun-incomplete
