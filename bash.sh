#!/bin/bash

snakemake -s Snakefile -p -j 16 --config dataset=${1} --configfile data/${1}/config.yaml --use-conda --use-singularity --nolock --rerun-incomplete
