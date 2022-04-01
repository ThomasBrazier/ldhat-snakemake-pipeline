#!/bin/bash
. /local/env/envsnakemake-6.0.5.sh
. /local/env/envsingularity-3.8.5.sh
. /local/env/envconda.sh
snakemake -s Snakefile -p -j 16 --config dataset=${1} bpen=5 --configfile data/${1}/config.yaml --use-conda --use-singularity --nolock --rerun-incomplete
snakemake -s Snakefile -p -j 16 --config dataset=${1} bpen=15 --configfile data/${1}/config.yaml --use-conda --use-singularity --nolock --rerun-incomplete
snakemake -s Snakefile -p -j 16 --config dataset=${1} bpen=25 --configfile data/${1}/config.yaml --use-conda --use-singularity --nolock --rerun-incomplete
