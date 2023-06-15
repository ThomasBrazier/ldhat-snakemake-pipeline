#!/bin/bash
#SBATCH --mail-user=thomas.brazier@univ-rennes1.fr
#SBATCH --mail-type=all
#SBATCH --mem=60GB
#SBATCH --time=20-60:00:00
#SBATCH --job-name=LDdataqual

# Load env
. /local/env/envsnakemake-6.0.5.sh
. /local/env/envsingularity-3.8.5.sh
. /local/env/envconda.sh

dataset=${1}
ncores=16

path="LDhat-snakemake-pipeline"
cd $path/ldhat-snakemake-pipeline


echo "Run pipeline"
snakemake -s data_quality.snake -p -j $ncores --configfile data/${dataset}/config.yaml --use-conda --use-singularity --nolock --rerun-incomplete --config dataset=${dataset}

