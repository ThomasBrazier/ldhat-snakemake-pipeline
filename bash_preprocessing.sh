#!/bin/bash
#SBATCH --mail-user=thomas.brazier@univ-rennes1.fr
#SBATCH --mail-type=all
#SBATCH --mem=60GB
#SBATCH --cpus-per-task=4
#SBATCH --time=25-60:00:00
#SBATCH --job-name=preprocess

# Load env
. /local/env/envsnakemake-6.0.5.sh
. /local/env/envsingularity-3.8.5.sh
. /local/env/envconda.sh

dataset=${1}
ncores=4

echo "Run pipeline"
snakemake -s data_preprocessing.snake -p -j $ncores --configfile data/${dataset}/config.yaml --use-conda --use-singularity --nolock --rerun-incomplete --config dataset=${dataset}
