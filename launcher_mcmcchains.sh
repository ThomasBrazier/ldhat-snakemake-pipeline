#!/bin/bash
#SBATCH --mail-user=thomas.brazier@univ-rennes1.fr
#SBATCH --mail-type=all
#SBATCH --mem=180GB
#SBATCH --cpus-per-task=10
#SBATCH --time=25-60:00:00
#SBATCH --job-name=ldhatMCMC

# Load env
. /local/env/envsnakemake-6.0.5.sh
. /local/env/envconda.sh

path="LDhat-snakemake-pipeline"
dataset=${1}
ncores=10

cd $path/ldhat-snakemake-pipeline

echo "Run pipeline"
snakemake -s data_mcmcchains.snake -p -j $ncores --configfile data/${dataset}/config.yaml --use-conda --use-singularity --nolock --rerun-incomplete --config dataset=${dataset}
