#!/bin/bash
#SBATCH --mail-user=user_name@mail.com
#SBATCH --mail-type=all
#SBATCH --mem=10GB
#SBATCH --cpus-per-task=1
#SBATCH --time=24-60:00:00
#SBATCH --job-name=LD-rec-map
#SBATCH --output=log/slurm-%A.out

# Load your environment (cluster specific, change accordingly)
. /local/env/envsnakemake-6.0.5.sh
. /local/env/envsingularity-3.8.5.sh
. /local/env/envconda.sh

# Print variables
dataset=${1}
chrom=${2}

# Launch snakemake
echo "Run the pipeline"

snakemake -s workflow/Snakefile -p -j $ncores --configfile data/${dataset}/config.yaml \
    --profile profiles/slurm \
    --use-conda --rerun-incomplete --printshellcmds \
    --config dataset=${dataset} chrom=${chrom}
