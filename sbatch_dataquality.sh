#!/bin/bash
#SBATCH --mail-user=thomas.brazier@univ-rennes1.fr
#SBATCH --mail-type=all
#SBATCH --time=12-60:00:00
## the only expected argument is <dataset>

## load the necessary environement (used for Genouest cluster)
. /local/env/envsnakemake-6.0.5.sh
. /local/env/envsingularity-3.8.5.sh

snakemake -s data_quality.snake -p -j 16 --config dataset=${1} --configfile data/${1}/config.yaml  --cluster-config cluster.yaml --cluster "sbatch --account={cluster.account} --nodes={cluster.node} --ntasks={cluster.n} --cpus-per-task={cluster.cpusPerTask} --time={cluster.time} --mem-per-cpu={cluster.memPerCpu} -p {cluster.p}" --use-conda --use-singularity --nolock --rerun-incomplete -F
