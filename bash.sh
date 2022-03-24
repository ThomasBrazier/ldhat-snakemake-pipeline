#!/bin/bash
conda activate snakemake

snakemake -s Snakefile -p -j 16 --config dataset=${1} --configfile data/${1}/config.yaml  --cluster-config cluster.yaml --cluster "sbatch --account={cluster.account} --nodes={cluster.node} --ntasks={cluster.n} --cpus-per-task={cluster.cpusPerTask} --time={cluster.time} --mem-per-cpu={cluster.memPerCpu} -p {cluster.p}" --use-conda --use-singularity --nolock --rerun-incomplete -F
