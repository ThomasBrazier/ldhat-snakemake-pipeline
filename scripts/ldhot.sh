#!/bin/bash
#SBATCH --job-name=Ldhot
#SBATCH --mail-type=END
#SBATCH --mail-user=elise.rolland@univ-rennes1.fr
#SBATCH --chdir=/omaha-beach/erolland
#SBATCH --output=/home/genouest/cnrs_umr6553/erolland/Ldhot_%j.out
#SBATCH --error=/home/genouest/cnrs_umr6553/erolland/Ldhot_%j.err

wd=/omaha-beach/erolland
wd_group=/groups/landrec # Working directory
prefix=Hordeum_vulgare_Dreissig2019 # The ID of the dataset
indiv=18 # Size of the population sampled

# Chromosome names MUST BE integers (e.g. 1) or integers+characters (e.g. 1A) in the vcf
# Chromosome names given in argument to the script
chrom=1
cluster=1

input_data=$wdgroup/data/polymorphism_data/$prefix # Input from data directory
#input=$wd/ldhat/$prefix/input/indiv_$indiv/chrom$chrom # Input directory
input=$wd/ldhat/$prefix/input/indiv_$indiv/cluster$cluster # Input directory
#output=$wd/ldhat/$prefix/output/indiv_$indiv/chrom$chrom # Output directory
output=$wd/ldhat/$prefix/output/indiv_$indiv/cluster$cluster # Output directory

cd $output

#-----------------------------------
# LDHOT, hotspot detection
#-----------------------------------


PathToLDhot=$wd_group/ldhat/LDhot

$PathToLDhot/ldhot --seq $prefix.chromosome$chrom.cluster$cluster.ldhat.sites \
        --loc $prefix.chromosome$chrom.cluster$cluster.ldhat.locs \
        --lk lk_ldpop.chromosome$chrom.cluster$cluster.txt \
        --res res.txt \
        --nsim 100 \
        --out $prefix.chromosome$chrom.cluster$cluster.sim

$PathToLDhot/ldhot_summary --res res.sim.txt \
        --hot "$prefix.chromosome$chrom.cluster$cluster.sim.hotspots.txt" \
        --out "$prefix.chromosome$chrom.cluster$cluster.sim"

