#!/bin/bash
#SBATCH --output=/home/genouest/cnrs_umr6553/erolland/Ldhat_%j.out
#SBATCH --error=/home/genouest/cnrs_umr6553/erolland/Ldhat_%j.err
#SBATCH --mail-type=END
#SBATCH --mail-user=elise.rolland@univ-rennes1.fr
#SBATCH --chdir=/scratch/erolland
#SBATCH --job-name=LDhat
#SBATCH --mem=30G
#SBATCH --tasks-per-node 10

wdgroup=/groups/landrec # Working directory
wd=/scratch/erolland
prefix=Hordeum_vulgare_Dreissig2019 # The ID of the dataset
indiv=30 # Size of the population sampled
pop=1

# Chromosome names MUST BE integers (e.g. 1) or integers+characters (e.g. 1A) in the vcf
# Chromosome names given in argument to the script
chrom=2 #The chromosome to sample

input_data=$wdgroup/data/polymorphism_data/$prefix # Input from data directory
#input=$wd/ldhat/$prefix/input/indiv_$indiv/population_$pop/chrom$chrom # Input directory
input=$wd/ldhat/$prefix/input/population_$pop/mask/chrom$chrom # Input directory
output=$wd/ldhat/$prefix/output/population_$pop/mask/chrom$chrom # Output directory

mkdir -p $input
mkdir -p $output

# PARALLEL
ncores=10 # Number of cores required on the cluster, passed in argument to the script

# ENVIRONMENT CONDA
. /local/env/envconda.sh
conda activate ~/env_smc
cd $wd

missing=0.1 # The proportion of missing data allowed

# Trimming Minor Allele Frequency
MAF=0.05 # MAF threshold; 0.05 is a standard value, commonly accepted

# Trimming maximum missing data per individual
maxmissing=0.9

. /softs/local/env/envvcftools-0.1.16.sh

# Keep only a subset of individuals in the trimmed dataset for analyses
vcftools --gzvcf $input_data/$prefix.SNP.vcf.gz --out $input/../$prefix.filter.SNP --recode --maf $MAF --max-missing $missing --keep ~/data/$prefix\_74.ind
mv $input/../$prefix.filter.SNP.recode.vcf $input/../$prefix.filter.SNP.vcf
gzip $input/../$prefix.filter.SNP.vcf
vcftools --gzvcf $input/../$prefix.filter.SNP.vcf.gz --missing-indv --out $input/../$prefix.ind_missing
awk '$5 > $maxmissing' $input/../$prefix.ind_missing.imiss | cut -f1 > $input/../lowDP.ind
vcftools --gzvcf $input/../$prefix.filter.SNP.vcf.gz --remove $input/../lowDP.ind --recode --out $input/../$prefix.filter.SNP # Filter individual with more than 0.5 missing data

vcftools --gzvcf $input/../$prefix.filter.SNP.vcf.gz --out $input/$prefix.chromosome$chrom.SNP --recode --chr $chrom --max-missing-count $indiv
mv $input/$prefix.chromosome$chrom.SNP.recode.vcf $input/$prefix.chromosome$chrom.SNP.vcf

# Compute the nucleotide diversity in a window of 100000 bp
vcftools --vcf $input/$prefix.chromosome$chrom.SNP.vcf --window-pi 100000 --out $input/$prefix.chromosome$chrom

# Compute the effective size of the chromosome
eff_size=$(awk '{sum+=$5} END {sum=(sum/NR)/(4*(10** -8)); print sum}' $input/$prefix.chromosome$chrom.windowed.pi)
eff_size=${eff_size/.*}
echo -e "\n $eff_size \n"

$wd/ldhat/shapeit/bin/shapeit --input-vcf $input/$prefix.chromosome$chrom.SNP.vcf  --output-max $input/$prefix.chromosome$chrom.phased --effective-size $eff_size --window 1 --thread $ncores --output-log $input/$prefix.chromosome$chrom.phased --force

# Convert output file from Shapeit to vcf
$wd/ldhat/shapeit/bin/shapeit -convert --input-haps $input/$prefix.chromosome$chrom.phased --output-vcf $input/$prefix.chromosome$chrom.phased.vcf --output-log $input/$prefix.chromosome$chrom.convert
gzip $input/$prefix.chromosome$chrom.phased.vcf

vcftools --gzvcf $input/$prefix.chromosome$chrom.phased.vcf.gz --out $input/$prefix.pop_$pop.chromosome$chrom.SNP --keep ~/data/cluster1_full.ind --recode
mv $input/$prefix.pop_$pop.chromosome$chrom.SNP.recode.vcf $input/$prefix.pop_$pop.chromosome$chrom.SNP.vcf
gzip $input/$prefix.pop_$pop.chromosome$chrom.SNP.vcf

. /local/env/envpython-3.7.6.sh


theta=0.01
N=$indiv
num_rh=10
max_rh=10
$wd/ldhat/ldpop/run/ldtable.py -n $N -th $theta -rh $num_rh,$max_rh --cores $ncores --log "$output/ldpop.chromosome$chrom.log" > $output/$prefix.pop_$pop.chromosome$chrom.lk_ldpop.txt
conda deactivate


cd $output
PathToLDhat=$wd/ldhat/LDhat-cluster

vcftools --gzvcf $input/$prefix.pop_$pop.chromosome$chrom.SNP.vcf.gz --ldhat --chr $chrom --out $output/$prefix.pop_$pop.chromosome$chrom

# For haploid correct the indication on the first line of the .sites
sed -e '/\-1/{N;d}' -i $output/$prefix.pop_$pop.chromosome$chrom.ldhat.sites
indiv2=$(($indiv*2))
sed "0,/$indiv2/{s//$indiv/}" -i $output/$prefix.pop_$pop.chromosome$chrom.ldhat.sites



iter=6000000
samp=5000

echo =========================================================================
echo  "           Running $prefix with $iter iterations"
echo  "           sampling every $samp iterations"
echo =========================================================================

# Compute directly in the output directory where input files must be
#cd $output

$PathToLDhat/interval -seq $prefix.pop_$pop.chromosome$chrom.ldhat.sites \
        -loc $prefix.pop_$pop.chromosome$chrom.ldhat.locs \
        -lk $prefix.pop_$pop.chromosome$chrom.lk_ldpop.txt  \
        -its "$iter" \
        -bpen 5 \
        -samp "$samp" \
        -prefix "$prefix.pop_$pop.chromosome$chrom."

$PathToLDhat/stat -input $prefix.pop_$pop.chromosome$chrom.rates.txt \
        -burn 200 -loc $prefix.pop_$pop.chromosome$chrom.ldhat.locs \
        -prefix "$prefix.pop_$pop.chromosome$chrom."


