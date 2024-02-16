#!/bin/bash
#SBATCH --output=/home/genouest/cnrs_umr6553/erolland/Ldhat_%j.out
#SBATCH --error=/home/genouest/cnrs_umr6553/erolland/Ldhat_%j.err
#SBATCH --mail-type=END
#SBATCH --mail-user=elise.rolland@univ-rennes1.fr
#SBATCH --chdir=/scratch/erolland
#SBATCH --job-name=LDhat
#SBATCH --mem=80G
#SBATCH --tasks-per-node 10

#==================================
# LDHAT PIPELINE
#==================================
# LDhat is a reference software to estimate fine-scale recombination rates based on polymorphism data (i.e. a reference genome and a vcf)
# Associated to LDhot for hotspot detection
# A complete framework

# This script is made for running a single analysis of LDhat to a single population
# with as input files a reference genome and a vcf file
# For a single chromosome

# The vcf file must contain one population, one chromosome

# LOAD ENVIRONMENT


# CREATE A DIRECTORY WHERE TO RUN ALL ANALYSES
# A directory named after the study name...
wdgroup=/groups/landrec # 
wd=/scratch/erolland # Working directory
prefix=Hordeum_vulgare_Dreissig2019 # The ID of the dataset
indiv=20 Size of the population sampled
pop=20_demo # Name of the population sampled

# Chromosome names MUST BE integers (e.g. 1) or integers+characters (e.g. 1A) in the vcf
# Chromosome names given in argument to the script
chrom=1 # The chromosome to sample

input_data=$wdgroup/data/polymorphism_data/$prefix # Input from data directory
input=$wd/ldhat/$prefix/input/population_$pop # Input directory
output=$wd/ldhat/$prefix/output/population_$pop # Output directory

mkdir -p $input
mkdir -p $output

# PARALLEL
ncores=10 # Number of cores required on the cluster, passed in argument to the script

# ENVIRONMENT CONDA
. /local/env/envconda.sh
conda activate ~/env_smc
cd $wd

#-----------------------------------
# TRIMMING THE VCF
#-----------------------------------
# All VCFs must be trimmed to the same quality filters to be comparable among different species and studies

# The trimmed VCF is then copied to the working directory

# Trimming missing data
missing=0.1 # The proportion of missing data allowed

# Trimming Minor Allele Frequency
MAF=0.05 # MAF threshold; 0.05 is a standard value, commonly accepted

. /softs/local/env/envvcftools-0.1.16.sh

# Keep only a subset of individuals in the trimmed dataset for analyses
vcftools --gzvcf $input_data/$prefix.SNP.vcf.gz --out $input/$prefix.filter.SNP --recode --maf $MAF --max-missing $missing --keep ~/data/$prefix\_74.ind
mv $input/$prefix.filter.SNP.recode.vcf $input/$prefix.filter.SNP.vcf
gzip $input/$prefix.filter.SNP.vcf
vcftools --gzvcf $input/$prefix.filter.SNP.vcf.gz --missing-indv --out $input/$prefix.ind_missing
awk '$5 > 0.7' $input/$prefix.ind_missing.imiss | cut -f1 > $input/lowDP.ind
vcftools --gzvcf $input/$prefix.filter.SNP.vcf.gz --remove $input/lowDP.ind --recode --out $input/$prefix.filter.SNP # Filter individual with more than 0.5 missing data
mv $input/$prefix.filter.SNP.recode.vcf $input/$prefix.filter.SNP.vcf
gzip -f $input/$prefix.filter.SNP.vcf


#----------------------------------
# DETECTING POPULATION STRUCTURE
#----------------------------------
# Population structure in the sampled population can be a source of biases in further inferences of recombination rates
# Hence, once the population is sampled from the global dataset, we need to test if the smaple present some genetic structure




#-----------------------------------
# PHASING VCF
#-----------------------------------
# "The method accommodates both phased or haplotype and unphased or genotype data, with arbitrary levels of missing data."
# Yet, estimates are better with phased hapltotypes.
# Unphased haplotypes can overestimate recombination rate, hence false postive hotspots may be detected (Booker et al., 2017 ; Shanfelter et al., 2019)

echo =========================================================================
echo "    Phase the VCF with SHAPEIT2                                        "
echo =========================================================================

# Split dataset by chromosome and remove SNP with missing data for all individuals
vcftools --gzvcf $input/$prefix.filter.SNP.vcf.gz --out $input/$prefix.chromosome$chrom.SNP --recode --chr $chrom --max-missing-count $indiv
mv $input/$prefix.chromosome$chrom.SNP.recode.vcf $input/$prefix.chromosome$chrom.SNP.vcf

# Compute the nucleotide diversity in a window of 100000 bp
vcftools --vcf $input/$prefix.chromosome$chrom.SNP.vcf --window-pi 100000 --out $input/$prefix.chromosome$chrom

# Compute the effective size of the chromosome
eff_size=$(awk '{sum+=$5} END {sum=(sum/NR)/(4*(10** -8)); print sum}' $input/$prefix.chromosome$chrom.windowed.pi)
eff_size=${eff_size/.*}
echo -e "\n $eff_size \n"

# Phasing with Shapeit
$wd/ldhat/shapeit/bin/shapeit --input-vcf $input/$prefix.chromosome$chrom.SNP.vcf  --output-max $input/$prefix.chromosome$chrom.phased --effective-size $eff_size --window 1 --thread $ncores --output-log $input/$prefix.chromosome$chrom.phased --force

# Convert output file from Shapeit to vcf
$wd/ldhat/shapeit/bin/shapeit -convert --input-haps $input/$prefix.chromosome$chrom.phased --output-vcf $input/$prefix.chromosome$chrom.phased.vcf --output-log $input/$prefix.chromosome$chrom.convert
gzip $input/$prefix.chromosome$chrom.phased.vcf

# Add the length of the contig to the vcf, need to compute the effective size
##samtools faidx $wd/data/genome/hordeum_vulgare/*/fasta/*.fna.gz
##bcftools reheader -f $wd/data/genome/hordeum_vulgare/*/fasta/*.fna.gz.fai $input/$prefix.chromosome$chrom.phased.vcf -o "$input/$prefix.chromosome$chrom.temp.vcf"
##mv $input/$prefix.chromosome$chrom.temp.vcf $input/$prefix.chromosome$chrom.phased.vcf
##sed -i '/CAJE/d' $input/$prefix.chromosome$chrom.phased.vcf


#-----------------------------------
# RANDOM SAMPLING OF THE POPULATION
#-----------------------------------
# Random sampling of N individuals/sequences in the dataset
# to make comparisons between samples of equal size
# Besides, consider selfers in sampling -> haploid sampling but twice the number of individuals

# Filtering individuals
vcftools --gzvcf $input/$prefix.chromosome$chrom.phased.vcf.gz --out $input/$prefix\_$indiv.SNP --keep $input/$prefix\_$indiv.ind --recode
#vcftools --gzvcf $input/$prefix.chromosome$chrom.phased.vcf.gz --out $input/$prefix\_$indiv.SNP --keep ~/data/cluster1_full.ind --recode

mv $input/$prefix\_$indiv.SNP.recode.vcf $input/$prefix.pop_$pop.chromosome$chrom.SNP.vcf
gzip $input/$prefix.pop_$pop.chromosome$chrom.SNP.vcf


#-----------------------------------
# PRE-ESTIMATE DEMOGRAPHY
#-----------------------------------
# Past demography is an important parameter that can biases inferences of LD-based recombination rates
# Hence we estimated an approximate past demography to pass to LDpop
# LDpop further compute a demography-aware look-up table for LDhat

echo =========================================================================
echo "    Estimate past demography with SMC++                                "
echo =========================================================================

#. /local/env/envsingularity-3.6.4.sh

# Convert your VCF(s) to the SMC++ input format with vcf2smc
# samples=$(cat $input/$prefix\_$indiv.ind | tr -s '\n' ',')
# singularity run smcpp_latest.sif vcf2smc $input/$prefix.chromosome$chrom.phased.vcf.gz $input/$prefix/$prefix.chromosome$chrom.smc.gz $chrom pop1:$samples --cores $ncores
# Index the vcf
#sample=$(sed -n -E '/ERX/p' $input/$prefix.chromosome$chrom.hap.vcf)
#sample=$(echo $sample | awk '{split($0, temp, " "); for(i=10; i < length(temp)+1; i++) {print temp[i]}}')
#sample=$(echo $sample | tr -s '[:space:]' ',')
# sample=$(cat $input/../$prefix\_$indiv.ind | tr -s '\n' ',')
#sample=${sample%?}

##bgzip "$input/$prefix.chromosome$chrom.hap.vcf"
##tabix -p vcf $input/$prefix.chromosome$chrom.hap.vcf.gz --csi

#vcftools --gzvcf $input/$prefix.chromosome$chrom.phased.vcf.gz --freq --out $input/$prefix.freq
#bcftools roh 


#bgzip "$input/$prefix.chromosome$chrom.mask.bed"
#tabix -p bed "$input/$prefix.chromosome$chrom.mask.bed.gz"
##singularity run --bind $wd:$wd smcpp_latest.sif vcf2smc "$input/$prefix.chromosome$chrom.hap.vcf.gz" "$input/$prefix.chromosome$chrom.smc.gz" $chrom pop1:$sample --cores $ncores --mask "$input/Hordeum_vulgare_Dreissig2019.chromosome1.mask.bed.gz" --length 502612092

mu=1.25e-8 # Per-generation mutation rate
#singularity run smcpp_latest.sif estimate -o  $input/smc/ $mu $input/$prefix.chromosome$chrom.smc.gz --cores $ncores 
##singularity run --bind $input:$input smcpp_latest.sif estimate -o  $input/smc/ $mu $input/$prefix.chromosome$chrom.smc.gz --cores $ncores 

##mv iterate.dat $input/smc/
# Visualize the results using plot
##singularity run --bind $input:$input smcpp_latest.sif plot $input/plot.pdf $input/smc/model.final.json -c
# Export posteriors
#singularity run smcpp_latest.sif posterior $input/smc/model.final.json $input/smc/posterior.smc \
#             $input/$prefix.chromosome$chrom.smc.gz --heatmap $input/smc/heatmap.png --colorbar

#-----------------------------------
# PRE-ESTIMATE THETA
#-----------------------------------


# The value of Theta
theta=$(awk -v mu="$mu" -v eff_size="$eff_size" 'BEGIN {print 4*mu*eff_size}')
echo -e "\n Theta $theta \n"

#-----------------------------------
# LDPOP
#-----------------------------------
# Use LDpop to create a look-up table taking into account the underlying demographic history
# https://github.com/popgenmethods/ldpop
# Use run/ldtable.py to create a lookup table.
# By default run/ldtable.py uses an exact algorithm to compute the likelihoods. To use a reasonable approximation that is much faster and scales to larger sample sizes, use the flag --approx.
. /local/env/envpython-3.7.6.sh

# Install LDpop
#git clone https://github.com/popgenmethods/ldpop.git
#cd ldpop
#pip install .
# Run run/ldtable.py --help for the manual
#N=20 # Number of individuals
N=$indiv
num_rh=10
max_rh=10

##N0=$( jq '.model.N0' $input/smc/model.final.json)
##coal_sizes=$( jq '.model.y' $input/smc/model.final.json | tr -d '[:space:][]') 
##coal_sizes=$(echo $coal_sizes | awk '{split($0, temp, ","); for(i=1; i < length(temp)+1; i++) {a=exp(temp[i]); print a}}')
##coal_sizes=$(echo $coal_sizes | tr -s '[:space:]' ',')
##coal_sizes=${coal_sizes%?}
##coal_sizes="$N0","$coal_sizes"
##coal_times=$( jq '.model.knots' $input/smc/model.final.json | tr -d '[:space:][]') #enlever premier element
##$wd/ldhat/ldpop/run/ldtable.py -n $N -th $theta -rh $num_rh,$max_rh -s $coal_sizes -t $coal_times --approx --cores $ncores --log "$output/ldpop.chromosome$chrom.cluster$cluster.log" > $output/lk_ldpop.chromosome$chrom.txt
$wd/ldhat/ldpop/run/ldtable.py -n $N -th $theta -rh $num_rh,$max_rh --cores $ncores  --log "$output/$prefix.pop_$pop.chromosome$chrom.ldpop.log" > "$output/$prefix.pop_$pop.chromosome$chrom.lk_ldpop.txt"

conda deactivate

#-----------------------------------
# LDHAT, estimating Rho
#-----------------------------------
# LDhat is a suite of programs used sequentially to estimate the recombination rate along the chromosome.
# LDhat is compiled locally
PathToLDhat=$wd/ldhat/LDhat-cluster

# INPUT FILES
# vcftools can convert VCFs to input files for LDhat with option
#         --ldhat
#         --ldhelmet
#         --ldhat-geno
#           These options output data in LDhat/LDhelmet format. This option requires the "--chr" filter option to also be used. The two first  options  output phased  data only, and therefore also implies "--phased" be used, leading to unphased individuals and genotypes being excluded. For LDhelmet, only snps will be considered, and therefore it implies "--remove-indels". The second option treats all of the data as unphased, and  therefore  outputs LDhat  files  in  genotype/unphased format. Two output files are generated with the suffixes ".ldhat.sites" and ".ldhat.locs", which correspond to the LDhat "sites" and "locs" input files respectively; for LDhelmet, the two files generated  have  the  suffixes  ".ldhelmet.snps"  and  ".ldhelmet.pos", which corresponds to the "SNPs" and "positions" files.
# vcftools --gzvcf $input/$prefix\_$indiv.SNP.vcf.gz --ldhat-geno --chr $chrom --out $output/$prefix.chromosome$chrom

vcftools --gzvcf $input/$prefix.pop_$pop.chromosome$chrom.SNP.vcf.gz --ldhat --chr $chrom --out $output/$prefix.pop_$pop.chromosome$chrom

# For haploid correct the indication on the first line of the .sites
sed -e '/\-1/{N;d}' -i $output/$prefix.pop_$pop.chromosome$chrom.ldhat.sites
indiv2=$(($indiv*2))
sed "0,/$indiv2/{s//$indiv/}" -i $output/$prefix.pop_$pop.chromosome$chrom.ldhat.sites

# sed '0,/.2/{s//\t1/}' -i $output/$prefix.chromosome$chrom.ldhat.sites 

# Analyses are performed in the current directory
cd $output

# CONVERT
# This program produce input files of the proper format for LDHat
# Generate sites.txt and locs.txt
#.$PathToLDhat/convert -seq <file_name> # Useless since vcftools can produce input files directly from a phased vcf

# PAIRWISE
# Pairwise can be used to produce a suitable estimate of theta first, that will be used in COMPLETE
#theta=0.01

# COMPLETE - Look-up table
# Generates a lookup table required for all analyses (except
# under certain circumstances, see pairwise). Input number of sequences,
# theta per site (assumes a two-allele symmetric mutation model), and the
# grid size for 4Ner (in terms of maximum value and number of points in
# the grid). The program calculates the coalescent likelihood of all possible
# two-locus haplotype configurations using the importance sampling method
# of Fearnhead and Donnelly [4].
# $PathToLDhat/complete -n $indiv -rhomax 100 -n_pts 101 \
# 			-theta $theta \
# 			-split $ncore \
# 			-element 0

# . /softs/local/env/envparallel-20190122.sh

# cat ~/group.txt |parallel -j $ncore $wd/ldhat/src/complete.sh

# cat new_lk* > lk_complete.txt
# rm new_lk*
# Look-up tables can be very long to generate (many days), hence their computation can be speed up
# By paralellising computation (split) and concatenating independant elements in a a single look-up table


# LKGEN
# Can be used to generate a custom likelihood look-up table from one generated in the LDhat website
# Get the number of samples in the vcf
# On completion, the output file, new_lk.txt should be renamed before future
# analyses are carried out. Note that lkgen can only be used to generate lookup
# tables for numbers of chromosomes less than that of the input lookup file.
# For diploids, count 2 chromosomes per individual
# e.g. nseq=100 means 5 diploid individuals
#. /softs/local/env/envbcftools-1.9.sh

# Get the number of samples using bcftools
#nseq=$(bcftools stats $input/$prefix\_$indiv.SNP.vcf.gz | grep samples | awk '{ print $6 }')
#nseq=10
# Multiply by 2 for diploids
# nseq=$(echo $((2*$nseq)))

#gunzip $PathToLDhat/lk_files/lk_n192_t0.001.gz
#lkfile=$PathToLDhat/lk_files/lk_n192_t0.001
#  lkgen is less computationnally intensive than complete
# ./lkgen -lk <file> -nseq <int>
# The new lk file is added in the current directory
#$PathToLDhat/lkgen -lk $lkfile -nseq $nseq

# INTERVAL
# Estimates a variable recombination rate using a Bayesian reversible-
# jump MCMC scheme [7] under the crossing-over model (only). As for
# pairwise, the sites, locs, and lookup table files are required.
# ./interval -seq <file_name> -loc <file_name> -lk <file_name>

#Interval options 
# iter=2000
# samp=100
 
# echo =========================================================================
# echo  "	     	  Running $prefix with $iter iterations"
# echo  "	     	  sampling every $samp iterations"
# echo =========================================================================

# Compute directly in the output directory where input files must be
#cd $output

# $PathToLDhat/interval -seq $prefix.chromosome$chrom.ldhat.sites \
#	 -loc $prefix.chromosome$chrom.ldhat.locs \
#	 -lk lk_complete.chromosome$chrom.txt  \
#	 -its "$iter" \
#	 -bpen 5 \
#	 -samp "$samp"

# $PathToLDhat/interval -seq $prefix.chromosome$chrom.ldhat.sites -loc $prefix.chromosome$chrom.ldhat.locs -lk lk_complete.txt -its "$iter" -bpen 5 -samp "$samp"

# cd $wd

# RHOMAP
# An alternative to INTERVAL
# As with interval, rhomap estimates a variable recombination
# rate using a Bayesian reversible-jump MCMC scheme. However, rhomap
# specifically fits a model with recombination hotspots on a background of
# low rate variation. In addition, the method rescales (flattens) the
# composite likelihood to improve estimation of posterior probabilities.

# Rhomap option
iter=1000000
samp=5000
burn=200000

echo =========================================================================
echo  "                 Running $prefix with $iter iterations"
echo  "                 sampling every $samp iterations"
echo =========================================================================


$PathToLDhat/rhomap -seq $prefix.pop_$pop.chromosome$chrom.ldhat.sites \
	-loc $prefix.pop_$pop.chromosome$chrom.ldhat.locs \
	-lk $prefix.pop_$pop.chromosome$chrom.lk_ldpop.txt \
	-its "$iter" \
	-samp "$samp" \
	-burn "$burn" \
	-prefix "$prefix.pop_$pop.chromosome$chrom."

# $PathToLDhat/rhomap -seq $prefix.chromosome$chrom.hap.ldhat.sites -loc $prefix.chromosome$chrom.hap.ldhat.locs -lk lk_ldpop.txt -its "$iter" -samp "$samp" -burn "$burn"

# STAT
# Summarises the output from interval in terms of the average, median, 2.5th
# percentile and 97.5th percentile of the estimated recombination rate between
# each pair of SNPs.


#-----------------------------------
# SUMMARIZING AND DISPLAYING RESULTS
#-----------------------------------

echo =========================================================================
echo "  Summarising final Ne r estimates --\n\t burning 5 first iterations "
echo =========================================================================

$PathToLDhat/stat -input $prefix.pop_$pop.chromosome$chrom.rates.txt \
	-burn 5 -loc $prefix.pop_$pop.chromosome$chrom.ldhat.locs \
	-prefix "$prefix.pop_$pop.chromosome$chrom."


#-----------------------------------
# LDHOT, hotspot detection
#-----------------------------------
# LDhot uses a LDhat recombination landscapes to infer hotspots
# i.e. regions of recombination higher thant the average background
PathToLDhot=$wd/ldhat/LDhot

$PathToLDhot/ldhot --seq $prefix.pop_$pop.chromosome$chrom.ldhat.sites \
	--loc $prefix.pop_$pop.chromosome$chrom.ldhat.locs \
	--lk $prefix.pop_$pop.chromosome$chrom.lk_ldpop.txt \
	--res "$prefix.pop_$pop.chromosome$chrom.res.txt" \
	--nsim 100 \
	--out $prefix.pop_$pop.chromosome$chrom

$PathToLDhot/ldhot_summary --res $prefix.pop_$pop.chromosome$chrom.res.txt \
	--hot "$prefix.pop_$pop.chromosome$chrom.hotspots.txt" \
	--out "$prefix.pop_$pop.chromosome$chrom"

