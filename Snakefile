"""
Pipeline to estimate fine-scale recombination maps from polymorphism data
"""

"""
Configuration of the analysis
i.e. dataset, name of chromosome, population to sample
"""

"""In addition to the configfile statement, config values can be overwritten via the command line"""
dataset=config["dataset"] # Name of your dataset directory and prefix of your vcf file
chrom=config["chrom"] # Name of the chromosome to analyse in your 'sample' dataset
K=str(config["K"])
pop=str(config["pop"])
bpen=config["interval.bpen"]


wdir=config["workingdir"] + dataset
wdirpop=config["workingdir"] + dataset + "/K" + K + ".pop" + pop


wildcard_constraints:
    wdir=wdir,
    wdirpop=wdirpop,
    dataset=dataset,
    chrom=chrom,
    K=K,
    pop=pop,
    bpen=bpen


rule all:
    """
    One ring to rule them all"
    """
    input:
        target = expand("{wdirpop}/ldhot/{dataset}.{chrom}.bpen{bpen}.hotspots.txt",wdirpop=wdirpop,dataset=dataset,chrom=chrom,bpen=bpen),
    shell:
        "echo 'Finished'"


rule poplist:
    """
    Sample individuals from a given genetic cluster (number of K to retain and name of the cluster) specified in config.yaml
    """
    output:
        "{wdirpop}/poplist"
    log:
        "{wdirpop}/logs/poplist.log"
    conda:
        "envs/vcftools.yaml"
    shell:
        """
        perl -ane '$r = 0; for my $i (1 .. $#F) {{$r = $i if $F[$i] > $F[$r];}} print $r + 1, " ";' < {wdir}/structure/faststructure.{K}.meanQ > {wdir}/structure/cluster.tmp
        sed -i -e "s/\\s\\+/\\n/g" {wdir}/structure/cluster.tmp
        paste {wdir}/indlist {wdir}/structure/cluster.tmp > {wdir}/structure/cluster.tmp2
        paste {wdir}/structure/cluster.tmp2 {wdir}/structure/faststructure.{K}.meanQ > {wdir}/structure/cluster
        awk -F' ' '{{if($2=={pop}) print $1}}' {wdir}/structure/cluster > {wdirpop}/poplist
        """


rule sampling_pop:
    """
    A first step to trim every 'sample' dataset to the same quality criteria
    """
    input:
        "{wdirpop}/poplist"
    output:
    	"{wdirpop}/{dataset}.pop.vcf.gz"
    log:
        "{wdirpop}/logs/{dataset}.sampling_pop.log"
    conda:
        "envs/vcftools.yaml"
    shell:
        """
        vcftools --gzvcf {wdir}/trimmed.vcf.gz --out {wdirpop}/out --recode --keep {wdirpop}/poplist --maf config[maf] --max-missing config[maxmissing] --min-alleles 2 --max-alleles 2
        mv {wdirpop}/out.recode.vcf {wdirpop}/{dataset}.pop.vcf
        bgzip -f {wdirpop}/{dataset}.pop.vcf
        """


rule effective_size:
    """
    Estimate population effective size
    """
    input:
        "{wdirpop}/{dataset}.pop.vcf.gz"
    output:
        "{wdirpop}/statistics/{dataset}.effective_size"
    log:
        "{wdirpop}/logs/{dataset}.effective_size.log"
    conda:
        "envs/vcftools.yaml"
    params:
        windowsize=100000
    shell:
        """
        vcftools --gzvcf {input} --window-pi {params.windowsize} --out {wdirpop}/statistics/vcftools
        eff_size=$(awk '{{sum+=$5}} END {{sum=(sum/NR)/(4*{config[mu]}); print sum}}' {wdirpop}/statistics/vcftools.windowed.pi)
        echo ${{eff_size/.*}} > {output}
        """


# TODO Take care of run of homozygosity (mask): raise WARNING
# Plink  search for ROH
# TODO How to treat selfing: https://github.com/popgenmethods/smcpp/issues/82
# TODO optimization: thinning, spline, knots, cores
rule demography:
    """
    Pre-estimate demography with smc++ to use it as a prior in likelihood
    Subset a small sample size is enough to achieve good accuracy (memory issues)
    Each chromosome is treated independently in vcf2smc
    then a consensus model is estimated on all chromosomes
    i.e. same demographic history for all chromosomes

    Terhorst J, Kamm JA, Song YS. Robust and scalable inference of
    population history from hundreds of unphased whole genomes.
    Nat Genet. 2017 Feb;49(2):303–9.
    """
    input:
        "{wdirpop}/statistics/{dataset}.effective_size"
    output:
        "{wdirpop}/smc/{dataset}.model.final.json"
    log:
        "{wdirpop}/logs/{dataset}.demography.log"
    conda:
        "envs/vcftools.yaml"
    shell:
        """
        tabix -f -p vcf {wdirpop}/{dataset}.pop.vcf.gz --csi
        shuf -n {config[subset]} --random-source=<(yes {config[seed]}) {wdirpop}/poplist > {wdirpop}/subsetpop
        individuals=$(cat {wdirpop}/subsetpop | tr '\n' ',' | sed "s/,$//g")
        chromosomes=$(zcat {wdirpop}/{dataset}.pop.vcf.gz | awk '{{ print $1 }}' | sort | uniq | grep -v '^#')
        for c in $chromosomes; do
        singularity exec --bind $PWD:/mnt smcpp.sif smc++ vcf2smc --ignore-missing /mnt/{wdirpop}/{dataset}.pop.vcf.gz /mnt/{wdirpop}/smc/vcf2smc.$c $c pop1:$individuals
        done
        paths=$(for i in $chromosomes; do echo /mnt/{wdirpop}/smc/vcf2smc.$i; done | tr '\n' ' ')
        singularity exec --bind $PWD:/mnt smcpp.sif smc++ estimate -o /mnt/{wdirpop}/smc/ {config[mu]} $paths --cores {config[cores]}
        mv {wdirpop}/smc/model.final.json {wdirpop}/smc/{dataset}.model.final.json
        singularity exec --bind $PWD:/mnt smcpp.sif smc++ plot /mnt/{wdirpop}/smc/plot_chromosome.pdf /mnt/{wdirpop}/smc/{dataset}.model.final.json -c
        singularity exec --bind $PWD:/mnt smcpp.sif smc++ posterior /mnt/{wdirpop}/smc/{dataset}.model.final.json /mnt/{wdirpop}/smc/{dataset}.posterior.smc $paths
        mv ~/iterate.dat {wdirpop}/smc/{dataset}.{chrom}.iterate.dat
        """


# From now, analyses are performed on individual chromosomes
rule split_chromosome:
    """
    Split the entire vcf into one vcf per chromosome
    Deal with missing data: remove SNPs with missing data for more than 90% of individuals
    """
    input:
        "{wdirpop}/smc/{dataset}.model.final.json"
    output:
        "{wdirpop}/{dataset}.chromosome.{chrom}.vcf.gz"
    log:
        "{wdirpop}/logs/{dataset}.split_chromosome.{chrom}.log"
    conda:
        "envs/vcftools.yaml"
    shell:
        """
        vcftools --gzvcf {wdirpop}/{dataset}.pop.vcf.gz --out {wdirpop}/out --recode --chr {chrom}
	vcftools --vcf {wdirpop}/out.recode.vcf --out out2 --recode --maf config[maf] --max-missing config[maxmissing]
        rm {wdirpop}/out.recode.vcf
	mv {wdirpop}/out2.recode.vcf {wdirpop}/{dataset}.chromosome.{chrom}.vcf
        bgzip -f {wdirpop}/{dataset}.chromosome.{chrom}.vcf
        """

rule phasing_vcf:
    """
    Phase the vcf with Shapeit2
    Unphased haplotypes can overestimate recombination rate, hence false postive hotspots may be detected
    https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html#output
    """
    input:
        "{wdirpop}/{dataset}.chromosome.{chrom}.vcf.gz"
    output:
        "{wdirpop}/{dataset}.chromosome.{chrom}.phased.vcf.gz"
    log:
        "{wdirpop}/logs/{dataset}.chromosome.{chrom}.phasing_vcf.log"
    conda:
        "envs/shapeit.yaml"
    shell:
        """
        shapeit --input-vcf {wdirpop}/{dataset}.chromosome.{chrom}.vcf.gz --output-max {wdirpop}/{dataset}.phased.chromosome.{chrom} --effective-size $(cat {wdirpop}/statistics/{dataset}.effective_size) --window 1 --thread {config[cores]} --output-log {wdirpop}/logs/{dataset}.chromosome.{chrom}.shapeit.log --force
        shapeit -convert --input-haps {wdirpop}/{dataset}.phased.chromosome.{chrom} --output-vcf {wdirpop}/{dataset}.chromosome.{chrom}.phased.vcf --output-log {wdirpop}/logs/{dataset}.chromosome.{chrom}.shapeit.convert.log
        # replace header in vcf to keep information of contig length
        zcat {wdir}/{dataset}.vcf.gz | grep '^#' > {wdirpop}/newheader
        cat {wdirpop}/newheader | grep -v '^#CHROM' > {wdirpop}/newheader2
        cat {wdirpop}/{dataset}.chromosome.{chrom}.phased.vcf | grep '^#CHROM' > {wdirpop}/colnames
        cat {wdirpop}/{dataset}.chromosome.{chrom}.phased.vcf | grep -v '^#' > {wdirpop}/newvcf
        cat {wdirpop}/newheader2 {wdirpop}/colnames {wdirpop}/newvcf > {wdirpop}/{dataset}.chromosome.{chrom}.phased.vcf
        bgzip -f {wdirpop}/{dataset}.chromosome.{chrom}.phased.vcf
        rm {wdirpop}/newheader {wdirpop}/newheader2 {wdirpop}/colnames {wdirpop}/newvcf
        """



# # TODO Problem: vcf must be diploid in pyrho
# # "ploidy set to 1 if using phased data and 2 for unphased genotype data"
# # https://knausb.github.io/vcfR_documentation/dip_to_hap.html
# rule pseudodiploid:
#     """
#     Pseudo-diploidization
#     to account for high selfing rates and homozygosity
#     Make pseudo-diploids for compatibility with pyrho
#     """
#     input:
#         "{wdirpop}/{dataset}.phased.vcf.gz"
#     output:
#         "{wdirpop}/{dataset}.haploid.vcf.gz"
#     log:
#         "{wdirpop}/logs/{dataset}.haploid.log"
#     conda:
#         "envs/vcftools.yaml"
#     shell:
#         """
#         """



# # TODO A step to reduce sample size to a random subset
# # Memory issues when sample size is high
# rule subset_pyrho:
#     """
#     Subset a random sample of individuals
#     """
#     input:
#         "{wdirpop}/{dataset}.chromosome.{chrom}.phased.vcf.gz"
#     output:
#         "{wdirpop}/{dataset}.chromosome.{chrom}.pyrho.vcf.gz"
#     log:
#         "{wdirpop}/logs/{dataset}.make_subset_pyrho.{chrom}.log"
#     conda:
#         "envs/vcftools.yaml"
#     shell:
#         """
#         shuf -n {config[subset]} --random-source=<(yes {config[seed]}) {wdirpop}/poplist > {wdirpop}/subsetpyrho
#         vcftools --gzvcf {wdirpop}/{dataset}.chromosome.{chrom}.phased.vcf.gz --out {wdirpop}/out --recode --keep {wdirpop}/subsetpyrho --maf config[maf] --max-missing config[maxmissing]
#         mv {wdirpop}/out.recode.vcf {wdirpop}/{dataset}.chromosome.{chrom}.pyrho.vcf
#         bgzip -f {wdirpop}/{dataset}.chromosome.{chrom}.pyrho.vcf
#         """


# rule make_table:
#     """
#     Generate demography-aware look up tables
#     """
#     input:
#         "{wdirpop}/{dataset}.chromosome.{chrom}.pyrho.vcf.gz"
#     output:
#         "{wdirpop}/pyrho/{dataset}.lookuptable.{chrom}.hdf"
#     log:
#         "{wdirpop}/logs/{dataset}.make_table.{chrom}.log"
#     conda:
#         "envs/jq.yaml"
#     shell:
#         """
#         # Extract the coalescent sizes and times from the model of SMC++
#         N0=$( jq '.model.N0' {wdirpop}/smc/{dataset}.model.final.json)
#         coal_sizes=$( jq '.model.y' {wdirpop}/smc/{dataset}.model.final.json | tr -d '[:space:][]')
#         coal_sizes=$(echo $coal_sizes | awk '{{split($0, temp, ","); for(i=1; i < length(temp)+1; i++) {{a=exp(temp[i]); print a}}}}')
#         coal_sizes=$(echo $coal_sizes | tr -s '[:space:]' ',')
#         coal_sizes=${{coal_sizes%?}}
#         coal_sizes="$N0","$coal_sizes"
#         coal_times=$( jq '.model.knots' {wdirpop}/smc/{dataset}.model.final.json | tr -d '[:space:][]')
#         n=$(zcat {wdirpop}/{dataset}.chromosome.{chrom}.pyrho.vcf.gz | grep ^#CHROM | awk '{{print NF-9}}')
#         n=$((2*$n))
#         N=$((2*$n))
#         singularity exec --bind $PWD:/mnt pyrho.sif pyrho make_table --samplesize $n --approx --moran_pop_size $N --numthreads {config[cores]} --mu {config[mu]} --outfile /mnt/{wdirpop}/pyrho/{dataset}.lookuptable.{chrom}.hdf  --smcpp_file /mnt/{wdirpop}/smc/plot_chromosome.csv
#         """


# rule hyperparam:
#     """
#     Search the best hyperparameters for the dataset
#     <ploidy> should be set to 1 if using phased data and 2 for unphased genotype data. Ploidies other than 1 or 2 are not currently supported
#     Select the best pair of parameters
#     Select windowsize with the best LogL2
#     then best bpen, if there is differences between LogL2
#     Otherwise take bpen=50
#     """
#     input:
#         "{wdirpop}/pyrho/{dataset}.lookuptable.{chrom}.hdf"
#     output:
#         "{wdirpop}/pyrho/{dataset}.hyperparam.{chrom}"
#     log:
#         "{wdirpop}/logs/{dataset}.hyperparam.{chrom}"
#     shell:
#         """
#         #N0=$( jq '.model.N0' {wdirpop}/smc/{dataset}.model.final.json)
#         #coal_sizes=$( jq '.model.y' {wdirpop}/smc/{dataset}.model.final.json | tr -d '[:space:][]')
#         #coal_sizes=$(echo $coal_sizes | awk '{{split($0, temp, ","); for(i=1; i < length(temp)+1; i++) {{a=exp(temp[i]); print a}}}}')
#         #coal_sizes=$(echo $coal_sizes | tr -s '[:space:]' ',')
#         #coal_sizes=${{coal_sizes%?}}
#         #coal_sizes="$N0","$coal_sizes"
#         #coal_times=$( jq '.model.knots' {wdirpop}/smc/{dataset}.model.final.json | tr -d '[:space:][]')
#         n=$(zcat {wdirpop}/{dataset}.chromosome.{chrom}.pyrho.vcf.gz | grep ^#CHROM | awk '{{print NF-9}}')
#         n=$((2*$n))
#         singularity exec --bind $PWD:/mnt pyrho.sif pyrho hyperparam --samplesize $n --blockpenalty {config[bpen]} --windowsize {config[windowsize]} --num_sims {config[num_sims]} --numthreads {config[cores]} --tablefile /mnt/{input} --mu {config[mu]} --ploidy {config[ploidy]} --smcpp_file /mnt/{wdirpop}/smc/plot_chromosome.csv --outfile /mnt/{output} --logfile /mnt/{wdirpop}/pyrho/hyperparam.{chrom}.log
#         """


# rule selectparam:
#     """
#     Select the best hyperparameters
#     """
#     input:
#         "{wdirpop}/pyrho/{dataset}.hyperparam.{chrom}"
#     output:
#         "{wdirpop}/pyrho/{dataset}.windowsize.{chrom}",
#         "{wdirpop}/pyrho/{dataset}.bpen.{chrom}"
#     log:
#         "{wdirpop}/logs/{dataset}.selectparam.{chrom}"
#     conda:
#         "envs/Renv.yaml"
#     shell:
#         """
#         Rscript hyperparam_selection.R {dataset} {chrom} {wdirpop}
#         """

# rule optimize:
#     """
#     Estimate fine-scale recombination rates with Pyrho
#     """
#     input:
#         "{wdirpop}/pyrho/{dataset}.windowsize.{chrom}",
#         "{wdirpop}/pyrho/{dataset}.bpen.{chrom}"
#     output:
#         "{wdirpop}/pyrho/{dataset}.optimize.{chrom}"
#     log:
#         "{wdirpop}/logs/{dataset}.optimize.{chrom}"
#     shell:
#         """
#         windowsize=$(cat {wdirpop}/pyrho/{dataset}.windowsize.{chrom})
#         bpen=$(cat {wdirpop}/pyrho/{dataset}.bpen.{chrom})
#         windowsize=50
#         bpen=50
#         singularity exec --bind $PWD:/mnt pyrho.sif pyrho optimize --vcffile /mnt/{wdirpop}/{dataset}.chromosome.{chrom}.pyrho.vcf.gz --windowsize $windowsize --blockpenalty $bpen --tablefile /mnt/{wdirpop}/pyrho/{dataset}.lookuptable.{chrom}.hdf --ploidy {config[ploidy]} --outfile /mnt/{output} --numthreads {config[cores]}
#         """


# rule compute_r2:
#     """
#     Compute the population-scaled recombination rate as infered by Pyrho
#     """
#     input:
#         "{wdirpop}/pyrho/{dataset}.optimize.{chrom}"
#     output:
#         "{wdirpop}/pyrho/{dataset}.r2.{chrom}"
#     log:
#         "{wdirpop}/logs/{dataset}.r2.{chrom}"
#     shell:
#         """
#         n=$(zcat {wdirpop}/{dataset}.chromosome.{chrom}.pyrho.vcf.gz | grep ^#CHROM | awk '{{print NF-9}}')
#         singularity exec --bind $PWD:/mnt pyrho.sif pyrho compute_r2 --tablefile /mnt/{wdirpop}/pyrho/{dataset}.lookuptable.{chrom}.hdf --samplesize $n --quantiles 0.025,0.25,0.5,0.75,0.975 --compute_mean > {wdirpop}/pyrho/{dataset}.r2.{chrom}
#         """


rule pseudodiploid:
    """
    Make pseudodiploids (phased haplotypes) to take into account homozygotes in high-selfing rates species
    """
    input:
        "{wdirpop}/{dataset}.chromosome.{chrom}.phased.vcf.gz"
    output:
        "{wdirpop}/{dataset}.chromosome.{chrom}.pseudodiploid.vcf.gz"
    log:
        "{wdirpop}/logs/{dataset}.chromosome.{chrom}.pseudodiploid.log"
    conda:
        "envs/Renv.yaml"
    shell:
        """
        if [ {config[pseudodiploid]} -eq 1 ]; then
            Rscript pseudodiploids.R {wdirpop} {chrom}
        else
            cp {input} {output}
        fi
	"""


rule gzpseudodiploid:
    """
    Make pseudodiploids (phased haplotypes) to take into account homozygotes in high-selfing rates species
    """
    input:
        "{wdirpop}/{dataset}.chromosome.{chrom}.pseudodiploid.vcf.gz"
    output:
        "{wdirpop}/{dataset}.chromosome.{chrom}.pseudodiploid.vcf.gz.csi"
    log:
        "{wdirpop}/logs/{dataset}.chromosome.{chrom}.gzpseudodiploid.log"
    conda:
        "envs/vcftools.yaml"
    shell:
        """
        gunzip {wdirpop}/{dataset}.chromosome.{chrom}.pseudodiploid.vcf.gz
        bgzip {wdirpop}/{dataset}.chromosome.{chrom}.pseudodiploid.vcf
        tabix -p vcf {wdirpop}/{dataset}.chromosome.{chrom}.pseudodiploid.vcf.gz --csi
        """


# TODO A step to reduce sample size to a random subset
# Memory issues when sample size is high
rule subset_ldhat:
    """
    Subset a random sample of individuals
    """
    input:
        "{wdirpop}/{dataset}.chromosome.{chrom}.pseudodiploid.vcf.gz"
    output:
        "{wdirpop}/{dataset}.chromosome.{chrom}.ldhat.vcf.gz"
    log:
        "{wdirpop}/logs/{dataset}.make_subset_ldhat.{chrom}.log"
    conda:
        "envs/vcftools.yaml"
    shell:
        """
        #zcat {wdirpop}/{dataset}.chromosome.{chrom}.pseudodiploid.vcf.gz | head -n 10000 | grep "#CHROM" | tr "\\t" "\\n" | tail -n +10 > {wdirpop}/poplistrandomsample.{chrom}
        #shuf -n {config[subset]} --random-source=<(yes {config[seed]}) -o {wdirpop}/subsetldhat.{chrom} {wdirpop}/poplistrandomsample.{chrom}
        # Set a random seed
	RANDOM=42
	vcftools --gzvcf {wdirpop}/{dataset}.chromosome.{chrom}.pseudodiploid.vcf.gz --out {wdirpop}/out --recode --max-indv {config[subset]} --maf {config[maf]} --max-missing {config[maxmissing]}
        mv {wdirpop}/out.recode.vcf {wdirpop}/{dataset}.chromosome.{chrom}.ldhat.vcf
	bgzip -f {wdirpop}/{dataset}.chromosome.{chrom}.ldhat.vcf
        """

# Use LDHat/LDhot to estimate fine scale recombination rate and detect hotspots
rule LDpop:
    """
    Generate a demography-aware look-up table
    """
    input:
        "{wdirpop}/{dataset}.chromosome.{chrom}.ldhat.vcf.gz"
    output:
        "{wdirpop}/ldhat/{dataset}.ldpop.{chrom}"
    log:
        "{wdirpop}/logs/{dataset}.ldpop.{chrom}.log"
    conda:
        "envs/jq.yaml"
    shell:
        """
        N0=$( jq '.model.N0' {wdirpop}/smc/{dataset}.model.final.json)
        coal_sizes=$( jq '.model.y' {wdirpop}/smc/{dataset}.model.final.json | tr -d '[:space:][]')
        coal_sizes=$(echo $coal_sizes | awk '{{split($0, temp, ","); for(i=1; i < length(temp)+1; i++) {{a=exp(temp[i]); print a}}}}')
        coal_sizes=$(echo $coal_sizes | tr -s '[:space:]' ',')
        coal_sizes=${{coal_sizes%?}}
        coal_sizes="$N0","$coal_sizes"
        coal_times=$( jq '.model.knots' {wdirpop}/smc/{dataset}.model.final.json | tr -d '[:space:][]')
        n=$(zcat {wdirpop}/{dataset}.chromosome.{chrom}.ldhat.vcf.gz | grep ^#CHROM | awk '{{print NF-9}}')
        n=$((2*$n))
	echo $n
        singularity exec --bind $PWD:/mnt pyrho.sif python3 /ldpop/run/ldtable.py -n $n -th {config[theta]} -s $coal_sizes -t $coal_times -rh 101,100 --approx --cores {config[cores]} --log . > {wdirpop}/ldhat/{dataset}.ldpop.{chrom}
        """

rule convert:
    """
    Produce input files for LDhat
    Generate sites.txt and locs.txt with vcftools
    --ldhat
    --ldhat-geno
    These options output data in LDhat format. This option requires the "--chr" filter option to also be used. The first option outputs phased data only, and therefore also implies "--phased" be used, leading to unphased individuals and genotypes being excluded. The second option treats all of the data as unphased, and therefore outputs LDhat files in genotype/unphased format. Two output files are generated with the suffixes ".ldhat.sites" and ".ldhat.locs", which correspond to the LDhat "sites" and "locs" input files respectively.
    """
    input:
        ldpop = "{wdirpop}/ldhat/{dataset}.ldpop.{chrom}",
        vcf = "{wdirpop}/{dataset}.chromosome.{chrom}.ldhat.vcf.gz"
    output:
        "{wdirpop}/ldhat/{dataset}.{chrom}.ldhat.sites",
        "{wdirpop}/ldhat/{dataset}.{chrom}.ldhat.locs"
    log:
        "{wdirpop}/logs/{dataset}.ldhatconvert.{chrom}.log"
    conda:
        "envs/vcftools.yaml"
    shell:
        """
        vcftools --gzvcf {input.vcf} --chr {chrom} --ldhat --out {wdirpop}/ldhat/{dataset}.{chrom}
        """

rule interval:
    """
    Estimate a recombination landscape with LDhat
    """
    input:
        "{wdirpop}/ldhat/{dataset}.{chrom}.ldhat.sites",
        "{wdirpop}/ldhat/{dataset}.{chrom}.ldhat.locs"
    output:
        "{wdirpop}/ldhat/{dataset}.{chrom}.bpen{bpen}.new_lk.txt",
        "{wdirpop}/ldhat/{dataset}.{chrom}.bpen{bpen}.bounds.txt",
        "{wdirpop}/ldhat/{dataset}.{chrom}.bpen{bpen}.rates.txt"
    log:
        "{wdirpop}/logs/{dataset}.ldhatinterval.{chrom}.bpen{bpen}.log"
    shell:
        """
        iter={config[interval.iter]}
        samp={config[interval.samp]}
        bpen={config[interval.bpen]}
        singularity exec --bind $PWD:/mnt ldhat.sif /LDhat/interval -seq /mnt/{wdirpop}/ldhat/{dataset}.{chrom}.ldhat.sites -loc /mnt/{wdirpop}/ldhat/{dataset}.{chrom}.ldhat.locs -lk /mnt/{wdirpop}/ldhat/{dataset}.ldpop.{chrom} -its $iter -bpen $bpen -samp $samp -prefix /mnt/{wdirpop}/ldhat/{dataset}.{chrom}.bpen{config[interval.bpen]}.
        """


rule rhomap:
    """
    Estimate a recombination landscape with recombination hotspots
    on a background of low rate variation
    """


rule stat:
    """
    Compute statistics on interval
    """
    input:
        "{wdirpop}/ldhat/{dataset}.{chrom}.bpen{bpen}.new_lk.txt",
        "{wdirpop}/ldhat/{dataset}.{chrom}.bpen{bpen}.bounds.txt",
        "{wdirpop}/ldhat/{dataset}.{chrom}.bpen{bpen}.rates.txt"
    output:
        "{wdirpop}/ldhat/{dataset}.{chrom}.bpen{bpen}.res.txt"
    log:
        "{wdirpop}/logs/{dataset}.ldhatstat.{chrom}.bpen{bpen}.log"
    shell:
        """
        burn={config[ldhat.burn]}
        singularity exec --bind $PWD:/mnt ldhat.sif /LDhat/stat -input /mnt/{wdirpop}/ldhat/{dataset}.{chrom}.bpen{config[interval.bpen]}.rates.txt -burn $burn -loc /mnt/{wdirpop}/ldhat/{dataset}.{chrom}.ldhat.locs -prefix /mnt/{wdirpop}/ldhat/{dataset}.{chrom}.bpen{config[interval.bpen]}.
        """


rule LDhot:
    """
    Infer recombination hotspots
    LDhot
    """
    input:
        "{wdirpop}/ldhat/{dataset}.{chrom}.bpen{bpen}.res.txt"
    output:
        "{wdirpop}/ldhot/{dataset}.{chrom}.bpen{bpen}.hotspots.txt"
    log:
        "{wdirpop}/logs/{dataset}.{chrom}.bpen{bpen}.ldhot.log"
    shell:
        """
        nsim={config[ldhot.nsim]}
        singularity exec --bind $PWD:/mnt ldhat.sif /LDhot/ldhot --seq /mnt/{wdirpop}/ldhat/{dataset}.{chrom}.ldhat.sites --loc /mnt/{wdirpop}/ldhat/{dataset}.{chrom}.ldhat.locs --lk /mnt/{wdirpop}/ldhat/{dataset}.ldpop.{chrom} --res /mnt/{input} --nsim 100 --out /mnt/{wdirpop}/ldhot/{dataset}.{chrom}.bpen{config[interval.bpen]}
        # Summarize the results
        singularity exec --bind $PWD:/mnt ldhat.sif /LDhot/ldhot_summary --res /mnt/{input} --hot /mnt/{wdirpop}/ldhot/{dataset}.{chrom}.bpen{config[interval.bpen]}.hotspots.txt --out /mnt/{wdirpop}/ldhot/{dataset}.{chrom}.bpen{config[interval.bpen]}
        """
