"""
Pipeline to estimate fine-scale recombination maps from polymorphism data
"""

"""
Configuration of the analysis
i.e. dataset, name of chromosome, population to sample
"""
configfile: "config.yaml"

"""In addition to the configfile statement, config values can be overwritten via the command line"""
dataset=config["dataset"] # Name of your dataset directory and prefix of your vcf file
chrom=config["chrom"] # Name of the chromosome to analyse in your 'sample' dataset
K=config["K"]
pop=config["pop"]


wdir=config["workingdir"] + dataset

wildcard_constraints:
    wdir=wdir,
    dataset=dataset,
    chrom=chrom,
    K=K,
    pop=pop


rule all:
    """
    One ring to rule them all"
    """
    input:
        target = expand("{wdir}/{dataset}.{chrom}.hotspots.txt",wdir=wdir,dataset=dataset,chrom=chrom),
    shell:
        "echo 'Finished'"


rule poplist:
    """
    Sample individuals from a given genetic cluster (number of K to retain and name of the cluster) specified in config.yaml
    """
    output:
        "{wdir}/poplist"
    log:
        "{wdir}/logs/poplist.log"
    conda:
        "envs/vcftools.yaml"
    shell:
        """
        perl -ane '$r = 0; for my $i (1 .. $#F) {{$r = $i if $F[$i] > $F[$r];}} print $r + 1, " ";' < {wdir}/structure/faststructure.{K}.meanQ > {wdir}/structure/cluster.tmp
        sed -i -e "s/\\s\\+/\\n/g" {wdir}/structure/cluster.tmp
        paste {wdir}/indlist {wdir}/structure/cluster.tmp > {wdir}/structure/cluster.tmp2
        paste {wdir}/structure/cluster.tmp2 {wdir}/structure/faststructure.{K}.meanQ > {wdir}/structure/cluster
        awk -F' ' '{{if($2=={pop}) print $1}}' {wdir}/structure/cluster > {wdir}/poplist
        """


rule sampling_pop:
    """
    A first step to trim every 'sample' dataset to the same quality criteria
    """
    input:
        "{wdir}/poplist"
    output:
    	"{wdir}/{dataset}.pop.vcf.gz"
    log:
        "{wdir}/logs/{dataset}.sampling_pop.log"
    conda:
        "envs/vcftools.yaml"
    shell:
        """
        vcftools --gzvcf {wdir}/{dataset}.trimmed.vcf.gz --out {wdir}/out --recode --keep {wdir}/poplist --maf config[maf] --max-missing config[maxmissing]
        mv {wdir}/out.recode.vcf {wdir}/{dataset}.pop.vcf
        bgzip -f {wdir}/{dataset}.pop.vcf
        rm {wdir}/out.log
        """


# TODO Problem: vcf must be diploid in pyrho
# "ploidy set to 1 if using phased data and 2 for unphased genotype data"
# https://knausb.github.io/vcfR_documentation/dip_to_hap.html
# rule haploid:
#     """
#     Haploidization of diploid vcfs
#     to account for high selfing rates and homozygosity
#     Make pseudo-diploids for compatibility with pyrho
#     """
#     input:
#         "{wdir}/{dataset}.chromosome.{chrom}.phased.vcf.gz"
#     output:
#         "{wdir}/{dataset}.chromosome.{chrom}.haploid.vcf.gz"
#     log:
#         "{wdir}/logs/{dataset}.chromosome.{chrom}.haploid.log"
#     conda:
#         "envs/vcftools.yaml"
#     shell:
#         """
#         # Create vcf containing pseudo-diploid individual
#         zcat {input} | sed 's/|[0-2]//g' > {wdir}/{dataset}.chromosome.{chrom}.haploid.vcf
#         bgzip {wdir}/{dataset}.chromosome.{chrom}.haploid.vcf
#         tabix -p vcf {output} --csi
#         """



rule effective_size:
    """
    Estimate population effective size
    """
    input:
        "{wdir}/{dataset}.pop.vcf.gz"
    output:
        "{wdir}/statistics/{dataset}.effective_size"
    log:
        "{wdir}/logs/{dataset}.effective_size.log"
    conda:
        "envs/vcftools.yaml"
    params:
        windowsize=100000
    shell:
        """
        vcftools --gzvcf {input} --window-pi {params.windowsize} --out {wdir}/statistics/vcftools
        eff_size=$(awk '{{sum+=$5}} END {{sum=(sum/NR)/(4*{config[mu]}); print sum}}' {wdir}/statistics/vcftools.windowed.pi)
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
        "{wdir}/statistics/{dataset}.effective_size"
    output:
        "{wdir}/smc/{dataset}.model.final.json"
    log:
        "{wdir}/logs/{dataset}.demography.log"
    conda:
        "envs/vcftools.yaml"
    shell:
        """
        tabix -f -p vcf {wdir}/{dataset}.pop.vcf.gz --csi
        shuf -n {config[subset]} --random-source=<(yes {config[seed]}) {wdir}/poplist > {wdir}/subsetpop
        individuals=$(cat {wdir}/subsetpop | tr '\n' ',' | sed "s/,$//g")
        chromosomes=$(zcat {wdir}/{dataset}.pop.vcf.gz | awk '{{ print $1 }}' | sort | uniq | grep -v '^#')
        for c in $chromosomes; do
        singularity exec --bind $PWD:/mnt smcpp.simg smc++ vcf2smc --ignore-missing /mnt/{wdir}/{dataset}.pop.vcf.gz /mnt/{wdir}/smc/vcf2smc.$c $c pop1:$individuals
        done
        paths=$(for i in $chromosomes; do echo /mnt/{wdir}/smc/vcf2smc.$i; done | tr '\n' ' ')
        singularity exec --bind $PWD:/mnt smcpp.simg smc++ estimate -o /mnt/{wdir}/smc/ {config[mu]} $paths --cores {config[cores]}
        mv {wdir}/smc/model.final.json {wdir}/smc/{dataset}.model.final.json
        singularity exec --bind $PWD:/mnt smcpp.simg smc++ plot /mnt/{wdir}/smc/plot_chromosome.pdf /mnt/{wdir}/smc/{dataset}.model.final.json -c
        singularity exec --bind $PWD:/mnt smcpp.simg smc++ posterior /mnt/{wdir}/smc/{dataset}.model.final.json /mnt/{wdir}/smc/{dataset}.posterior.smc $paths
        mv ~/iterate.dat {wdir}/smc/{dataset}.{chrom}.iterate.dat
        """


# From now, analyses are performed on individual chromosomes
rule split_chromosome:
    """
    Split the entire vcf into one vcf per chromosome
    Deal with missing data: remove SNPs with missing data for more than 90% of individuals
    """
    input:
        "{wdir}/smc/{dataset}.model.final.json"
    output:
        "{wdir}/{dataset}.chromosome.{chrom}.vcf.gz"
    log:
        "{wdir}/logs/{dataset}.split_chromosome.{chrom}.log"
    conda:
        "envs/vcftools.yaml"
    shell:
        """
        vcftools --gzvcf {wdir}/{dataset}.pop.vcf.gz --out {wdir}/out --recode --chr {chrom} --maf config[maf] --max-missing config[maxmissing]
        mv {wdir}/out.recode.vcf {wdir}/{dataset}.chromosome.{chrom}.vcf
        bgzip -f {wdir}/{dataset}.chromosome.{chrom}.vcf
        rm {wdir}/out.log
        """

rule phasing_vcf:
    """
    Phase the vcf with Shapeit2
    Unphased haplotypes can overestimate recombination rate, hence false postive hotspots may be detected
    https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html#output
    """
    input:
        "{wdir}/{dataset}.chromosome.{chrom}.vcf.gz"
    output:
        "{wdir}/{dataset}.chromosome.{chrom}.phased.vcf.gz"
    log:
        "{wdir}/logs/{dataset}.chromosome.{chrom}.phasing_vcf.log"
    conda:
        "envs/shapeit.yaml"
    shell:
        """
        shapeit --input-vcf {wdir}/{dataset}.chromosome.{chrom}.vcf.gz --output-max {wdir}/{dataset}.phased.chromosome.{chrom} --effective-size $(cat {wdir}/statistics/{dataset}.effective_size) --window 1 --thread {config[cores]} --output-log {wdir}/logs/{dataset}.chromosome.{chrom}.shapeit.log --force
        shapeit -convert --input-haps {wdir}/{dataset}.phased.chromosome.{chrom} --output-vcf {wdir}/{dataset}.chromosome.{chrom}.phased.vcf --output-log {wdir}/logs/{dataset}.chromosome.{chrom}.shapeit.convert.log
        # replace header in vcf to keep information of contig length
        zcat {wdir}/{dataset}.vcf.gz | grep '^#' > {wdir}/newheader
        cat {wdir}/newheader | grep -v '^#CHROM' > {wdir}/newheader2
        cat {wdir}/{dataset}.chromosome.{chrom}.phased.vcf | grep '^#CHROM' > {wdir}/colnames
        cat {wdir}/{dataset}.chromosome.{chrom}.phased.vcf | grep -v '^#' > {wdir}/newvcf
        cat {wdir}/newheader2 {wdir}/colnames {wdir}/newvcf > {wdir}/{dataset}.chromosome.{chrom}.phased.vcf
        bgzip -f {wdir}/{dataset}.chromosome.{chrom}.phased.vcf
        rm {wdir}/newheader {wdir}/newheader2 {wdir}/colnames {wdir}/newvcf
        """


# TODO A step to reduce sample size to a random subset
# Memory issues when sample size is high
rule subset_pyrho:
    """
    Subset a random sample of individuals
    """
    input:
        "{wdir}/{dataset}.chromosome.{chrom}.phased.vcf.gz"
    output:
        "{wdir}/{dataset}.chromosome.{chrom}.pyrho.vcf.gz"
    log:
        "{wdir}/logs/{dataset}.make_subset_pyrho.{chrom}.log"
    conda:
        "envs/vcftools.yaml"
    shell:
        """
        shuf -n {config[subset]} --random-source=<(yes {config[seed]}) {wdir}/poplist > {wdir}/subsetpyrho
        vcftools --gzvcf {wdir}/{dataset}.chromosome.{chrom}.phased.vcf.gz --out {wdir}/out --recode --keep {wdir}/subsetpyrho --maf config[maf] --max-missing config[maxmissing]
        mv {wdir}/out.recode.vcf {wdir}/{dataset}.chromosome.{chrom}.pyrho.vcf
        bgzip -f {wdir}/{dataset}.chromosome.{chrom}.pyrho.vcf
        rm {wdir}/out.log
        """


rule make_table:
    """
    Generate demography-aware look up tables
    """
    input:
        "{wdir}/{dataset}.chromosome.{chrom}.pyrho.vcf.gz"
    output:
        "{wdir}/pyrho/{dataset}.lookuptable.{chrom}.hdf"
    log:
        "{wdir}/{dataset}.make_table.{chrom}.log"
    conda:
        "envs/jq.yaml"
    shell:
        """
        # Extract the coalescent sizes and times from the model of SMC++
        N0=$( jq '.model.N0' {wdir}/smc/{dataset}.model.final.json)
        coal_sizes=$( jq '.model.y' {wdir}/smc/{dataset}.model.final.json | tr -d '[:space:][]')
        coal_sizes=$(echo $coal_sizes | awk '{{split($0, temp, ","); for(i=1; i < length(temp)+1; i++) {{a=exp(temp[i]); print a}}}}')
        coal_sizes=$(echo $coal_sizes | tr -s '[:space:]' ',')
        coal_sizes=${{coal_sizes%?}}
        coal_sizes="$N0","$coal_sizes"
        coal_times=$( jq '.model.knots' {wdir}/smc/{dataset}.model.final.json | tr -d '[:space:][]')
        n=$(zcat {wdir}/{dataset}.chromosome.{chrom}.pyrho.vcf.gz | grep ^#CHROM | awk '{{print NF-9}}')
        n=$((2*$n))
        N=$((2*$n))
        singularity exec --bind $PWD:/mnt pyrho.simg pyrho make_table --samplesize $n --approx --moran_pop_size $N --numthreads {config[cores]} --mu {config[mu]} --outfile /mnt/{wdir}/pyrho/{dataset}.lookuptable.{chrom}.hdf  --smcpp_file /mnt/{wdir}/smc/plot_chromosome.csv
        """


rule hyperparam:
    """
    Search the best hyperparameters for the dataset
    <ploidy> should be set to 1 if using phased data and 2 for unphased genotype data. Ploidies other than 1 or 2 are not currently supported
    Select the best pair of parameters
    Select windowsize with the best LogL2
    then best bpen, if there is differences between LogL2
    Otherwise take bpen=50
    """
    input:
        "{wdir}/pyrho/{dataset}.lookuptable.{chrom}.hdf"
    output:
        "{wdir}/pyrho/{dataset}.hyperparam.{chrom}"
    log:
        "{wdir}/logs/{dataset}.hyperparam.{chrom}"
    shell:
        """
        #N0=$( jq '.model.N0' {wdir}/smc/{dataset}.model.final.json)
        #coal_sizes=$( jq '.model.y' {wdir}/smc/{dataset}.model.final.json | tr -d '[:space:][]')
        #coal_sizes=$(echo $coal_sizes | awk '{{split($0, temp, ","); for(i=1; i < length(temp)+1; i++) {{a=exp(temp[i]); print a}}}}')
        #coal_sizes=$(echo $coal_sizes | tr -s '[:space:]' ',')
        #coal_sizes=${{coal_sizes%?}}
        #coal_sizes="$N0","$coal_sizes"
        #coal_times=$( jq '.model.knots' {wdir}/smc/{dataset}.model.final.json | tr -d '[:space:][]')
        n=$(zcat {wdir}/{dataset}.chromosome.{chrom}.pyrho.vcf.gz | grep ^#CHROM | awk '{{print NF-9}}')
        n=$((2*$n))
        singularity exec --bind $PWD:/mnt pyrho.simg pyrho hyperparam --samplesize $n --blockpenalty {config[bpen]} --windowsize {config[windowsize]} --num_sims {config[num_sims]} --numthreads {config[cores]} --tablefile /mnt/{input} --mu {config[mu]} --ploidy {config[ploidy]} --smcpp_file /mnt/{wdir}/smc/plot_chromosome.csv --outfile /mnt/{output} --logfile /mnt/{wdir}/pyrho/hyperparam.{chrom}.log
        """


rule selectparam:
    """
    Select the best hyperparameters
    """
    input:
        "{wdir}/pyrho/{dataset}.hyperparam.{chrom}"
    output:
        "{wdir}/pyrho/{dataset}.windowsize.{chrom}",
        "{wdir}/pyrho/{dataset}.bpen.{chrom}"
    log:
        "{wdir}/logs/{dataset}.selectparam.{chrom}"
    conda:
        "envs/Renv.yaml"
    shell:
        """
        Rscript hyperparam_selection.R {dataset} {chrom}
        """

rule optimize:
    """
    Estimate fine-scale recombination rates with Pyrho
    """
    input:
        "{wdir}/pyrho/{dataset}.windowsize.{chrom}",
        "{wdir}/pyrho/{dataset}.bpen.{chrom}"
    output:
        "{wdir}/pyrho/{dataset}.optimize.{chrom}"
    log:
        "{wdir}/logs/{dataset}.optimize.{chrom}"
    shell:
        """
        windowsize=$(cat {wdir}/pyrho/{dataset}.windowsize.{chrom})
        bpen=$(cat {wdir}/pyrho/{dataset}.bpen.{chrom})
        windowsize=50
        bpen=50
        singularity exec --bind $PWD:/mnt pyrho.simg pyrho optimize --vcffile /mnt/{wdir}/{dataset}.chromosome.{chrom}.pyrho.vcf.gz --windowsize $windowsize --blockpenalty $bpen --tablefile /mnt/{wdir}/pyrho/{dataset}.lookuptable.{chrom}.hdf --ploidy {config[ploidy]} --outfile /mnt/{output} --numthreads {config[cores]}
        """


rule compute_r2:
    """
    Compute the population-scaled recombination rate as infered by Pyrho
    """
    input:
        "{wdir}/pyrho/{dataset}.optimize.{chrom}"
    output:
        "{wdir}/pyrho/{dataset}.r2.{chrom}"
    log:
        "{wdir}/logs/{dataset}.r2.{chrom}"
    shell:
        """
        n=$(zcat {wdir}/{dataset}.chromosome.{chrom}.pyrho.vcf.gz | grep ^#CHROM | awk '{{print NF-9}}')
        singularity exec --bind $PWD:/mnt pyrho.simg pyrho compute_r2 --tablefile /mnt/{wdir}/pyrho/{dataset}.lookuptable.{chrom}.hdf --samplesize $n --quantiles 0.025,0.25,0.5,0.75,0.975 --compute_mean > {wdir}/pyrho/{dataset}.r2.{chrom}
        """


# Use LDHat/LDhot to estimate fine scale recombination rate and detect hotspots
rule LDpop:
    """
    Generate a demography-aware look-up table
    """
    input:
        "{wdir}/pyrho/{dataset}.r2.{chrom}"
    output:
        "{wdir}/ldhat/{dataset}.ldpop.{chrom}"
    log:
        "{wdir}/logs/{dataset}.ldpop.{chrom}.log"
    conda:
        "envs/jq.yaml"
    shell:
        """
        N0=$( jq '.model.N0' {wdir}/smc/{dataset}.model.final.json)
        coal_sizes=$( jq '.model.y' {wdir}/smc/{dataset}.model.final.json | tr -d '[:space:][]')
        coal_sizes=$(echo $coal_sizes | awk '{{split($0, temp, ","); for(i=1; i < length(temp)+1; i++) {{a=exp(temp[i]); print a}}}}')
        coal_sizes=$(echo $coal_sizes | tr -s '[:space:]' ',')
        coal_sizes=${{coal_sizes%?}}
        coal_sizes="$N0","$coal_sizes"
        coal_times=$( jq '.model.knots' {wdir}/smc/{dataset}.model.final.json | tr -d '[:space:][]')
        n=$(zcat {wdir}/{dataset}.chromosome.{chrom}.pyrho.vcf.gz | grep ^#CHROM | awk '{{print NF-9}}')
        n=$((2*$n))
        singularity exec --bind $PWD:/mnt pyrho.simg python3 /ldpop/run/ldtable.py -n $n -th .001 -s $coal_sizes -t $coal_times -rh 101,100 --approx --cores {config[cores]} --log . > {wdir}/ldhat/{dataset}.ldpop.{chrom}
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
        ldpop = "{wdir}/ldhat/{dataset}.ldpop.{chrom}",
        vcf = "{wdir}/{dataset}.chromosome.{chrom}.pyrho.vcf.gz"
    output:
        "{wdir}/ldhat/{dataset}.{chrom}.ldhat.sites",
        "{wdir}/ldhat/{dataset}.{chrom}.ldhat.locs"
    log:
        "{wdir}/{dataset}.ldhatconvert.{chrom}.log"
    conda:
        "envs/vcftools.yaml"
    shell:
        """
        vcftools --gzvcf {wdir}/{dataset}.chromosome.{chrom}.pyhro.vcf.gz --chr {chrom} --ldhat-geno --out {wdir}/ldhat/{dataset}.{chrom}
        """

rule interval:
    """
    Estimate a recombination landscape with LDhat
    """
    input:
        "{wdir}/ldhat/{dataset}.{chrom}.ldhat.sites",
        "{wdir}/ldhat/{dataset}.{chrom}.ldhat.locs"
    output:
        "{wdir}/{dataset}.{chrom}.new_lk.txt",
        "{wdir}/{dataset}.{chrom}.bounds.txt",
        "{wdir}/{dataset}.{chrom}.rates.txt"
    log:
        "{wdir}/{dataset}.ldhatinterval.{chrom}.log"
    shell:
        """
        iter={config[interval.iter]}
        samp={config[interval.samp]}
        bpen={config[interval.bpen]}
        singularity exec --bind $PWD:/mnt ldhat.simg /LDhat/interval -seq /mnt/{wdir}/ldhat/{dataset}.{chrom}.ldhat.sites -loc /mnt/{wdir}/ldhat/{dataset}.{chrom}.ldhat.locs -lk /mnt/{wdir}/ldhat/{dataset}.ldpop.{chrom} -its $iter -bpen $bpen -samp $samp
        mv ~/new_lk.txt {wdir}/ldhat/{dataset}.{chrom}.new_lk.txt
        mv ~/bounds.txt {wdir}/ldhat/{dataset}.{chrom}.bounds.txt
        mv ~/rates.txt {wdir}/ldhat/{dataset}.{chrom}.rates.txt
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
        "{wdir}/{dataset}.{chrom}.new_lk.txt",
        "{wdir}/{dataset}.{chrom}.bounds.txt",
        "{wdir}/{dataset}.{chrom}.rates.txt"
    output:
        "{wdir}/{dataset}.{chrom}.res.txt"
    log:
        "{wdir}/{dataset}.ldhatstat.{chrom}.log"
    shell:
        """
        burn={config[ldhat.burn]}
        singularity exec --bind $PWD:/mnt ldhat.simg /LDhat/stat -input {wdir}/ldhat/{dataset}.{chrom}.rates.txt -burn $burn -loc {wdir}/ldhat/{dataset}.{chrom}.ldhat.locs  -prefix {wdir}/ldhat/{dataset}.{chrom}
        """


rule LDhot:
    """
    Infer recombination hotspots
    LDhot
    """
    input:
        "{wdir}/{dataset}.{chrom}.res.txt"
    output:
        "{wdir}/{dataset}.{chrom}.hotspots.txt"
    log:
        "{wdir}/logs/{dataset}.{chrom}.ldhot.log"
    shell:
        """
        nsim={config[ldhot.nsim]}
        singularity exec --bind $PWD:/mnt ldhat.simg /LDhot/ldhot --seq {wdir}/ldhat/{dataset}.{chrom}.ldhat.sites --loc {wdir}/ldhat/{dataset}.{chrom}.ldhat.locs --lk {wdir}/ldhat/{dataset}.ldpop.{chrom} --res {input} --nsim 100 --out {wdir}/ldhot/{dataset}.{chrom}
        # Summarize the results
        singularity exec --bind $PWD:/mnt ldhat.simg /LDhot/ldhot_summary --res {input} --hot {wdir}/{dataset}.{chrom}.hotspots.txt --out {wdir}/ldhot/{dataset}.{chrom}
        """
