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


wdir='data/' + dataset

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
        target = expand("{wdir}/pyrho/{dataset}.table.{chrom}.hdf",wdir=wdir,dataset=dataset,chrom=chrom)
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


# TODO Take care of run of homozygosity (mask): raise WARNING
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
        "{wdir}/{dataset}.pop.vcf.gz"
    output:
        "{wdir}/smc/{dataset}.model.final.json"
    log:
        "{wdir}/logs/{dataset}.demography.log"
    conda:
        "envs/vcftools.yaml"
    shell:
        """
        tabix -f -p vcf {input} --csi
        shuf -n 40 --random-source=<(yes 42) {wdir}/poplist > {wdir}/subsetpop
        individuals=$(cat {wdir}/subsetpop | tr '\n' ',' | sed "s/,$//g")
        chromosomes=$(zcat {input} | awk '{{ print $1 }}' | sort | uniq | grep -v '^#')
        for c in $chromosomes; do
        singularity exec --bind $PWD:/mnt smcpp.simg smc++ vcf2smc --ignore-missing /mnt/{input} /mnt/{wdir}/smc/vcf2smc.$c $c pop1:$individuals
        done
        paths=$(for i in $chromosomes; do echo /mnt/{wdir}/smc/vcf2smc.$i; done | tr '\n' ' ')
        singularity exec --bind $PWD:/mnt smcpp.simg smc++ estimate -o /mnt/{wdir}/smc/ {config[mu]} $paths --cores {config[cores]}
        mv {wdir}/smc/model.final.json {wdir}/smc/{dataset}.model.final.json
        singularity exec --bind $PWD:/mnt smcpp.simg smc++ plot /mnt/{wdir}/smc/plot_chromosome.pdf /mnt/{wdir}/smc/{dataset}.model.final.json -c
        singularity exec --bind $PWD:/mnt smcpp.simg smc++ posterior /mnt/{wdir}/smc/{dataset}.model.final.json /mnt/{wdir}/smc/{dataset}.posterior.smc $paths
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


rule effective_size:
    """
    Estimate population effective size
    """
    input:
        "{wdir}/{dataset}.chromosome.{chrom}.vcf.gz"
    output:
        "{wdir}/statistics/{dataset}.effective_size.chromosome.{chrom}"
    log:
        "{wdir}/logs/{dataset}.effective_size.{chrom}.log"
    conda:
        "envs/vcftools.yaml"
    params:
        windowsize=100000
    shell:
        """
        vcftools --gzvcf {input} --window-pi {params.windowsize} --out {wdir}/statistics/vcftools.{chrom}
        eff_size=$(awk '{{sum+=$5}} END {{sum=(sum/NR)/(4*{config[mu]}); print sum}}' {wdir}/statistics/vcftools.{chrom}.windowed.pi)
        echo ${{eff_size/.*}} > {output}
        """

rule phasing_vcf:
    """
    Phase the vcf with Shapeit2
    Unphased haplotypes can overestimate recombination rate, hence false postive hotspots may be detected
    https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html#output
    """
    input:
        "{wdir}/statistics/{dataset}.effective_size.chromosome.{chrom}"
    output:
        "{wdir}/{dataset}.chromosome.{chrom}.phased.vcf.gz"
    log:
        "{wdir}/logs/{dataset}.phasing_vcf.chromosome.{chrom}.log"
    conda:
        "envs/shapeit.yaml"
    shell:
        """
        shapeit --input-vcf {wdir}/{dataset}.chromosome.{chrom}.vcf.gz --output-max {wdir}/{dataset}.chromosome.{chrom}.phased --effective-size $(cat {input}) --window 1 --thread {config[cores]} --output-log {wdir}/logs/{dataset}.shapeit.chromosome.{chrom}.log --force
        shapeit -convert --input-haps {wdir}/{dataset}.chromosome.{chrom}.phased --output-vcf {wdir}/{dataset}.chromosome.{chrom}.phased.vcf --output-log {wdir}/logs/{dataset}.shapeit.convert.chromosome.{chrom}.log
        # replace header in vcf to keep information of contig length
        zcat {wdir}/{dataset}.vcf.gz | grep '^#' > {wdir}/newheader
        cat {wdir}/newheader | grep -v '^#CHROM' > {wdir}/newheader2
        cat {wdir}/{dataset}.chromosome.{chrom}.phased.vcf | grep '^#CHROM' > {wdir}/colnames
        cat {wdir}/{dataset}.chromosome.{chrom}.phased.vcf | grep -v '^#' > {wdir}/newvcf
        cat {wdir}/newheader2 {wdir}/colnames {wdir}/newvcf > {wdir}/{dataset}.chromosome.{chrom}.phased.vcf
        bgzip -f {wdir}/{dataset}.chromosome.{chrom}.phased.vcf
        rm {wdir}/newheader {wdir}/newheader2 {wdir}/colnames {wdir}/newvcf
        """

# TODO haploidization
# https://github.com/popgenmethods/smcpp/issues/82
# rule haploidization:
#     """
#     Create pseudo-diploid vcf but with haploid genotypes
#     Take first allele, regardless of phasing
#     """
#     input:
#         "{wdir}/{dataset}.chromosome.{chrom}.phased.vcf.gz"
#     output:
#         "{wdir}/{dataset}.chromosome.{chrom}.haploid.vcf.gz"
#     log:
#         "{wdir}/logs/{dataset}.haploidization.{chrom}.log"
#     conda:
#         "envs/vcftools.yaml"
#     shell:
#         """

#         zcat {input} | sed 's/|.\t.|/|/g' > {wdir}/{dataset}.chromosome.{chrom}.haploid.vcf
#         sed -r -i '/^#CHROM/s/\S+//g' {wdir}/{dataset}.chromosome.{chrom}.haploid.vcf
#         bgzip -f {wdir}/{dataset}.chromosome.{chrom}.haploid.vcf
#         """


rule make_table:
    """
    Generate demography-aware look up tables
    """
    input:
        "{wdir}/{dataset}.chromosome.{chrom}.phased.vcf.gz"
    output:
        "{wdir}/pyrho/{dataset}.table.{chrom}.hdf"
    log:
        "{wdir}/{dataset}.make_table.{chrom}.log"
    conda:
        "envs/pyrho.yaml"
    shell:
        """
        # Extract the coalescent sizes and times from the model of SMC++
        N0=$( jq '.model.N0' {wdir}/smc/{dataset}.model.final.json)
        coal_sizes=$( jq '.model.y' {wdir}/smc/{dataset}.model.final.json | tr -d '[:space:][]')
        coal_sizes=$(echo $coal_sizes | awk '{split($0, temp, ","); for(i=1; i < length(temp)+1; i++) {a=exp(temp[i]); print a}}')
        coal_sizes=$(echo $coal_sizes | tr -s '[:space:]' ',')
        coal_sizes=${coal_sizes%?}
        coal_sizes="$N0","$coal_sizes"
        coal_times=$( jq '.model.knots' {wdir}/smc/{dataset}.model.final.json | tr -d '[:space:][]')
        n=$(awk '{if ($1 == "#CHROM"){print NF-9; exit}}' {input})
        N=1.5*$n
        pyrho make_table --samplesize $n --approx --moran_pop_size $N \
        --numthreads config[cores] --mu config[mu] --outfile {output} \
        --popsizes $coal_sizes --epochtimes $coal_times
        """


rule hyperparam:
    """
    Search the best hyperparameters for the dataset
    """


rule optimize:
    """
    Estimate fine-scale recombination rates with Pyrho
    """



rule compute_r2:
    """
    Compute the population-scaled recombination rate as infered by Pyrho
    """
