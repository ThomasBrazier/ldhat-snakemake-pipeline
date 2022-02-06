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
        target = expand("{wdir}/statistics/{dataset}.effective_size.chromosome.{chrom}",wdir=wdir,dataset=dataset,chrom=chrom)
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
        vcf = "{wdir}/{dataset}.trimmed.vcf.gz",
        poplist = "{wdir}/poplist"
    output:
    	"{wdir}/{dataset}.pop.vcf.gz"
    log:
        "{wdir}/logs/{dataset}.sampling_pop.log"
    conda:
        "envs/vcftools.yaml"
    shell:
        """
        vcftools --gzvcf {input.vcf} --out {wdir}/out --recode --keep {wdir}/poplist --maf config[maf] --max-missing config[maxmissing]
        mv {wdir}/out.recode.vcf {wdir}/{dataset}.pop.vcf
        bgzip -f {wdir}/{dataset}.pop.vcf
        rm {wdir}/out.log
        """


# From now, analyses are performed on individual chromosomes
rule split_chromosome:
    """
    Split the entire vcf into one vcf per chromosome
    Deal with missing data: remove SNPs with missing data for more than 90% of individuals
    """
    input:
        "{wdir}/{dataset}.pop.vcf.gz"
    output:
        "{wdir}/{dataset}.chromosome.{chrom}.vcf.gz"
    log:
        "{wdir}/logs/{dataset}.split_chromosome.{chrom}.log"
    conda:
        "envs/vcftools.yaml"
    shell:
        """
        vcftools --gzvcf {input} --out {wdir}/out --recode --chr {chrom} --maf config[maf] --max-missing config[maxmissing]
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
        vcftools --gzvcf {input} --window-pi {params.windowsize} --out {wdir}/statistics/vcftools.{chrom} && eff_size=$(awk '{{sum+=$5}} END {{sum=(sum/NR)/(4*{config[mu]}); print sum}}') {wdir}/statistics/vcftools.{chrom}.windowed.pi && echo ${{eff_size/.*}} > {output}
        """


rule phasing_vcf:
    """
    Phase the vcf with Shapeit2
    Unphased haplotypes can overestimate recombination rate, hence false postive hotspots may be detected
    https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html#output
    """
    input:
        vcf = "{wdir}/{dataset}.chromosome.{chrom}.vcf.gz"
        effsize = "{wdir}/statistics/{dataset}.effective_size.{chrom}"
    output:
        vcf = "{wdir}/{dataset}.chromosome.{chrom}.phased.vcf.gz"
    log:
        "{wdir}/logs/{dataset}.phasing_vcf.chromosome.{chrom}.log"
    conda:
        "envs/vcftools.yaml"
        "envs/shapeit.yaml"
    shell:
        """
        shapeit --input-vcf {input.vcf} --output-max {wdir}/{dataset}.chromosome.{chrom}.phased --effective-size $(cat {input.effsize}) --window 1 --thread {config[cores]} --output-log {wdir}/logs/{dataset}.shapeit.chromosome.{chrom}.log --force
        shapeit -convert --input-haps {wdir}/{dataset}.chromosome.{chrom}.phased.haps --output-vcf {wdir}/{dataset}.chromosome.{chrom}.vcf --output-log {wdir}/logs/{dataset}.shapeit.convert.chromosome.{chrom}.log
        bgzip {wdir}/{dataset}.chromosome.{chrom}.vcf
        """



rule haploidization:
    """
    Create pseudo-diploid vcf but with haploid genotypes
    """
    input:
        "{wdir}/{dataset}.chromosome.{chrom}.phased.vcf.gz"
    output:
        "{wdir}/{dataset}.chromosome.{chrom}.haploid.vcf.gz"
    log:
        "{wdir}/logs/{dataset}.haploidization.{chrom}.log"
    conda:
        "envs/vcftools.yaml"
    shell:
        """
        zcat {input} | sed 's/|.\t.|/|/g' > {wdir}/{dataset}.chromosome.{chrom}.haploid.vcf
        sed -r -i '/^#CHROM/s/\S+//g' {wdir}/{dataset}.chromosome.{chrom}.haploid.vcf
        bgzip {wdir}/{dataset}.chromosome.{chrom}.haploid.vcf
        """

