"""
Pipeline to estimate fine-scale recombination maps from polymorphism data
LDhelmet is the central element
"""

"""In addition to the configfile statement, config values can be overwritten via the command line"""
dataset=config["dataset"] # Name of your dataset directory and prefix of your vcf file
pop=config["pop"] # Specify the population once you have defined it with genetic clustering
chrom=config["chrom"] # Name of the chromosome to analyse in your 'sample' dataset


"""
Configuration of the analysis
i.e. dataset, name of chromosome, population to sample
"""
configfile: "{dataset}/config.yaml"


rule all:
    """
    One ring to rule them all"
    """
    input:
        ""


# The global dataset is trimmed for SNPs
# and individuals if a 'samplelist' file is provided
# Otherwise all individuals are kept
rule trimming_vcf:
    """
    A first step to trim every 'sample' dataset to the same quality criteria
    """
    input:
        "{sample}/{sample}.vcf.gz"
    output:
    	"{sample}/{sample}.trimmed.vcf.gz"
    log:
        "{sample}/logs/trimming_vcf.log"
    conda:
        "envs/vcftools.yaml"
    params:
        maf=config["maf"]
	    missing=config["missing"]
	    maxmissing=config["maxmissing"]
    shell:
        "vcftools -gzvcf {input} --out {output} --recode --maf {maf} --max-missing {maxmissing} && mv {sample}/{sample}.recode.vcf {sample}/{sample}.vcf && bgzip {sample}/{sample}.vcf"


rule population_struture:
    """
    Detects population structure and infer genetic clusters
    Produces diagnostic plots to select the right population
    """

# checkpoints trigger a re-evaluation of the Snakefile rules in light of new information
checkpoint choose_population:
    """
    Choose a population from the genetic cluster analysis
    """
    input:
        "{sample}/clusters"
    output:
        "{sample}/pop"
    log:
        "{sample}/logs/choose_population.log"
    shell:
        # if 'pop', does not exist yet, ask user which population to sample and save the result in a file 'pop'
        # 'pop' contains only a single cluster name
        "[ ! -z '${pop}' ] && echo ${pop} > {pop} || echo 'Set a population to sample in the 'pop' file' "


# input function for the rule sample_population
def is_popexist(wildcards):
# decision based on the existence of a non empty 'pop' file
    with checkpoints.choose_population.get(sample=wildcards.sample).output[0].open() as f:
        if len(f.read()) > 0:
            return "{sample}/pop"


# Sample population is conditional on the 'pop' file that identifies the population to sample
# Or the population name passed in command-line argument
rule sample_population:
    """
    Sample a population in the vcf, depending on results of population structure
    Needs intervention of the user, manual selection of the best population
    Return a list of individuals to sample and a vcf subset
    """
    input:
        pop = is_popexist
        vcf = "{sample}/{sample}.trimmed.vcf.gz"
        clusters = "{sample}/{sample}.clusters"
    output:
        poplist = "{sample}/poplist" # A list of individuals to sample in a given population
        vcf = "{sample}/{sample}.population.vcf.gz"
    log:
        "{sample}/logs/sample_population.log"
    shell:
        # Select individuals in genetic clusters for the chosen cluster
        "awk '$1 ~ /${input.pop}/ { print $0 }' {input.clusters} > {output.poplist}"
        # Subset the vcf
        "vcftools -gzvcf {input.vcf} --out {sample}/{sample} --recode --maf {maf} --max-missing {maxmissing} && mv {sample}/{sample}.recode.vcf {sample}/{sample}.population.vcf && bgzip {sample}/{sample}.population.vcf"


# From now, analyses are performed on individual chromosomes
rule split_by_chromosome:
    """
    Split the entire vcf into one vcf per chromosome
    Deal with missing data: remove SNPs with missing data for more than 90% of individuals
    """
    input:
        "{sample}/{sample}.population.vcf.gz"
    output:
        "{sample}/{sample}.chromosome.{chr}"
    log:
        "{sample}/logs/split_by_chromosome.{chr}.log"
    conda:
        "envs/vcftools.yaml"
    params:
        maxmissing=config["maxmissing"]
    shell:
        "vcftools --gzvcf {input} --out {output} --recode --chr {chr} --max-missing {maxmissing} && mv {output}.recode.vcf {output}.vcf && bgzip {output}.vcf"


rule effective_size:
    """
    Estimate population effective size
    """
    input:
        "{sample}/{sample}.chromosome.{chr}.vcf.gz"
    output:
        out = "{sample}/{sample}.chromosome.{chr}"
        effsize = "{sample}/effective_size.chromosome.{chr}"
    log:
        "{sample}/logs/effective_size.chromosome.{chr}.log"
    conda:
        "envs/vcftools.yaml"
    params:
        chr=config["chromosome"]
        windowsize=100000
        mu=config["mu"]
    shell:
        "vcftools --vcf {input} --window-pi {params.windowsize} --out {output.out} && eff_size=$(awk '{sum+=$5} END {sum=(sum/NR)/(4*{mu})); print sum}' {output.out}.windowed.pi && echo ${eff_size/.*} > {output.effsize}"


rule phasing_vcf:
    """
    Phase the vcf with Shapeit2
    Unphased haplotypes can overestimate recombination rate, hence false postive hotspots may be detected
    """
    input:
        vcf = "{sample}/{sample}.chromosome.{chr}.vcf.gz"
        effsize = "{sample}/effective_size.chromosome.{chr}"
    output:
        vcf = "{sample}/{sample}.chromosome.{chr}.phased.vcf"
        shapeit_log = "{sample}/logs/{sample}.chromosome.{chr}.shapeit.log"
    log:
        "{sample}/logs/phasing_vcf.chromosome.{chr}.log"
    conda:
        "envs/vcftools.yaml"
        "envs/shapeit.yaml"
    threads:
        workflow.cores * 0.9
    shell:
        "shapeit --input-vcf {input.vcf} --output-max {output.vcf} --effective-size $(cat {input.effsize}) --window 1 --thread {threads} --output-log {output.shapeit_log} --force"
        "shapeit -convert --input-haps {output} --output-vcf {output.vcf} --output-log {output.shapeit_log}"


rule haploidization:
    """
    Create pseudo-diploid vcf but with haploid genotypes
    """
    input:
        "{sample}/{sample}.chromosome.{chr}.phased.vcf"
    output:
        "{sample}/{sample}.chromosome.{chr}.haploid.vcf"
    log:
        "{sample}/logs/haploidization.chromosome.{chr}.log"
    conda:
        "envs/vcftools.yaml"
    shell:
        "sed 's/|.\t.|/|/g' {input} > {output}"
        "sed -r -i '/^#CHROM/s/\S+//g' {output}"
        "rm {input}"
        "bgzip {output}"


rule demography:
    """
    Pre-estimate demography with smc++ to use it as a prior in likelihood
    """
    input:
        vcf = "{sample}/{sample}.chromosome.{chr}.haploid.vcf.gz"
        poplist = "{sample}/poplist" # Provide a list of individuals to process in a 'poplist' file
    output:
        "{sample}/{sample}.chromosome.{chr}.smc"
    log:
        "{sample}/logs/demography.chromosome.{chr}.log"
    params:
        mu=config["mu"]
        knots=config["knots"]
    conda:
        "envs/vcftools.yaml"
        "envs/smcpp.yaml"
    threads:
        workflow.cores * 0.9
    shell:
        "tabix -p vcf {input.vcf} --csi"
        "individuals=$(cat {input.poplist} | tr '\n' ',')"
        "smc++ vcf2smc {input.vcf} {output} {chr} Pop1:$individuals"
        "smc++ estimate -o {sample}/smc/ {mu} {output} --knots {knots} --cores {threads}"
        "mv iterate.dat {sample}/smc/"
        "smc++ plot {sample}/plot_chromosome{chr}.pdf {sample}/smc/model.final.{chr}.json -c"
        "smc++ posterior {sample}/smc/model.final.{chr}.json {sample}/smc/posterior.chromosome.{chr}.smc {output} --heatmap {sample}/smc/heatmap.chromosome.{chr}.png --colorbar"


rule lookup_table:
    """
    Generate demography-aware look up tables
    """
    log:
        "{sample}/ldpop.chromosome{chr}.log"
    conda:
        "pyrho.yaml"
    shell:
        "pyrho make_table"

rule hyperparam:
    """
    Search the best hyperparameters for the dataset
    """


rule optimize:
    """
    Estimate fine-scale recombination rates with Pyrho
    """
    conda:
        "envs/pyrho.yaml"

rule compute_r2:
    """
    Compute the population-scaled recombination rate as infered by Pyrho
    """


rule ldhot:
    """
    Infer recombination hotspots
    """
