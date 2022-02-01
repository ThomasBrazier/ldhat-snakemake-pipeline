"""
Pipeline to estimate fine-scale recombination maps from polymorphism data
LDhelmet is the central element
"""

"""
Configuration of the analysis
i.e. dataset, name of chromosome, population to sample
"""
configfile: "config.yaml"


"""In addition to the configfile statement, config values can be overwritten via the command line"""
dataset=config["dataset"] # Name of your dataset directory and prefix of your vcf file
pop=config["pop"] # Specify the population once you have defined it with genetic clustering
chrom=config["chrom"] # Name of the chromosome to analyse in your 'sample' dataset

wdir='data/' + dataset

wildcard_constraints:
    wdir=wdir,
    dataset=dataset

rule all:
    """
    One ring to rule them all"
    """
    input:
        vcf = expand('{wdir}/{dataset}.trimmed.vcf.gz',wdir=wdir,dataset=dataset)
    shell:
        "echo 'Finished'"

# The global dataset is trimmed for SNPs
# and individuals if a 'samplelist' file is provided
# Otherwise all individuals are kept
rule trimming_vcf:
    """
    A first step to trim every 'sample' dataset to the same quality criteria
    """
    input:
        "{wdir}/{dataset}.vcf.gz"
    output:
    	"{wdir}/{dataset}.trimmed.vcf.gz"
    log:
        "{wdir}/logs/{dataset}.trimming_vcf.log"
    conda:
        "envs/vcftools.yaml"
    shell:
        """
        if [ -f "{wdir}/samplelist" ];
        then
            vcftools --gzvcf {input} --out {wdir}/out --recode --keep {wdir}/samplelist --maf config[maf] --max-missing config[maxmissing]
            mv {wdir}/out.recode.vcf {wdir}/{dataset}.trimmed.vcf
            bgzip -f {wdir}/{dataset}.trimmed.vcf
            rm {wdir}/out.log
        else
            vcftools --gzvcf {input} --out {wdir}/out --recode --maf config[maf] --maxmissing config[maxmissing]
            mv {wdir}/out.recode.vcf {wdir}/{dataset}.trimmed.vcf
            bgzip {wdir}/{dataset}.trimmed.vcf
            rm {wdir}/out.log
        fi
        """

rule population_struture:
    """
    Detects population structure and infer genetic clusters
    Produces diagnostic plots and summary statistics
    to help selecting the best population for recombination map
    """
    input:
        "{wdir}/{dataset}.trimmed.vcf.gz"
    output:
        dapc_plot = "{wdir}/popstructure/{dataset}.DAPC.pdf"
        summary_stats = "{wdir}/popstructure/{dataset}.statistics.csv}"
    log:
        "{wdir}/logs/{dataset}.population_structure.log"
    conda:
        "envs/Renv.yaml"
    shell:
        "Rscript popstructure.R {input}"

