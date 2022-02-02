"""
Pipeline to estimate fine-scale recombination maps from polymorphism data
Preprocesing data - Making a consistent population dataset
"""

"""
Configuration of the analysis
i.e. dataset, name of chromosome, population to sample
"""
configfile: "config.yaml"


"""In addition to the configfile statement, config values can be overwritten via the command line"""
dataset=config["dataset"] # Name of your dataset directory and prefix of your vcf file

wdir='data/' + dataset

wildcard_constraints:
    wdir=wdir,
    dataset=dataset

rule all:
    """
    One ring to rule them all"
    """
    input:
        expand("{wdir}/poplist", wdir=wdir)
    shell:
        "echo 'Preprocessing: Finished to infer population structure'"


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


# Use vcftools to export to .ped format with the --plink switch.  Then convert the .ped to .bed using Plink.  Faststructure will read the .bed format files (there are three files for each project)
rule vcf2structure:
    """
    Convert the vcf file to bed format for FastStructure
    """
    input:
        "{wdir}/{dataset}.trimmed.vcf.gz"
    output:
        bed = "{wdir}/{dataset}.bed",
        bam = "{wdir}/{dataset}.fam",
        bim = "{wdir}/{dataset}.bim"
    log:
        "{wdir}/logs/{dataset}.vcf2structure.log"
    conda:
        "envs/vcftools.yaml"
    shell:
        """
        plink --vcf {input} --out {wdir}/{dataset}
        """


rule faststructure:
    """
    FastStructure
    Detects population structure and infer genetic clusters
    Produces diagnostic plots and summary statistics
    to help selecting the best population for recombination map
    """
    input:
        bed = expand("{wdir}/{dataset}.bed", wdir=wdir, dataset=dataset)
    output:
        expand("{wdir}/structure/faststructure.{k}.log", wdir=wdir, k=range(1,7)),
        expand("{wdir}/structure/faststructure.{k}.meanP", wdir=wdir, k=range(1,7)),
        expand("{wdir}/structure/faststructure.{k}.varP", wdir=wdir, k=range(1,7)),
        expand("{wdir}/structure/faststructure.{k}.meanQ", wdir=wdir, k=range(1,7)),
        expand("{wdir}/structure/faststructure.{k}.varQ", wdir=wdir, k=range(1,7))
    log:
        expand("{wdir}/logs/faststructure.{k}.log", wdir=wdir, k=range(1,7))
    conda:
        "envs/faststructure.yaml"
    shell:
        """
        for k in {{1..7}}
        do
            python structure.py -K k --input={input.bed} --output={wdir}/structure/faststructure
        done
        """


rule poplist:
    """
    Identify the best K and classify individuals
    """
    input:
        expand("{wdir}/structure/faststructure.{k}.log", wdir=wdir, k=range(1,7)),
        expand("{wdir}/structure/faststructure.{k}.meanP", wdir=wdir, k=range(1,7)),
        expand("{wdir}/structure/faststructure.{k}.varP", wdir=wdir, k=range(1,7)),
        expand("{wdir}/structure/faststructure.{k}.meanQ", wdir=wdir, k=range(1,7)),
        expand("{wdir}/structure/faststructure.{k}.varQ", wdir=wdir, k=range(1,7))
    output:
        "{wdir}/poplist"
    shell:
        "echo 'poplist' > {wdir}/poplist"

