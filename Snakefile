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

wdir='data/' + dataset

wildcard_constraints:
    wdir=wdir,
    dataset=dataset

rule all:
    """
    One ring to rule them all"
    """
    input:
        vcf = expand('{wdir}/{dataset}.pop.vcf.gz',wdir=wdir,dataset=dataset)
    shell:
        "echo 'Finished'"

# The global dataset is trimmed for SNPs
# and individuals if a 'samplelist' file is provided
# Otherwise all individuals are kept
rule sampling_pop:
    """
    A first step to trim every 'sample' dataset to the same quality criteria
    """
    input:
        vcf = "{wdir}/{dataset}.vcf.gz"
        poplist = "{wdir}/poplist"
    output:
    	"{wdir}/{dataset}.pop.vcf.gz"
    log:
        "{wdir}/logs/{dataset}.sampling_pop.log"
    conda:
        "envs/vcftools.yaml"
    shell:
        """
        if [ -f "{input.poplist}" ];
        then
            vcftools --gzvcf {input.vcf} --out {wdir}/out --recode --keep {input.poplist} --maf config[maf] --max-missing config[maxmissing]
            mv {wdir}/out.recode.vcf {wdir}/{dataset}.pop.vcf
            bgzip -f {wdir}/{dataset}.pop.vcf
            rm {wdir}/out.log
        else
            echo "<poplist> is a mandatory input file"
        fi
        """
