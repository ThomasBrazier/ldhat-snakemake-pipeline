"""
Pipeline to estimate fine-scale recombination maps from polymorphism data
Test for convergence of LDhat MCMC chain

Run only one chromosome
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
    One ring to rule them all
    """
    input:
        target = expand("{wdirpop}/MCMC/{dataset}.{chrom}.bpen{bpen}.ldhat_MCMC.html", wdirpop=wdirpop, dataset=dataset, chrom=chrom, bpen=bpen)
    shell:
        "echo 'MCMC convergence assessment: finished'"



rule MCMC_report:
    """
    Produce a Rmarkdown/pdf report
    """
    input:
        "{wdirpop}/ldhat/{dataset}.{chrom}.bpen{bpen}.rates.txt.gz",
        "{wdirpop}/ldhat/{dataset}.{chrom}.bpen{bpen}.res.txt.gz",
        "{wdirpop}/ldhat/{dataset}.{chromosome}.ldhat.locs"
    output:
        "{wdirpop}/MCMC/{dataset}.{chrom}.bpen{bpen}.ldhat_MCMC.html"
    conda:
        "envs/Renv.yaml"
    shell:
        """
        Rscript ldhat_MCMC.R {wdirpop} {dataset} {chrom} {bpen}
        mv ldhat_MCMC.html {wdirpop}/MCMC/{dataset}.{chrom}.bpen{bpen}.ldhat_MCMC.html
        """
