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
