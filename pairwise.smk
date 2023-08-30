"""
Snakemake pipeline to estimate fine-scale recombination maps from polymorphism data
"""

"""
Configuration of the analysis
i.e. dataset, name of chromosome, population to sample
"""
include: "rules/config.smk"

rule all:
    """
    One ring to rule them all"
    """
    input:
        target = expand("{wdirpop}/{dataset}.{chrom}.bpen{bpen}.outfile.txt",wdirpop=wdirpop,dataset=dataset,chrom=chrom,bpen=bpen)
    shell:
        "echo 'Finished'"



include: "rules/sampling.smk"
include: "rules/phasing.smk"

rule lkgen:
    """
    Generate a complete look-up table from a preexisting one
    https://github.com/auton1/LDhat/tree/master/lk_files
    or generate a new one if explicitly asked in config.yaml
    """
    input:
        "{wdirpop}/{dataset}.chromosome.{chrom}.ldhat.vcf.gz"
    output:
        "{wdirpop}/ldhat/{dataset}.lookup.{chrom}.new_lk.txt"
    log:
        "{wdirpop}/logs/{dataset}.lookup.{chrom}.log"
    shell:
        """
	n=$(zcat {wdirpop}/{dataset}.chromosome.{chrom}.ldhat.vcf.gz | grep ^#CHROM | awk '{{print NF-9}}')
        n=$((2*$n))
        echo $n
	if [ "{config[completelk]}" == "no" ]
	then
	echo "Runnning lkgen"
	singularity exec --bind $PWD:/data ldhat.sif lkgen -prefix /data/{wdirpop}/ldhat/{dataset}.lookup.{chrom}. -lk /data/lk_files/lk_n100_t{config[theta]} -nseq $n
	else
	echo "Generate a new look-up table"
	singularity exec --bind $PWD:/data ldhat.sif complete -n $n -rhomax 100 -n_pts 101 -theta {config[theta]} -prefix /data/{wdirpop}/ldhat/{dataset}.lookup.{chrom}.
	fi
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
        lookup = "{wdirpop}/ldhat/{dataset}.lookup.{chrom}.new_lk.txt",
        vcf = "{wdirpop}/{dataset}.chromosome.{chrom}.ldhat.vcf.gz"
    output:
        "{wdirpop}/ldhat/{dataset}.{chrom}.{bpen}.ldhat.sites",
        "{wdirpop}/ldhat/{dataset}.{chrom}.{bpen}.ldhat.locs"
    log:
        "{wdirpop}/logs/{dataset}.ldhatconvert.{chrom}.bpen{bpen}.log"
    conda:
        "envs/vcftools.yaml"
    shell:
        """
        vcftools --gzvcf {input.vcf} --chr {chrom} --ldhat --out {wdirpop}/ldhat/{dataset}.{chrom}.{bpen}
        """


    rule pairwise:
        """
        Estimate pairwise recombination rates and theta per site with LDhat
        """
        input:
            "{wdirpop}/ldhat/{dataset}.{chrom}.{bpen}.ldhat.sites",
            "{wdirpop}/ldhat/{dataset}.{chrom}.{bpen}.ldhat.locs"
        output:
            "{wdirpop}/pairwise/{dataset}.{chrom}.bpen{bpen}.outfile.txt",
            "{wdirpop}/pairwise/{dataset}.{chrom}.bpen{bpen}.fits.txt",
            "{wdirpop}/pairwise/{dataset}.{chrom}.bpen{bpen}.window_out.txt",
            "{wdirpop}/pairwise/{dataset}.{chrom}.bpen{bpen}.output.txt"
        log:
            "{wdirpop}/logs/{dataset}.ldhatinterval.{chrom}.bpen{bpen}.log"
        shell:
            """
            iter={config[interval.iter]}
            samp={config[interval.samp]}
            bpen={config[interval.bpen]}
            singularity exec --bind $PWD:/data ldhat.sif pairwise -seq /data/{wdirpop}/ldhat/{dataset}.{chrom}.{bpen}.ldhat.sites -loc /data/{wdirpop}/ldhat/{dataset}.{chrom}.{bpen}.ldhat.locs -lk /data/{wdirpop}/ldhat/{dataset}.lookup.{chrom}.new_lk.txt -prefix /data/{wdirpop}/pairwise/{dataset}.{chrom}.bpen{bpen}. | tee {wdirpop}/pairwise/{dataset}.{chrom}.bpen{bpen}.output.txt
            """
