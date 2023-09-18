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
        target = expand("{wdirpop}/pairwise/{dataset}.{chrom}.outfile.txt",wdirpop=wdirpop,dataset=dataset,chrom=chrom)
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


#rule convert:
#    """
#    Produce input files for LDhat
#    Generate sites.txt and locs.txt with vcftools
#    --ldhat
#    --ldhat-geno
#    These options output data in LDhat format. This option requires the "--chr" filter option to also be used. The first option outputs phased data only, and therefore also implies "--phased" be used, leading to unphased individuals and genotypes being excluded. The second option treats all of the data as unphased, and therefore outputs LDhat files in genotype/unphased format. Two output files are generated with the suffixes ".ldhat.sites" and ".ldhat.locs", which correspond to the LDhat "sites" and "locs" input files respectively.
#    """
#    input:
#        lookup = "{wdirpop}/ldhat/{dataset}.lookup.{chrom}.new_lk.txt",
#        vcf = "{wdirpop}/{dataset}.chromosome.{chrom}.ldhat.vcf.gz"
#    output:
#        "{wdirpop}/ldhat/{dataset}.{chrom}.{bpen}.ldhat.sites",
#        "{wdirpop}/ldhat/{dataset}.{chrom}.{bpen}.ldhat.locs"
#    log:
#        "{wdirpop}/logs/{dataset}.ldhatconvert.{chrom}.bpen{bpen}.log"
#    conda:
#        "envs/vcftools.yaml"
#    shell:
#        """
#        vcftools --gzvcf {input.vcf} --chr {chrom} --ldhat --out {wdirpop}/ldhat/{dataset}.{chrom}.{bpen}
#        """


#    rule pairwise:
#        """
#        Estimate pairwise recombination rates and theta per site with LDhat
#        """
#        input:
#            "{wdirpop}/ldhat/{dataset}.{chrom}.{bpen}.ldhat.sites",
#            "{wdirpop}/ldhat/{dataset}.{chrom}.{bpen}.ldhat.locs"
#        output:
#            "{wdirpop}/pairwise/{dataset}.{chrom}.bpen{bpen}.outfile.txt",
#            "{wdirpop}/pairwise/{dataset}.{chrom}.bpen{bpen}.fits.txt",
#            "{wdirpop}/pairwise/{dataset}.{chrom}.bpen{bpen}.window_out.txt",
#            "{wdirpop}/pairwise/{dataset}.{chrom}.bpen{bpen}.output.txt"
#        log:
#            "{wdirpop}/logs/{dataset}.ldhatinterval.{chrom}.bpen{bpen}.log"
#        shell:
#            """
#            iter={config[interval.iter]}
#            samp={config[interval.samp]}
#            bpen={config[interval.bpen]}
#            echo -e "0\n0\n0\n0\n0" | singularity exec --bind $PWD:/data ldhat.sif pairwise -seq /data/{wdirpop}/ldhat/{dataset}.{chrom}.{bpen}.ldhat.sites -loc /data/{wdirpop}/ldhat/{dataset}.{chrom}.{bpen}.ldhat.locs -lk /data/{wdirpop}/ldhat/{dataset}.lookup.{chrom}.new_lk.txt -prefix /data/{wdirpop}/pairwise/{dataset}.{chrom}.bpen{bpen}. | tee {wdirpop}/pairwise/{dataset}.{chrom}.bpen{bpen}.output.txt
#            """


rule split_dataset:
    """
    Split the dataset in pieces
    """
    input:
        "{wdirpop}/ldhat/{dataset}.lookup.{chrom}.new_lk.txt"
    output:
        "{wdirpop}/ldhat/{dataset}.{chrom}/nbatch"
    conda:
        "envs/vcftools.yaml"
    log:
        "{wdirpop}/logs/{dataset}.split_dataset.{chrom}.log"
    shell:
        """
 	    # The first line splits up the snps into chunks of whatever size you want (-l) and then the next line loops over each file and subsets the vcf according.
        bcftools query -f'%CHROM\t%POS\n' {wdirpop}/{dataset}.chromosome.{chrom}.ldhat.vcf.gz > {wdirpop}/{dataset}.{chrom}.positions 
        python split_dataset.py {wdirpop}/{dataset}.{chrom}.positions {wdirpop}/ldhat/{dataset}.{chrom} {config[cut_size]} {config[cut_overlap]}
        nbatch=$(cat {wdirpop}/ldhat/{dataset}.{chrom}/nbatch_split)
        for i in $(seq $nbatch); do
        vcftools --gzvcf {wdirpop}/{dataset}.chromosome.{chrom}.ldhat.vcf.gz --chr {chrom} --positions {wdirpop}/ldhat/{dataset}.{chrom}/batch_$i.pos --recode --out {wdirpop}/ldhat/{dataset}.{chrom}/batch_$i
        done
        gzip -f {wdirpop}/ldhat/{dataset}.{chrom}/batch_*.recode.vcf
        echo $nbatch > {wdirpop}/ldhat/{dataset}.{chrom}/nbatch
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
        "{wdirpop}/ldhat/{dataset}.{chrom}/nbatch"
    output:
        "{wdirpop}/ldhat/{dataset}.{chrom}/convert.done"
    log:
        "{wdirpop}/logs/{dataset}.ldhatconvert.{chrom}.log"
    conda:
        "envs/vcftools.yaml"
    shell:
        """
        nbatch=$(cat {wdirpop}/ldhat/{dataset}.{chrom}/nbatch)
        echo "nbatch = $nbatch"
        rm  {wdirpop}/ldhat/{dataset}.{chrom}/convert.done || true
        touch {wdirpop}/ldhat/{dataset}.{chrom}/convert.done
        for i in $(seq $nbatch); do
        vcftools --gzvcf {wdirpop}/ldhat/{dataset}.{chrom}/batch_$i.recode.vcf.gz --chr {chrom} --ldhat --out {wdirpop}/ldhat/{dataset}.{chrom}/batch_$i | tee -a {wdirpop}/ldhat/{dataset}.{chrom}/convert.done
        done
 	    echo "Done" > {wdirpop}/ldhat/{dataset}.{chrom}/convert.done
        """


rule pairwise:
    """
    Estimate pairwise recombination rates and theta per site with LDhat
    """
    input:
        "{wdirpop}/ldhat/{dataset}.{chrom}.{bpen}.ldhat.sites",
        "{wdirpop}/ldhat/{dataset}.{chrom}.{bpen}.ldhat.locs",
        "{wdirpop}/ldhat/{dataset}.{chrom}/convert.done"
    output:
        "{wdirpop}/pairwise/{dataset}.{chrom}.outfile.txt",
        "{wdirpop}/pairwise/{dataset}.{chrom}.fits.txt",
        "{wdirpop}/pairwise/{dataset}.{chrom}.window_out.txt",
        "{wdirpop}/pairwise/{dataset}.{chrom}.output.txt"
    log:
        "{wdirpop}/logs/{dataset}.pairwise.{chrom}.log"
    shell:
        """
        #echo -e "0\n0\n0\n0\n0" | singularity exec --bind $PWD:/data ldhat.sif pairwise -seq /data/{wdirpop}/ldhat/{dataset}.{chrom}.{bpen}.ldhat.sites -loc /data/{wdirpop}/ldhat/{dataset}.{chrom}.{bpen}.ldhat.locs -lk /data/{wdirpop}/ldhat/{dataset}.lookup.{chrom}.new_lk.txt -prefix /data/{wdirpop}/pairwise/{dataset}.{chrom}.bpen{bpen}. | tee {wdirpop}/pairwise/{dataset}.{chrom}.bpen{bpen}.output.txt
        nbatch=$(cat {wdirpop}/ldhat/{dataset}.{chrom}/nbatch)
        echo "nbatch = $nbatch"
        for i in $(seq $nbatch)
        do
        echo -e "0\n0\n0\n0\n0" | singularity exec --bind $PWD:/data ldhat.sif pairwise -seq /data/{wdirpop}/ldhat/{dataset}.{chrom}/batch_$i.ldhat.sites -loc /data/{wdirpop}/ldhat/{dataset}.{chrom}/batch_$i.ldhat.locs -lk /data/{wdirpop}/ldhat/{dataset}.lookup.{chrom}.new_lk.txt -prefix /data/{wdirpop}/pairwise/{dataset}.{chrom}/{dataset}.{chrom}.batch_$i. | tee {wdirpop}/pairwise/{dataset}.{chrom}.output.txt
        done
        cat {wdirpop}/pairwise/{dataset}.{chrom}/*.outfile.txt > {wdirpop}/pairwise/{dataset}.{chrom}.outfile.txt
        cat {wdirpop}/pairwise/{dataset}.{chrom}/*.fits.txt > {wdirpop}/pairwise/{dataset}.{chrom}.fits.txt
        cat {wdirpop}/pairwise/{dataset}.{chrom}/*.window_out.txt > {wdirpop}/pairwise/{dataset}.{chrom}.window_out.txt
        cat {wdirpop}/pairwise/{dataset}.{chrom}/*.output.txt > {wdirpop}/pairwise/{dataset}.{chrom}.output.txt
        """

