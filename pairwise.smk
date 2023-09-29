"""
Snakemake pipeline to estimate fine-scale recombination maps from polymorphism data
"""

"""
Configuration of the analysis
i.e. dataset, name of chromosome, population to sample
"""
include: "config.smk"

rule all:
    """
    One ring to rule them all"
    """
    input:
        target = expand("{wdirpop}/pairwise/{dataset}.{chrom}.rho.tsv",wdirpop=wdirpop,dataset=dataset,chrom=chrom)
    shell:
        "echo 'Finished'"



include: "sampling.smk"
include: "phasing.smk"

rule lkgen:
    """
    Generate a complete look-up table from a preexisting one
    https://github.com/auton1/LDhat/tree/master/lk_files
    or generate a new one if explicitly asked in config.yaml
    """
    input:
        "{wdirpop}/{dataset}.chromosome.{chrom}.ldhat.vcf.gz"
    output:
        "{wdirpop}/pairwise/{dataset}.lookup.{chrom}.new_lk.txt"
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
	singularity exec --bind $PWD:/data ldhat.sif lkgen -prefix /data/{wdirpop}/pairwise/{dataset}.lookup.{chrom}. -lk /data/lk_files/lk_n100_t{config[theta]} -nseq $n
	else
	echo "Generate a new look-up table"
	singularity exec --bind $PWD:/data ldhat.sif complete -n $n -rhomax {rhomax} -n_pts {n_pts} -theta {config[theta]} -prefix /data/{wdirpop}/pairwise/{dataset}.lookup.{chrom}.
	fi
	"""


rule split_dataset:
    """
    Split the dataset in pieces
    """
    input:
        "{wdirpop}/pairwise/{dataset}.lookup.{chrom}.new_lk.txt"
    output:
        "{wdirpop}/pairwise/{dataset}.{chrom}/nbatch"
    conda:
        "envs/vcftools.yaml"
    log:
        "{wdirpop}/logs/{dataset}.split_dataset.{chrom}.log"
    shell:
        """
 	# The first line splits up the snps into chunks of whatever size you want (-l) and then the next line loops over each file and subsets the vcf according.
        bcftools query -f'%CHROM\t%POS\n' {wdirpop}/{dataset}.chromosome.{chrom}.ldhat.vcf.gz > {wdirpop}/{dataset}.{chrom}.positions 
        python split_dataset.py {wdirpop}/{dataset}.{chrom}.positions {wdirpop}/pairwise/{dataset}.{chrom} {config[cut_size]} {config[cut_overlap]}
        nbatch=$(cat {wdirpop}/pairwise/{dataset}.{chrom}/nbatch_split)
        for i in $(seq $nbatch); do
        vcftools --gzvcf {wdirpop}/{dataset}.chromosome.{chrom}.ldhat.vcf.gz --chr {chrom} --positions {wdirpop}/pairwise/{dataset}.{chrom}/batch_$i.pos --recode --out {wdirpop}/pairwise/{dataset}.{chrom}/batch_$i
        done
        gzip -f {wdirpop}/pairwise/{dataset}.{chrom}/batch_*.recode.vcf
        echo $nbatch > {wdirpop}/pairwise/{dataset}.{chrom}/nbatch
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
        "{wdirpop}/pairwise/{dataset}.{chrom}/nbatch"
    output:
        "{wdirpop}/pairwise/convert.{dataset}.{chrom}.done"
    log:
        "{wdirpop}/logs/convert.{dataset}.{chrom}.log"
    conda:
        "envs/vcftools.yaml"
    shell:
        """
        nbatch=$(cat {wdirpop}/pairwise/{dataset}.{chrom}/nbatch)
        echo "nbatch = $nbatch"
        rm {wdirpop}/pairwise/convert.{dataset}.{chrom}.done || true
        for i in $(seq $nbatch); do
        vcftools --gzvcf {wdirpop}/pairwise/{dataset}.{chrom}/batch_$i.recode.vcf.gz --chr {chrom} --ldhat --out {wdirpop}/pairwise/{chrom}.batch_$i | tee -a {wdirpop}/pairwise/convert.{chrom}.log
        done
        echo "Done" > {wdirpop}/pairwise/convert.{dataset}.{chrom}.done
        """


rule pairwise:
    """
    Estimate pairwise recombination rates and theta per site with LDhat
    """
    input:
        "{wdirpop}/pairwise/convert.{dataset}.{chrom}.done"
    output:
        "{wdirpop}/pairwise/{dataset}.{chrom}.outfile.txt",
        "{wdirpop}/pairwise/{dataset}.{chrom}.fit.txt",
        "{wdirpop}/pairwise/{dataset}.{chrom}.output.txt"
    threads: workflow.cores
    log:
        "{wdirpop}/logs/pairwise.{dataset}.{chrom}.log"
    shell:
        """
        nbatch=$(cat {wdirpop}/pairwise/{dataset}.{chrom}/nbatch)
        echo "nbatch = $nbatch"
        for i in $(seq $nbatch)
        do
        echo -e "{config[pairwiseGrid]}\n{config[pairwiseSliding]}\n{config[pairwiseRmin]}\n{config[pairwiseMoment]}\n{config[pairwiseTest]}" | singularity exec --bind $PWD:/data ldhat.sif pairwise -seq /data/{wdirpop}/pairwise/{chrom}.batch_$i.ldhat.sites -loc /data/{wdirpop}/pairwise/{chrom}.batch_$i.ldhat.locs -lk /data/{wdirpop}/pairwise/{dataset}.lookup.{chrom}.new_lk.txt -prefix /data/{wdirpop}/pairwise/{dataset}.{chrom}.batch_$i. | tee {wdirpop}/pairwise/{dataset}.{chrom}.output.txt
        done
        cat {wdirpop}/pairwise/{dataset}.{chrom}.batch_*.outfile.txt > {wdirpop}/pairwise/{dataset}.{chrom}.outfile.txt
        cat {wdirpop}/pairwise/{dataset}.{chrom}.batch_*.fit.txt > {wdirpop}/pairwise/{dataset}.{chrom}.fit.txt
        """


rule summary:
    """
    Generate summary tables of pairwise results
    """
    input:
        "{wdirpop}/pairwise/{dataset}.{chrom}.outfile.txt",
        "{wdirpop}/pairwise/{dataset}.{chrom}.fit.txt",
        "{wdirpop}/pairwise/{dataset}.{chrom}.output.txt"
    output:
        "{wdirpop}/pairwise/{dataset}.{chrom}.rho.tsv"
    shell:
        """
        echo "dataset\tchromosome\tbatch\tstart\tend\ttheta\trho\tlk" > {wdirpop}/pairwise/{dataset}.{chrom}.rho.tsv
        nbatch=$(cat {wdirpop}/pairwise/{dataset}.{chrom}/nbatch)
        echo "nbatch = $nbatch"
        for i in $(seq $nbatch); do
        start=$(cat {wdirpop}/pairwise/{dataset}.{chrom}/batch_$i.pos | head -n 1 | awk '{{ print $2 }}')
        end=$(cat {wdirpop}/pairwise/{dataset}.{chrom}/batch_$i.pos | tail -n 1 | awk '{{ print $2 }}')
        theta=$(cat {wdirpop}/pairwise/{dataset}.{chrom}.batch_$i.outfile.txt | grep Theta | awk '{{ print $3 }}')
        rho=$(cat {wdirpop}/pairwise/{dataset}.{chrom}.batch_$i.outfile.txt | grep Maximum | awk '{{ print $5 }}')
        lk=$(cat {wdirpop}/pairwise/{dataset}.{chrom}.batch_$i.outfile.txt | grep Maximum | awk '{{ print $9 }}')
        echo "{dataset}\t{chrom}\t$i\t$start\t$end\t$theta\t$rho\t$lk" >> {wdirpop}/pairwise/{dataset}.{chrom}.rho.tsv
        done
        """



## TODO Lk matrix
## TODO Summary report
