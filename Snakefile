"""
Snakemake pipeline to estimate fine-scale recombination maps from polymorphism data
"""

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


rule all:
    """
    One ring to rule them all"
    """
    input:
        target = expand("{wdirpop}/ldhot/{dataset}.{chrom}.bpen{bpen}.hotspots.txt.gz",wdirpop=wdirpop,dataset=dataset,chrom=chrom,bpen=bpen),
    shell:
        "echo 'Finished'"


rule poplist:
    """
    Sample individuals from a given genetic cluster (number of K to retain and name of the cluster) specified in config.yaml
    """
    output:
        "{wdirpop}/poplist",
        temporary("{wdirpop}/structure/cluster.tmp"),
        temporary("{wdirpop}/structure/cluster.tmp2"),
        temporary("{wdirpop}/structure/cluster")
    log:
        "{wdirpop}/logs/poplist.log"
    conda:
        "envs/vcftools.yaml"
    shell:
        """
        perl -ane '$r = 0; for my $i (1 .. $#F) {{$r = $i if $F[$i] > $F[$r];}} print $r + 1, " ";' < {wdir}/structure/faststructure.{K}.meanQ > {wdirpop}/structure/cluster.tmp
        sed -i -e "s/\\s\\+/\\n/g" {wdirpop}/structure/cluster.tmp
        paste {wdir}/indlist {wdirpop}/structure/cluster.tmp > {wdirpop}/structure/cluster.tmp2
        paste {wdirpop}/structure/cluster.tmp2 {wdir}/structure/faststructure.{K}.meanQ > {wdirpop}/structure/cluster
        awk -F' ' '{{if($2=={pop}) print $1}}' {wdirpop}/structure/cluster > {wdirpop}/poplist
        """


rule sampling_pop:
    """
    A first step to trim every 'sample' dataset to the same quality criteria
    Filtering:
    - keep only biallelic
    - perform HWE test on site to remove excess of heterozygotes
    - maf
    - missing data
    - Quality score
    """
    input:
        "{wdirpop}/poplist"
    output:
    	"{wdirpop}/{dataset}.pop.vcf.gz"
    log:
        "{wdirpop}/logs/{dataset}.sampling_pop.log"
    conda:
        "envs/vcftools.yaml"
    shell:
        """
        # Filter chromosomes and keep only bi-allelic alelles
        if [ {config[minQ]} -eq 0 ]
        then
        vcftools --gzvcf {wdir}/{dataset}.vcf.gz --out {wdirpop}/out --recode --keep {wdirpop}/poplist --maf {config[maf]} --max-missing {config[maxmissing]} --min-alleles 2 --max-alleles 2
	else
        vcftools --gzvcf {wdir}/{dataset}.vcf.gz --out {wdirpop}/out --recode --keep {wdirpop}/poplist --maf {config[maf]} --max-missing {config[maxmissing]} --min-alleles 2 --max-alleles 2 --minQ {config[minQ]}
	fi
        mv {wdirpop}/out.recode.vcf {wdirpop}/{dataset}.pop.vcf
	bgzip -f {wdirpop}/{dataset}.pop.vcf
        bcftools norm -d all {wdirpop}/{dataset}.pop.vcf.gz -o {wdirpop}/{dataset}.pop.vcf
	bgzip -f {wdirpop}/{dataset}.pop.vcf
        tabix --csi {wdirpop}/{dataset}.pop.vcf.gz
        """


rule effective_size:
    """
    Estimate population effective size
    """
    input:
        "{wdirpop}/{dataset}.pop.vcf.gz"
    output:
        "{wdirpop}/statistics/{dataset}.effective_size"
    log:
        "{wdirpop}/logs/{dataset}.effective_size.log"
    conda:
        "envs/vcftools.yaml"
    params:
        windowsize=100000
    shell:
        """
        vcftools --gzvcf {input} --window-pi {params.windowsize} --out {wdirpop}/statistics/vcftools
        eff_size=$(awk '{{sum+=$5}} END {{sum=(sum/NR)/(4*{config[mu]}); print sum}}' {wdirpop}/statistics/vcftools.windowed.pi)
        echo ${{eff_size/.*}} > {output}
        """


# From now, analyses are performed on individual chromosomes
rule split_chromosome:
    """
    Split the entire vcf into one vcf per chromosome
    Deal with missing data: remove SNPs with missing data for more than 90% of individuals
    """
    input:
        "{wdirpop}/statistics/{dataset}.effective_size"
    output:
        "{wdirpop}/{dataset}.chromosome.{chrom}.vcf.gz"
    log:
        "{wdirpop}/logs/{dataset}.split_chromosome.{chrom}.log"
    conda:
        "envs/vcftools.yaml"
    shell:
        """
        vcftools --gzvcf {wdirpop}/{dataset}.pop.vcf.gz --out {wdirpop}/out --recode --chr {chrom} --max-missing {config[maxmissing]} --min-alleles 2 --max-alleles 2
        mv {wdirpop}/out.recode.vcf {wdirpop}/{dataset}.chromosome.{chrom}.vcf
        bgzip -f {wdirpop}/{dataset}.chromosome.{chrom}.vcf
        #tabix -f -p vcf {wdirpop}/{dataset}.chromosome.{chrom}.vcf.gz
        tabix -f --csi -p vcf {wdirpop}/{dataset}.chromosome.{chrom}.vcf.gz
        """


rule phasing_vcf:
    """
    Phase the vcf with Shapeit2
    Unphased haplotypes can overestimate recombination rate, hence false postive hotspots may be detected
    https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html#output
    """
    input:
        "{wdirpop}/{dataset}.chromosome.{chrom}.vcf.gz"
    output:
        "{wdirpop}/{dataset}.chromosome.{chrom}.phased.vcf.gz"
    log:
        "{wdirpop}/logs/{dataset}.chromosome.{chrom}.phasing_vcf.log"
    conda:
        "envs/shapeit.yaml"
    shell:
        """
	    # Remove --thread {config[cores]} if causing errors
        shapeit --input-vcf {wdirpop}/{dataset}.chromosome.{chrom}.vcf.gz --output-max {wdirpop}/{dataset}.phased.chromosome.{chrom} --effective-size $(cat {wdirpop}/statistics/{dataset}.effective_size) --window {config[shapeitWindow]} --thread {config[cores]} --output-log {wdirpop}/logs/{dataset}.chromosome.{chrom}.shapeit.log --force
        shapeit -convert --input-haps {wdirpop}/{dataset}.phased.chromosome.{chrom} --output-vcf {wdirpop}/{dataset}.chromosome.{chrom}.phased.vcf --output-log {wdirpop}/logs/{dataset}.chromosome.{chrom}.shapeit.convert.log
        # replace header in vcf to keep information of contig length
        zcat {wdir}/{dataset}.vcf.gz | grep '^#' > {wdirpop}/newheader
        cat {wdirpop}/newheader | grep -v '^#CHROM' > {wdirpop}/newheader2
        cat {wdirpop}/{dataset}.chromosome.{chrom}.phased.vcf | grep '^#CHROM' > {wdirpop}/colnames
        cat {wdirpop}/{dataset}.chromosome.{chrom}.phased.vcf | grep -v '^#' > {wdirpop}/newvcf
        cat {wdirpop}/newheader2 {wdirpop}/colnames {wdirpop}/newvcf > {wdirpop}/{dataset}.chromosome.{chrom}.phased.vcf
        bgzip -f {wdirpop}/{dataset}.chromosome.{chrom}.phased.vcf
        rm {wdirpop}/newheader {wdirpop}/newheader2 {wdirpop}/colnames {wdirpop}/newvcf
        """


rule pseudodiploid:
    """
    Make pseudodiploids (phased haplotypes) to take into account homozygotes in high-selfing rates species
    Method 1. Keep only one haplotype per individual and reconstruct diploids
    Method 2. Keep both haplotypes and resample among individuals
    """
    input:
        "{wdirpop}/{dataset}.chromosome.{chrom}.phased.vcf.gz"
    output:
        "{wdirpop}/{dataset}.chromosome.{chrom}.pseudodiploid.vcf.gz"
    log:
        "{wdirpop}/logs/{dataset}.chromosome.{chrom}.pseudodiploid.log"
    conda:
        "envs/Renv.yaml"
    shell:
        """
        if [ {config[pseudodiploid]} -eq 1 ]; then
            Rscript pseudodiploids.R {wdirpop} {chrom} 1
        else
	        if [ {config[pseudodiploid]} -eq 2 ]; then
	            Rscript pseudodiploids.R {wdirpop} {chrom} 2
            else
	            cp {input} {output}
	        fi
        fi
        """


rule gzpseudodiploid:
    """
    Make pseudodiploids (phased haplotypes) to take into account homozygotes in high-selfing rates species
    """
    input:
        "{wdirpop}/{dataset}.chromosome.{chrom}.pseudodiploid.vcf.gz"
    output:
        "{wdirpop}/{dataset}.chromosome.{chrom}.pseudodiploid.vcf.gz.csi"
    log:
        "{wdirpop}/logs/{dataset}.chromosome.{chrom}.gzpseudodiploid.log"
    conda:
        "envs/vcftools.yaml"
    shell:
        """
        gunzip {wdirpop}/{dataset}.chromosome.{chrom}.pseudodiploid.vcf.gz
        bgzip {wdirpop}/{dataset}.chromosome.{chrom}.pseudodiploid.vcf
        tabix -p vcf {wdirpop}/{dataset}.chromosome.{chrom}.pseudodiploid.vcf.gz --csi
        """


# TODO A step to reduce sample size to a random subset
# Memory issues when sample size is high
rule subset_ldhat:
    """
    Subset a random sample of individuals
    """
    input:
        "{wdirpop}/{dataset}.chromosome.{chrom}.pseudodiploid.vcf.gz"
    output:
        "{wdirpop}/{dataset}.chromosome.{chrom}.ldhat.vcf.gz"
    log:
        "{wdirpop}/logs/{dataset}.make_subset_ldhat.{chrom}.log"
    conda:
        "envs/vcftools.yaml"
    shell:
        """
	RANDOM=42
	vcftools --gzvcf {wdirpop}/{dataset}.chromosome.{chrom}.pseudodiploid.vcf.gz --out {wdirpop}/out --recode --max-indv {config[subset]} --maf {config[maf]} --max-missing {config[maxmissing]}
        mv {wdirpop}/out.recode.vcf {wdirpop}/{dataset}.chromosome.{chrom}.ldhat.vcf
	bgzip -f {wdirpop}/{dataset}.chromosome.{chrom}.ldhat.vcf
        """

#rule complete:
#    """
#    Generate a complete look-up table
#    """
#    input:
#        "{wdirpop}/{dataset}.chromosome.{chrom}.ldhat.vcf.gz"
#    output:
#        "{wdirpop}/ldhat/{dataset}.lookup.{chrom}"
#    log:
#        "{wdirpop}/logs/{dataset}.lookup.{chrom}.log"
#    shell:
#        """
#        n=$(zcat {wdirpop}/{dataset}.chromosome.{chrom}.ldhat.vcf.gz | grep ^#CHROM | awk '{{print NF-9}}')
#        n=$((2*$n))
#	 echo $n
#        singularity exec --bind $PWD:/data ldhat.sif complete -n $n -rhomax 100 -n_pts 101 -theta {config[theta} -prefix {wdirpop}/ldhat/{dataset}.lookup.{chrom}
#	 """


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


if config["large_sample"] == "yes":
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
            for i in $(seq $nbatch); do
            vcftools --gzvcf {wdirpop}/ldhat/{dataset}.{chrom}/batch_$i.recode.vcf.gz --chr {chrom} --ldhat --out {wdirpop}/ldhat/{dataset}.{chrom}/batch_$i
            done
	    echo "Done" > {wdirpop}/ldhat/{dataset}.{chrom}/convert.done
            """


    rule interval_split:
        """
        Estimate a recombination landscape with LDhat interval
        """
        input:
            "{wdirpop}/ldhat/{dataset}.{chrom}/convert.done"
        output:
            "{wdirpop}/ldhat/{dataset}.{chrom}/interval_bpen{bpen}.done"
        log:
            "{wdirpop}/logs/{dataset}.ldhatinterval.{chrom}.bpen{bpen}.log"
        conda:
            "envs/vcftools.yaml"
        shell:
            """
            iter={config[interval.iter]}
            samp={config[interval.samp]}
            bpen={config[interval.bpen]}
            nbatch=$(cat {wdirpop}/ldhat/{dataset}.{chrom}/nbatch)
            echo "nbatch = $nbatch"
            for i in $(seq $nbatch); do
            singularity exec --bind $PWD:/data ldhat.sif interval -seq /data/{wdirpop}/ldhat/{dataset}.{chrom}/batch_$i.ldhat.sites -loc /data/{wdirpop}/ldhat/{dataset}.{chrom}/batch_$i.ldhat.locs -lk /data/{wdirpop}/ldhat/{dataset}.lookup.{chrom}.new_lk.txt -its $iter -bpen $bpen -samp $samp -prefix /data/{wdirpop}/ldhat/{dataset}.{chrom}/bpen{bpen}.batch_$i.
            done
            echo "Done" > {wdirpop}/ldhat/{dataset}.{chrom}/interval_bpen{bpen}.done
            """


    rule stat_split:
        """
        Compute statistics on interval
        """
        input:
            "{wdirpop}/ldhat/{dataset}.{chrom}/interval_bpen{bpen}.done"
        output:
            "{wdirpop}/ldhat/{dataset}.{chrom}/stat_bpen{bpen}.done"
        log:
            "{wdirpop}/logs/{dataset}.ldhatstat.{chrom}.bpen{bpen}.log"
        conda:
            "envs/vcftools.yaml"
        shell:
            """
            burn={config[ldhat.burn]}
            nbatch=$(cat {wdirpop}/ldhat/{dataset}.{chrom}/nbatch)
            echo "nbatch = $nbatch"
            for i in $(seq $nbatch); do
            singularity exec --bind $PWD:/data ldhat.sif stat -input /data/{wdirpop}/ldhat/{dataset}.{chrom}/bpen{bpen}.batch_$i.rates.txt -burn $burn -loc /data/{wdirpop}/ldhat/{dataset}.{chrom}/batch_$i.ldhat.locs -prefix /data/{wdirpop}/ldhat/{dataset}.{chrom}/bpen{bpen}.batch_$i.
	    done
	    echo "Done" > {wdirpop}/ldhat/{dataset}.{chrom}/stat_bpen{bpen}.done
            """


    rule concatenate:
        """
        Concatenate .res.txt files
        """
        input:
            "{wdirpop}/ldhat/{dataset}.{chrom}/stat_bpen{bpen}.done"
        output:
            "{wdirpop}/ldhat/{dataset}.{chrom}.bpen{bpen}.res.txt"
        log:
            "{wdirpop}/logs/{dataset}.concatenate.{chrom}.bpen{bpen}.log"
        shell:
            """
            nbatch=$(cat {wdirpop}/ldhat/{dataset}.{chrom}/nbatch)
            echo "Concatenate into one file"
            overlap={config[cut_overlap]}
            chunk={config[cut_size]}
            bigchunk=$(echo $(( $chunk-$overlap/2 )))
            smalloverlap=$(echo $(( $overlap/2 )))
            echo $bigchunk
            echo $smalloverlap
            cd {wdirpop}/ldhat/{dataset}.{chrom}/
            echo $PWD
            echo "First chunk"
            echo $nbatch
            cat bpen{bpen}.batch_1.res.txt | grep -v "\-1\.00" | grep -v "Loci" | head -n $bigchunk > bpen{bpen}.res_noheader.txt
            for i in $(seq 2 $(( $nbatch-1 )))
            do
            echo $i
            cat bpen{bpen}.batch_$i.res.txt | grep -v "\-1" | grep -v "Loci" | head -n $bigchunk | tail -n +$(( $smalloverlap+1 )) >> bpen{bpen}.res_noheader.txt
            done
            echo "End of loop on split files."
            cat bpen{bpen}.batch_$nbatch.res.txt | grep -v "\-1" | grep -v "Loci" | tail -n +$(( $smalloverlap+1 )) >> bpen{bpen}.res_noheader.txt
            cd ../../../../..
            echo "Loci	Mean_rho	Median	L95	U95" > {wdirpop}/ldhat/{dataset}.{chrom}/bpen{bpen}.header
            Loci="-1.000"
            MeanRho=$(cat {wdirpop}/ldhat/{dataset}.{chrom}/bpen{bpen}.res_noheader.txt | awk '{{s+=$2}} END {{printf "%.0f", s}}')
            Median=$(cat {wdirpop}/ldhat/{dataset}.{chrom}/bpen{bpen}.res_noheader.txt | awk '{{s+=$3}} END {{printf "%.0f", s}}')
            L95=$(cat {wdirpop}/ldhat/{dataset}.{chrom}/bpen{bpen}.res_noheader.txt | awk '{{s+=$4}} END {{printf "%.0f", s}}')
            U95=$(cat {wdirpop}/ldhat/{dataset}.{chrom}/bpen{bpen}.res_noheader.txt | awk '{{s+=$5}} END {{printf "%.0f", s}}')
            echo "$Loci     $MeanRho        $Median $L95    $U95"
            echo "$Loci	$MeanRho	$Median	$L95	$U95" >> {wdirpop}/ldhat/{dataset}.{chrom}/bpen{bpen}.header
            cat {wdirpop}/ldhat/{dataset}.{chrom}/bpen{bpen}.header > {wdirpop}/ldhat/{dataset}.{chrom}/bpen{bpen}.res.txt
            cat {wdirpop}/ldhat/{dataset}.{chrom}/bpen{bpen}.res_noheader.txt >> {wdirpop}/ldhat/{dataset}.{chrom}/bpen{bpen}.res.txt
            cp {wdirpop}/ldhat/{dataset}.{chrom}/bpen{bpen}.res.txt {wdirpop}/ldhat/{dataset}.{chrom}.bpen{bpen}.res.txt
            bash concat_ldhat_rates.sh {wdirpop} {dataset} {chrom} {bpen}
            """


    rule LDhot_convert:
        """
        Reconcile ldhat.locs/ldhat.sites and res.txt for LDhot input
        """
        input:
            expand("{wdirpop}/ldhat/{dataset}.{chrom}.bpen{bpen}.res.txt", wdirpop=wdirpop, dataset=dataset, chrom=chrom, bpen=bpen)
        output:
            "{wdirpop}/ldhat/{dataset}.{chrom}.ldhat.sites",
            "{wdirpop}/ldhat/{dataset}.{chrom}.ldhat.locs"
        log:
            "{wdirpop}/logs/{dataset}.ldhatconvert.{chrom}.log"
        conda:
            "envs/vcftools.yaml"
        shell:
            """
            vcftools --gzvcf {wdirpop}/{dataset}.chromosome.{chrom}.ldhat.vcf.gz --chr {chrom} --ldhat --out {wdirpop}/ldhat/{dataset}.{chrom}
            """
elif config["large_sample"] == "no":
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
            "{wdirpop}/ldhat/{dataset}.{chrom}.ldhat.sites",
            "{wdirpop}/ldhat/{dataset}.{chrom}.ldhat.locs"
        log:
            "{wdirpop}/logs/{dataset}.ldhatconvert.{chrom}.bpen{bpen}.log"
        conda:
            "envs/vcftools.yaml"
        shell:
            """
            vcftools --gzvcf {input.vcf} --chr {chrom} --ldhat --out {wdirpop}/ldhat/{dataset}.{chrom}
            """

    rule interval:
        """
        Estimate a recombination landscape with LDhat
        """
        input:
            "{wdirpop}/ldhat/{dataset}.{chrom}.ldhat.sites",
            "{wdirpop}/ldhat/{dataset}.{chrom}.ldhat.locs"
        output:
            "{wdirpop}/ldhat/{dataset}.{chrom}.bpen{bpen}.new_lk.txt",
            "{wdirpop}/ldhat/{dataset}.{chrom}.bpen{bpen}.bounds.txt",
            "{wdirpop}/ldhat/{dataset}.{chrom}.bpen{bpen}.rates.txt",
            temporary("{wdirpop}/ldhat/{dataset}.{chrom}.bpen{bpen}.type_table.txt")
        log:
            "{wdirpop}/logs/{dataset}.ldhatinterval.{chrom}.bpen{bpen}.log"
        shell:
            """
            iter={config[interval.iter]}
            samp={config[interval.samp]}
            bpen={config[interval.bpen]}
            singularity exec --bind $PWD:/data ldhat.sif interval -seq /data/{wdirpop}/ldhat/{dataset}.{chrom}.ldhat.sites -loc /data/{wdirpop}/ldhat/{dataset}.{chrom}.ldhat.locs -lk /data/{wdirpop}/ldhat/{dataset}.lookup.{chrom}.new_lk.txt -its $iter -bpen $bpen -samp $samp -prefix /data/{wdirpop}/ldhat/{dataset}.{chrom}.bpen{bpen}.
            """

    rule stat:
        """
        Compute statistics on interval
        """
        input:
            "{wdirpop}/ldhat/{dataset}.{chrom}.bpen{bpen}.new_lk.txt",
            "{wdirpop}/ldhat/{dataset}.{chrom}.bpen{bpen}.bounds.txt",
            "{wdirpop}/ldhat/{dataset}.{chrom}.bpen{bpen}.rates.txt"
        output:
            "{wdirpop}/ldhat/{dataset}.{chrom}.bpen{bpen}.res.txt",
	    "{wdirpop}/ldhat/{dataset}.{chrom}.bpen{bpen}.bounds.txt.gz",
            "{wdirpop}/ldhat/{dataset}.{chrom}.bpen{bpen}.rates.txt.gz"
        log:
            "{wdirpop}/logs/{dataset}.ldhatstat.{chrom}.bpen{bpen}.log"
        shell:
            """
            burn={config[ldhat.burn]}
            singularity exec --bind $PWD:/data ldhat.sif stat -input /data/{wdirpop}/ldhat/{dataset}.{chrom}.bpen{config[bpen]}.rates.txt -burn $burn -loc /data/{wdirpop}/ldhat/{dataset}.{chrom}.ldhat.locs -prefix /data/{wdirpop}/ldhat/{dataset}.{chrom}.bpen{bpen}.
            # Compress intermediary files
	    gzip -f {wdirpop}/ldhat/{dataset}.{chrom}.bpen{bpen}.rates.txt
            gzip -f {wdirpop}/ldhat/{dataset}.{chrom}.bpen{bpen}.bounds.txt
	    """


rule LDhot:
    """
    Infer recombination hotspots
    LDhot
    """
    input:
        "{wdirpop}/ldhat/{dataset}.{chrom}.ldhat.sites",
        "{wdirpop}/ldhat/{dataset}.{chrom}.ldhat.locs"
    output:
        "{wdirpop}/ldhot/{dataset}.{chrom}.bpen{bpen}.hotspots.txt.gz",
	"{wdirpop}/ldhat/{dataset}.{chrom}.bpen{bpen}.res.txt.gz",
        temporary("{wdirpop}/ldhot/{dataset}.{chrom}.bpen{bpen}.log")
    log:
        "{wdirpop}/logs/{dataset}.{chrom}.bpen{bpen}.ldhot.log"
    shell:
        """
        nsim={config[ldhot.nsim]}
        singularity exec --bind $PWD:/data ldhot.sif ldhot --seq /data/{wdirpop}/ldhat/{dataset}.{chrom}.ldhat.sites --loc /data/{wdirpop}/ldhat/{dataset}.{chrom}.ldhat.locs --lk /data/{wdirpop}/ldhat/{dataset}.lookup.{chrom}.new_lk.txt --res /data/{wdirpop}/ldhat/{dataset}.{chrom}.bpen{bpen}.res.txt --nsim $nsim --out /data/{wdirpop}/ldhot/{dataset}.{chrom}.bpen{bpen} --hotdist {config[ldhot.hotdist]} --seed {config[seed]}
        singularity exec --bind $PWD:/data ldhot.sif ldhot_summary --res /data/{wdirpop}/ldhat/{dataset}.{chrom}.bpen{bpen}.res.txt --hot /data/{wdirpop}/ldhot/{dataset}.{chrom}.bpen{bpen}.hotspots.txt --out /data/{wdirpop}/ldhot/{dataset}.{chrom}.bpen{bpen} --sig {config[ldhot.sig]} --sigjoin {config[ldhot.sigjoin]} 
	gzip -f {wdirpop}/ldhat/{dataset}.{chrom}.bpen{bpen}.res.txt
	gzip -f {wdirpop}/ldhot/{dataset}.{chrom}.bpen{bpen}.hot_summary.txt
	gzip -f {wdirpop}/ldhot/{dataset}.{chrom}.bpen{bpen}.hotspots.txt
        """

