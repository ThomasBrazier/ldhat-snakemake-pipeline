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
    threads: workflow.cores
    conda:
        "envs/vcftools.yaml"
    shell:
        """
        # Filter chromosomes and keep only bi-allelic alelles
        if [ {config[minQ]} -eq 0 ]
        then
        vcftools --gzvcf {wdir}/{dataset}.vcf.gz --out {wdirpop}/{dataset}.pop --recode --keep {wdirpop}/poplist --maf {config[maf]} --max-missing {config[maxmissing]} --min-alleles 2 --max-alleles 2
        else
        vcftools --gzvcf {wdir}/{dataset}.vcf.gz --out {wdirpop}/{dataset}.pop --recode --keep {wdirpop}/poplist --maf {config[maf]} --max-missing {config[maxmissing]} --min-alleles 2 --max-alleles 2 --minQ {config[minQ]}
        fi
        mv {wdirpop}/{dataset}.pop.recode.vcf {wdirpop}/{dataset}.pop.vcf
        bgzip -f {wdirpop}/{dataset}.pop.vcf
        bcftools norm -d all {wdirpop}/{dataset}.pop.vcf.gz -o {wdirpop}/{dataset}.pop.vcf
        bgzip -f {wdirpop}/{dataset}.pop.vcf
        tabix -f --csi {wdirpop}/{dataset}.pop.vcf.gz
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
    threads: workflow.cores
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
    threads: workflow.cores
    log:
        "{wdirpop}/logs/{dataset}.split_chromosome.{chrom}.log"
    conda:
        "envs/vcftools.yaml"
    shell:
        """
        vcftools --gzvcf {wdirpop}/{dataset}.pop.vcf.gz --out {wdirpop}/{dataset}.chromosome.{chrom} --recode --chr {chrom} --max-missing {config[maxmissing]} --min-alleles 2 --max-alleles 2 --hwe {config[hwe]}
        mv {wdirpop}/{dataset}.chromosome.{chrom}.recode.vcf {wdirpop}/{dataset}.chromosome.{chrom}.vcf
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
    threads: workflow.cores
    conda:
        "envs/shapeit.yaml"
    shell:
        """
	    # Remove --thread {config[cores]} if causing errors
        shapeit --input-vcf {wdirpop}/{dataset}.chromosome.{chrom}.vcf.gz --output-max {wdirpop}/{dataset}.phased.chromosome.{chrom} --effective-size $(cat {wdirpop}/statistics/{dataset}.effective_size) --window {config[shapeitWindow]} --thread {config[cores]} --output-log {wdirpop}/logs/{dataset}.chromosome.{chrom}.shapeit.log --force
        #shapeit --input-vcf {wdirpop}/{dataset}.chromosome.{chrom}.vcf.gz --output-max {wdirpop}/{dataset}.phased.chromosome.{chrom} --effective-size $(cat {wdirpop}/statistics/{dataset}.effective_size) --window {config[shapeitWindow]} --output-log {wdirpop}/logs/{dataset}.chromosome.{chrom}.shapeit.log --force
        shapeit -convert --input-haps {wdirpop}/{dataset}.phased.chromosome.{chrom} --output-vcf {wdirpop}/{dataset}.chromosome.{chrom}.phased.vcf --output-log {wdirpop}/logs/{dataset}.chromosome.{chrom}.shapeit.convert.log
        # replace header in vcf to keep information of contig length
        zcat {wdir}/{dataset}.vcf.gz | head -n 1000 | grep '^#' > {wdirpop}/newheader || true
        cat {wdirpop}/newheader | grep -v '^#CHROM' > {wdirpop}/newheader2 || true
        cat {wdirpop}/{dataset}.chromosome.{chrom}.phased.vcf | grep '^#CHROM' > {wdirpop}/colnames || true
        cat {wdirpop}/{dataset}.chromosome.{chrom}.phased.vcf | grep -v '^#' > {wdirpop}/newvcf || true
        cat {wdirpop}/newheader2 {wdirpop}/colnames {wdirpop}/newvcf > {wdirpop}/{dataset}.chromosome.{chrom}.phased.vcf
        bgzip -f {wdirpop}/{dataset}.chromosome.{chrom}.phased.vcf
        #rm {wdirpop}/newheader {wdirpop}/newheader2 {wdirpop}/colnames {wdirpop}/newvcf
        """

rule roh:
    """
    Detect Runs of Homozygosity with PLINK and generate a report of ROHs
    https://zzz.bwh.harvard.edu/plink/ibdibs.shtml#homo
    A simple screen for runs of homozygous genotypes within any one individual is provided by the commands --homozyg-snp and --homozyg-kb which define the run in terms of the required number of homozygous SNPs spanning a certain kb distance, e.g.
    The algorithm is as follows: Take a window of X SNPs and slide this across the genome. At each window position determine whether this window looks 'homozygous' enough (yes/no) (i.e. allowing for some number of hets or missing calls). Then, for each SNP, calculate the proportion of 'homozygous' windows that overlap that position. Call segments based on this metric, e.g. based on a threshold for the average.
    The exact window size and thresholds, relative to the SNP density and expected size of homozygous segments, etc, is obviously important: sensible default values are supplied for the context of dense SNP maps, scanning for large segments. In general, this approach will ensure that otherwise long runs of homozygosity are not broken by the occassional heterozygote. (For more accurate detection of smaller segments, one might consider approaches that also take population parameters such as allele frequency and recombination rate into account, in a HMM approach for example: but for now, PLINK only supports this basic detection of long, homozygous segments).
    """
    input:
        "{wdirpop}/{dataset}.chromosome.{chrom}.phased.vcf.gz"
    output:
        "{wdirpop}/mask/{dataset}.chromosome.{chrom}.hom",
        "{wdirpop}/mask/{dataset}.chromosome.{chrom}.hom.summary",
        "{wdirpop}/mask/{dataset}.chromosome.{chrom}.hom.indiv",
        "{wdirpop}/mask/{dataset}.chromosome.{chrom}.log",
        "{wdirpop}/mask/{dataset}.chromosome.{chrom}.nosex"
    log:
        "{wdirpop}/logs/{dataset}.chromosome.{chrom}.plink_roh.log"
    threads: workflow.cores
    conda:
        "envs/vcftools.yaml"
    shell:
        """
        mkdir -p {wdirpop}/mask
        plink --vcf {wdirpop}/{dataset}.chromosome.{chrom}.phased.vcf.gz --double-id --homozyg --allow-extra-chr --out {wdirpop}/mask/{dataset}.chromosome.{chrom}
        """


rule mask_low_snp_density:
    """
    Detect regions of low SNP density based on a sliding window approach.
    """
    input:
        "{wdirpop}/{dataset}.chromosome.{chrom}.phased.vcf.gz",
        "{wdirpop}/mask/{dataset}.chromosome.{chrom}.hom",
        "{wdirpop}/mask/{dataset}.chromosome.{chrom}.hom.summary",
        "{wdirpop}/mask/{dataset}.chromosome.{chrom}.hom.indiv",
        "{wdirpop}/mask/{dataset}.chromosome.{chrom}.log",
        "{wdirpop}/mask/{dataset}.chromosome.{chrom}.nosex"
    output:
        "{wdirpop}/mask/{dataset}.chromosome.{chrom}.snpden"
    log:
        "{wdirpop}/logs/{dataset}.chromosome.{chrom}.snp_dens.log"
    threads: workflow.cores
    conda:
        "envs/vcftools.yaml"
    shell:
        """
        vcftools --gzvcf {wdirpop}/{dataset}.chromosome.{chrom}.phased.vcf.gz --SNPdensity {config[snpdens.binsize]} --out {wdirpop}/mask/{dataset}.chromosome.{chrom}
        """
	

rule pseudodiploid:
    """
    Make pseudodiploids (phased haplotypes) to take into account homozygotes in high-selfing rates species
    Method 1. Keep only one haplotype per individual and reconstruct diploids
    Method 2. Keep both haplotypes and resample among individuals
    """
    input:
        "{wdirpop}/{dataset}.chromosome.{chrom}.phased.vcf.gz",
        "{wdirpop}/mask/{dataset}.chromosome.{chrom}.snpden"
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
	            cp {wdirpop}/{dataset}.chromosome.{chrom}.phased.vcf.gz {output}
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
    threads: workflow.cores
    conda:
        "envs/vcftools.yaml"
    shell:
        """
        gunzip {wdirpop}/{dataset}.chromosome.{chrom}.pseudodiploid.vcf.gz
        bgzip {wdirpop}/{dataset}.chromosome.{chrom}.pseudodiploid.vcf
        tabix -f -p vcf {wdirpop}/{dataset}.chromosome.{chrom}.pseudodiploid.vcf.gz --csi
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
    threads: workflow.cores
    conda:
        "envs/vcftools.yaml"
    shell:
        """
        RANDOM={config[seed]}
        vcftools --gzvcf {wdirpop}/{dataset}.chromosome.{chrom}.pseudodiploid.vcf.gz --out {wdirpop}/out --recode --max-indv {config[subset]} --maf {config[maf]} --max-missing {config[maxmissing]}
        mv {wdirpop}/out.recode.vcf {wdirpop}/{dataset}.chromosome.{chrom}.ldhat.vcf
        bgzip -f {wdirpop}/{dataset}.chromosome.{chrom}.ldhat.vcf
        """


rule smcpp:
    """
    Estimate a rough demography of the population
    Used to estimate a demography-aware look-up table with LDpop
    """
    input:
        "{wdirpop}/{dataset}.chromosome.{chrom}.ldhat.vcf.gz",
        "{wdirpop}/mask/{dataset}.chromosome.{chrom}.snpden"
    output:
        "{wdirpop}/smcpp/{dataset}.{chrom}/plot.pdf",
        "{wdirpop}/smcpp/{dataset}.{chrom}/model.final.json",
        "{wdirpop}/smcpp/{dataset}.{chrom}.smc.gz",
        "{wdirpop}/smcpp/{dataset}.{chrom}/plot.csv"
    log:
        "{wdirpop}/logs/{dataset}.lookup.{chrom}.log"
    threads: workflow.cores
    conda:
        "envs/vcftools.yaml"
    shell:
        """
        mkdir -p {wdirpop}/smcpp
        zcat {wdirpop}/{dataset}.chromosome.{chrom}.ldhat.vcf.gz | grep '#CHROM' | cut -f 10- | tr '\t' '\n' > {wdirpop}/indlist
        samples=$(cat {wdirpop}/indlist | awk '{{printf("%s,",$0)}}' | sed 's/,\s*$//')
        tabix --csi -f {wdirpop}/{dataset}.chromosome.{chrom}.ldhat.vcf.gz
        singularity exec --bind $PWD:/mnt smcpp.sif smc++ vcf2smc /mnt/{wdirpop}/{dataset}.chromosome.{chrom}.ldhat.vcf.gz /mnt/{wdirpop}/smcpp/{dataset}.{chrom}.smc.gz {chrom} Pop1:$samples -c {config[smcpp.cutoff]}
        # Fit the model using estimate:
        singularity exec --bind $PWD:/mnt smcpp.sif smc++ estimate -o /mnt/{wdirpop}/smcpp/{dataset}.{chrom}/ {config[mu]} /mnt/{wdirpop}/smcpp/{dataset}.{chrom}.smc.gz
        # The model.final.json output file contains fields named rho and N0. rho is the estimated population-scaled recombination rate per base-pair. To convert it to units of generations, multiply by 2 * N0.
        # Visualize the results using plot:
        # -c produces a CSV-formatted table containing the data used to generate the plot.
        singularity exec --bind $PWD:/mnt smcpp.sif smc++ plot -c /mnt/{wdirpop}/smcpp/{dataset}.{chrom}/plot.pdf /mnt/{wdirpop}/smcpp/{dataset}.{chrom}/model.final.json
        # A useful diagnostic for understanding the final output of SMC++ are the sequence of intermediate estimates .model.iter<k>.json which are saved by --estimate in the --output directory. By plotting these, you can get a sense of whether the optimizer is overfitting and requires additional regularization. 
		"""


rule LDpop:
    """
    Generate a complete demography-aware look-up table
    """
    input:
        "{wdirpop}/{dataset}.chromosome.{chrom}.ldhat.vcf.gz",
        "{wdirpop}/smcpp/{dataset}.{chrom}/plot.csv"
    output:
        "{wdirpop}/ldhat/{dataset}.lookup.{chrom}.new_lk.txt"
    log:
        "{wdirpop}/logs/{dataset}.lookup.{chrom}.log"
    threads: workflow.cores
    conda:
        "envs/pyrho.yaml"
    shell:
        """
        Rscript scripts/smcpp_estimates.R {wdirpop}/smcpp/{dataset}.{chrom}/plot.csv {wdirpop}/smcpp/{dataset}.{chrom}/ {config[mu]} {config[theta]}
        s=$(cat {wdirpop}/smcpp/{dataset}.{chrom}/Ne.txt) # coalescent scaled population sizes (s0=present size, sD=ancient size), e.g 100,.1,1
        t=$(cat {wdirpop}/smcpp/{dataset}.{chrom}/times.txt) # times of size changes from present backwards. Must be increasing positive reals. e.g. .5,.58
        #rh=$() # grid of rhos (twice the recomb rate). The grid has num_rh uniformly spaced points from 0 to max_rh, inclusive.
        n=$(zcat {wdirpop}/{dataset}.chromosome.{chrom}.ldhat.vcf.gz | grep ^#CHROM | awk '{{print NF-9}}') # Sample size, number of sequences/haplotypes
        n=$(( 2 * $n )) # Number of haplotypes
        echo "Generating look-up table for $n samples"
        ldpop/run/ldtable.py -n $n -th {config[theta]} -s $s -t $t -rh 101,100 --cores {config[cores]} --approx > {wdirpop}/ldhat/{dataset}.lookup.{chrom}.new_lk.txt
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


    rule interval_stat_split:
        """
        Estimate a recombination landscape with LDhat interval
        Remove temporary files at each iteration
        """
        input:
            "{wdirpop}/ldhat/{dataset}.{chrom}/convert.done"
        output:
            "{wdirpop}/ldhat/{dataset}.{chrom}/stat_bpen{bpen}.done"
        threads: workflow.cores
        log:
            "{wdirpop}/logs/{dataset}.ldhatstat.{chrom}.bpen{bpen}.log"
        conda:
            "envs/vcftools.yaml"
        shell:
            """
            iter={config[interval.iter]}
            samp={config[interval.samp]}
            bpen={config[interval.bpen]}
            burn={config[ldhat.burn]}
            nbatch=$(cat {wdirpop}/ldhat/{dataset}.{chrom}/nbatch)
            echo "nbatch = $nbatch"
            for i in $(seq $nbatch)
            do
            singularity exec --bind $PWD:/data ldhat.sif interval -seq /data/{wdirpop}/ldhat/{dataset}.{chrom}/batch_$i.ldhat.sites -loc /data/{wdirpop}/ldhat/{dataset}.{chrom}/batch_$i.ldhat.locs -lk /data/{wdirpop}/ldhat/{dataset}.lookup.{chrom}.new_lk.txt -its $iter -bpen $bpen -samp $samp -prefix /data/{wdirpop}/ldhat/{dataset}.{chrom}/bpen{bpen}.batch_$i.
            singularity exec --bind $PWD:/data ldhat.sif stat -input /data/{wdirpop}/ldhat/{dataset}.{chrom}/bpen{bpen}.batch_$i.rates.txt -burn $burn -loc /data/{wdirpop}/ldhat/{dataset}.{chrom}/batch_$i.ldhat.locs -prefix /data/{wdirpop}/ldhat/{dataset}.{chrom}/bpen{bpen}.batch_$i.
            rm {wdirpop}/ldhat/{dataset}.{chrom}/bpen{bpen}.batch_$i.new_lk.txt {wdirpop}/ldhat/{dataset}.{chrom}/bpen{bpen}.batch_$i.type_table.txt
            done
            gzip {wdirpop}/ldhat/{dataset}.{chrom}/bpen{bpen}.batch_*.res.txt
            gzip {wdirpop}/ldhat/{dataset}.{chrom}/bpen{bpen}.batch_*.rates.txt
            echo "Done" > {wdirpop}/ldhat/{dataset}.{chrom}/stat_bpen{bpen}.done
            """


    rule concatenate:
        """
        Concatenate .res.txt files
        """
        input:
            "{wdirpop}/ldhat/{dataset}.{chrom}/stat_bpen{bpen}.done"
        output:
            "{wdirpop}/ldhat/{dataset}.{chrom}.bpen{bpen}.res.txt",
            "{wdirpop}/ldhat/{dataset}.{chrom}.bpen{bpen}.rates.txt.gz"
        conda:
            "envs/vcftools.yaml"
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
            zcat bpen{bpen}.batch_1.res.txt.gz | tail -n +3 | head -n $bigchunk > bpen{bpen}.res_noheader.txt || true
            for i in $(seq 2 $(( $nbatch-1 )))
            do
            echo $i
            zcat bpen{bpen}.batch_$i.res.txt.gz | tail -n +3 | head -n $bigchunk | tail -n +$(( $smalloverlap+1 )) >> bpen{bpen}.res_noheader.txt || true
            done
            echo "End of loop on split files."
            zcat bpen{bpen}.batch_$nbatch.res.txt.gz | tail -n +3 | tail -n +$(( $smalloverlap+1 )) >> bpen{bpen}.res_noheader.txt || true
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
            echo "Concat LDhat rates"
            bash concat_ldhat_rates.sh {wdirpop} {dataset} {chrom} {bpen} $bigchunk $smalloverlap $chunk
            echo "Gzip files"
            gzip {wdirpop}/ldhat/{dataset}.{chrom}/bpen{bpen}.rates.txt
            mv {wdirpop}/ldhat/{dataset}.{chrom}/bpen{bpen}.rates.txt.gz {wdirpop}/ldhat/{dataset}.{chrom}.bpen{bpen}.rates.txt.gz
            """


    rule LDhot_convert:
        """
        Reconcile ldhat.locs/ldhat.sites and res.txt for LDhot input
        """
        input:
            expand("{wdirpop}/ldhat/{dataset}.{chrom}.bpen{bpen}.res.txt", wdirpop=wdirpop, dataset=dataset, chrom=chrom, bpen=bpen)
        output:
            "{wdirpop}/ldhat/{dataset}.{chrom}.{bpen}.ldhat.sites",
            "{wdirpop}/ldhat/{dataset}.{chrom}.{bpen}.ldhat.locs"
        log:
            "{wdirpop}/logs/{dataset}.ldhatconvert.{chrom}.{bpen}.log"
        conda:
            "envs/vcftools.yaml"
        shell:
            """
            vcftools --gzvcf {wdirpop}/{dataset}.chromosome.{chrom}.ldhat.vcf.gz --chr {chrom} --ldhat --out {wdirpop}/ldhat/{dataset}.{chrom}.{bpen}
            # Test data integrity
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

    rule interval:
        """
        Estimate a recombination landscape with LDhat
        """
        input:
            "{wdirpop}/ldhat/{dataset}.{chrom}.{bpen}.ldhat.sites",
            "{wdirpop}/ldhat/{dataset}.{chrom}.{bpen}.ldhat.locs"
        output:
            "{wdirpop}/ldhat/{dataset}.{chrom}.bpen{bpen}.new_lk.txt",
            "{wdirpop}/ldhat/{dataset}.{chrom}.bpen{bpen}.bounds.txt",
            "{wdirpop}/ldhat/{dataset}.{chrom}.bpen{bpen}.rates.txt",
            temporary("{wdirpop}/ldhat/{dataset}.{chrom}.bpen{bpen}.type_table.txt")
        log:
            "{wdirpop}/logs/{dataset}.ldhatinterval.{chrom}.bpen{bpen}.log"
        threads: workflow.cores
        shell:
            """
            iter={config[interval.iter]}
            samp={config[interval.samp]}
            bpen={config[interval.bpen]}
            singularity exec --bind $PWD:/data ldhat.sif interval -seq /data/{wdirpop}/ldhat/{dataset}.{chrom}.{bpen}.ldhat.sites -loc /data/{wdirpop}/ldhat/{dataset}.{chrom}.{bpen}.ldhat.locs -lk /data/{wdirpop}/ldhat/{dataset}.lookup.{chrom}.new_lk.txt -its $iter -bpen $bpen -samp $samp -prefix /data/{wdirpop}/ldhat/{dataset}.{chrom}.bpen{bpen}.
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
        threads: workflow.cores
        shell:
            """
            burn={config[ldhat.burn]}
            singularity exec --bind $PWD:/data ldhat.sif stat -input /data/{wdirpop}/ldhat/{dataset}.{chrom}.bpen{bpen}.rates.txt -burn $burn -loc /data/{wdirpop}/ldhat/{dataset}.{chrom}.{bpen}.ldhat.locs -prefix /data/{wdirpop}/ldhat/{dataset}.{chrom}.bpen{bpen}.
            # Compress intermediary files
            gzip -f {wdirpop}/ldhat/{dataset}.{chrom}.bpen{bpen}.rates.txt
            gzip -f {wdirpop}/ldhat/{dataset}.{chrom}.bpen{bpen}.bounds.txt
            """



rule MCMC_report:
    """
    Produce a Rmarkdown/pdf report
    """
    input:
        "{wdirpop}/ldhat/{dataset}.{chrom}.bpen{bpen}.rates.txt.gz",
        "{wdirpop}/ldhat/{dataset}.{chrom}.bpen{bpen}.res.txt",
        "{wdirpop}/ldhat/{dataset}.{chrom}.{bpen}.ldhat.locs"
    output:
        "{wdirpop}/ldhat/{dataset}.{chrom}.bpen{bpen}.res.txt.gz",
        "{wdirpop}/MCMC/{dataset}.{chrom}.bpen{bpen}.ldhat_MCMC.html"
    conda:
        "envs/Renv.yaml"
    shell:
        """
        gzip -c {wdirpop}/ldhat/{dataset}.{chrom}.bpen{bpen}.res.txt > {wdirpop}/ldhat/{dataset}.{chrom}.bpen{bpen}.res.txt.gz
        Rscript --vanilla ldhat_MCMC.R {wdirpop} {dataset} {chrom} {bpen}
        #mv ldhat_MCMC.html {wdirpop}/MCMC/{dataset}.{chrom}.bpen{bpen}.ldhat_MCMC.html
        """


rule Rmd_report:
    """
    Produce a Rmarkdown/pdf report with vcfR - common summary statistics
    """
    input:
        "{wdirpop}/MCMC/{dataset}.{chrom}.bpen{bpen}.ldhat_MCMC.html"
    output:
        "{wdirpop}/{dataset}.{chrom}.bpen{bpen}.quality.html",
        "{wdirpop}/{dataset}.{chrom}.bpen{bpen}.yaml"
    threads: workflow.cores
    conda:
        "envs/Renv.yaml"
    shell:
        """ 
        Rscript vcf_qualityreport_chrom.R {dataset} {chrom} {wdirpop} {bpen}
        #mv vcf_qualityreport_chrom.html {wdirpop}/{dataset}.{chrom}.bpen{bpen}.quality.html
        # Copy the .yaml config
        cp {wdir}/config.yaml {wdirpop}/{dataset}.{chrom}.bpen{bpen}.yaml
        """


rule LDhot:
    """
    Infer recombination hotspots
    LDhot
    """
    input:
        "{wdirpop}/ldhat/{dataset}.{chrom}.{bpen}.ldhat.sites",
        "{wdirpop}/ldhat/{dataset}.{chrom}.{bpen}.ldhat.locs",
        "{wdirpop}/ldhat/{dataset}.{chrom}.bpen{bpen}.res.txt",
        "{wdirpop}/{dataset}.{chrom}.bpen{bpen}.quality.html"
    output:
        "{wdirpop}/ldhot/{dataset}.{chrom}.bpen{bpen}.hotspots.txt.gz",
        temporary("{wdirpop}/ldhot/{dataset}.{chrom}.bpen{bpen}.log")
    threads: workflow.cores
    log:
        "{wdirpop}/logs/{dataset}.{chrom}.bpen{bpen}.ldhot.log"
    shell:
        """
        nsim={config[ldhot.nsim]}
        singularity exec --bind $PWD:/data ldhot.sif ldhot --seq /data/{wdirpop}/ldhat/{dataset}.{chrom}.{bpen}.ldhat.sites --loc /data/{wdirpop}/ldhat/{dataset}.{chrom}.{bpen}.ldhat.locs --lk /data/{wdirpop}/ldhat/{dataset}.lookup.{chrom}.new_lk.txt --res /data/{wdirpop}/ldhat/{dataset}.{chrom}.bpen{bpen}.res.txt --nsim $nsim --out /data/{wdirpop}/ldhot/{dataset}.{chrom}.bpen{bpen} --hotdist {config[ldhot.hotdist]} --seed {config[ldhotseed]}
        singularity exec --bind $PWD:/data ldhot.sif ldhot_summary --res /data/{wdirpop}/ldhat/{dataset}.{chrom}.bpen{bpen}.res.txt --hot /data/{wdirpop}/ldhot/{dataset}.{chrom}.bpen{bpen}.hotspots.txt --out /data/{wdirpop}/ldhot/{dataset}.{chrom}.bpen{bpen} --sig {config[ldhot.sig]} --sigjoin {config[ldhot.sigjoin]} 
	    gzip -f {wdirpop}/ldhot/{dataset}.{chrom}.bpen{bpen}.hot_summary.txt
	    gzip -f {wdirpop}/ldhot/{dataset}.{chrom}.bpen{bpen}.hotspots.txt
        """

