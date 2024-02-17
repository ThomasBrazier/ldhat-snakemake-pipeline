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
        "../envs/vcftools.yaml"
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
        "../envs/vcftools.yaml"
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
        "../envs/pyrho.yaml"
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
