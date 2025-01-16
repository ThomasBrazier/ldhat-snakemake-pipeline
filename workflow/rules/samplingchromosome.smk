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
        "../envs/vcftools.yaml"
    shell:
        """
        vcftools --gzvcf {wdirpop}/{dataset}.pop.vcf.gz --out {wdirpop}/{dataset}.chromosome.{chrom} --recode --chr {chrom} --max-missing {config[maxmissing]} --min-alleles 2 --max-alleles 2 --hwe {config[hwe]}
        mv {wdirpop}/{dataset}.chromosome.{chrom}.recode.vcf {wdirpop}/{dataset}.chromosome.{chrom}.vcf
        bgzip -f {wdirpop}/{dataset}.chromosome.{chrom}.vcf
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
        "../envs/shapeit.yaml"
    shell:
        """
        # Remove --thread {config[cores]} if causing errors
        # shapeit --input-vcf {wdirpop}/{dataset}.chromosome.{chrom}.vcf.gz --output-max {wdirpop}/{dataset}.phased.chromosome.{chrom} --effective-size $(cat {wdirpop}/statistics/{dataset}.effective_size) --window {config[shapeitWindow]} --thread {config[cores]} --output-log {wdirpop}/logs/{dataset}.chromosome.{chrom}.shapeit.log --force
        shapeit --input-vcf {wdirpop}/{dataset}.chromosome.{chrom}.vcf.gz --output-max {wdirpop}/{dataset}.phased.chromosome.{chrom} --effective-size $(cat {wdirpop}/statistics/{dataset}.effective_size) --window {config[shapeitWindow]} --output-log {wdirpop}/logs/{dataset}.chromosome.{chrom}.shapeit.log --force
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
        "../envs/vcftools.yaml"
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
        "../envs/vcftools.yaml"
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
        vcf = "{wdirpop}/{dataset}.chromosome.{chrom}.phased.vcf.gz",
        snpden = "{wdirpop}/mask/{dataset}.chromosome.{chrom}.snpden"
    output:
        "{wdirpop}/{dataset}.chromosome.{chrom}.pseudodiploid.vcf.gz"
    log:
        "{wdirpop}/logs/{dataset}.chromosome.{chrom}.pseudodiploid.log"
    conda:
        "../envs/Renv.yaml"
    shell:
        """
        if [ {config[pseudodiploid]} -eq 1 ]; then
            Rscript workflow/scripts/pseudodiploids.R {wdirpop} {chrom} 1
        else
            if [ {config[pseudodiploid]} -eq 2 ]; then
                Rscript workflow/scripts/pseudodiploids.R {wdirpop} {chrom} 2
            else
                cp {input.vcf} {output}
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
        "../envs/vcftools.yaml"
    shell:
        """
        gunzip {wdirpop}/{dataset}.chromosome.{chrom}.pseudodiploid.vcf.gz
        bgzip {wdirpop}/{dataset}.chromosome.{chrom}.pseudodiploid.vcf
        tabix -f -p vcf {wdirpop}/{dataset}.chromosome.{chrom}.pseudodiploid.vcf.gz --csi
        """