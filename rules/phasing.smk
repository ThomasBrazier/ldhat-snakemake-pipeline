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
        RANDOM={config[seed]}
        vcftools --gzvcf {wdirpop}/{dataset}.chromosome.{chrom}.pseudodiploid.vcf.gz --out {wdirpop}/out --recode --max-indv {config[subset]} --maf {config[maf]} --max-missing {config[maxmissing]}
        mv {wdirpop}/out.recode.vcf {wdirpop}/{dataset}.chromosome.{chrom}.ldhat.vcf
        bgzip -f {wdirpop}/{dataset}.chromosome.{chrom}.ldhat.vcf
        """
