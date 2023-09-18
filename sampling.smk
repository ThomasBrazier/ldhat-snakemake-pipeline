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
