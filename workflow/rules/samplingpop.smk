rule poplist:
    """
    Sample individuals from a given genetic cluster (number of K to retain and name of the cluster) specified in config.yaml
    """
    input:
        expand("{wdir}/{dataset}.vcf.gz", wdir=wdir, dataset=dataset),
        expand("{wdir}/structure/faststructure.{k}.meanP", wdir=wdir, k=k_structure),
        expand("{wdir}/structure/{dataset}.popstatistics.{k}.pop1.sites.pi", wdir=wdir, dataset=dataset, k=k_structure),
        expand("{wdir}/structure/{dataset}.popstatistics.{k}.pop1.imiss", wdir=wdir, dataset=dataset, k=k_structure),
        expand("{wdir}/structure/{dataset}.popstatistics.{k}.pop1.lmiss", wdir=wdir, dataset=dataset, k=k_structure),
        expand("{wdir}/structure/{dataset}.popstatistics.{k}.pop1.het", wdir=wdir, dataset=dataset, k=k_structure),
	    expand("{wdir}/{dataset}.popstatistics.{k}", wdir=wdir, dataset=dataset, k=k_structure)
    output:
        "{wdirpop}/poplist",
        temporary("{wdirpop}/structure/cluster.tmp"),
        temporary("{wdirpop}/structure/cluster.tmp2"),
        temporary("{wdirpop}/structure/cluster")
    log:
        "{wdirpop}/logs/poplist.log"
    conda:
        "../envs/vcftools.yaml"
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
        "{wdirpop}/poplist",
        {wdir}/{dataset}.vcf.gz
    output:
        "{wdirpop}/{dataset}.pop.vcf.gz"
    log:
        "{wdirpop}/logs/{dataset}.sampling_pop.log"
    threads: workflow.cores
    conda:
        "../envs/vcftools.yaml"
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
        "../envs/vcftools.yaml"
    params:
        windowsize=100000
    shell:
        """
        vcftools --gzvcf {input} --window-pi {params.windowsize} --out {wdirpop}/statistics/vcftools
        eff_size=$(awk '{{sum+=$5}} END {{sum=(sum/NR)/(4*{config[mu]}); print sum}}' {wdirpop}/statistics/vcftools.windowed.pi)
        echo ${{eff_size/.*}} > {output}
        """

