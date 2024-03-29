# The global dataset is trimmed for SNPs
# and individuals if a 'samplelist' file is provided
# Otherwise all individuals are kept
rule trimming_vcf:
    """
    A first step to trim every 'sample' dataset to the same quality criteria
    """
    input:
        expand("{wdir}/{dataset}.vcf.gz", wdir=wdir, dataset=dataset)
    output:
    	"{wdir}/trimmed.vcf.gz"
    log:
        "{wdir}/logs/trimming_vcf.log"
    conda:
        "../envs/vcftools.yaml"
    shell:
        """
        if [ -f "{wdir}/samplelist" ];
        then
            if [ {config[minQ]} -eq 0 ]
            then
            vcftools --gzvcf {input} --out {wdir}/out --recode --remove-indels --keep {wdir}/samplelist --maf {config[maf]} --max-missing {config[maxmissing]} --min-alleles 2 --max-alleles 2 
            else
            vcftools --gzvcf {input} --out {wdir}/out --recode --remove-indels --keep {wdir}/samplelist --maf {config[maf]} --max-missing {config[maxmissing]} --min-alleles 2 --max-alleles 2 --minQ {config[minQ]}
            fi
            mv {wdir}/out.recode.vcf {wdir}/trimmed.vcf
            bgzip -f {wdir}/trimmed.vcf
            #rm {wdir}/out.log
        else
            if [ {config[minQ]} -eq 0 ]
            then
            vcftools --gzvcf {input} --out {wdir}/out --recode --remove-indels --maf {config[maf]} --max-missing {config[maxmissing]} --min-alleles 2 --max-alleles 2
            else
            vcftools --gzvcf {input} --out {wdir}/out --recode --remove-indels --maf {config[maf]} --max-missing {config[maxmissing]} --min-alleles 2 --max-alleles 2 --minQ {config[minQ]}
            fi
            mv {wdir}/out.recode.vcf {wdir}/trimmed.vcf
            bgzip {wdir}/trimmed.vcf
            #rm {wdir}/out.log
        fi
        # Collapse duplicate SNPs with identical positions
        bcftools norm -d all {wdir}/trimmed.vcf.gz -o {wdir}/trimmed.vcf
        rm {wdir}/trimmed.vcf.gz
        # Subset N random SNPs (default=100,000)
        N=100000
        cat {wdir}/trimmed.vcf | grep '^#' > {wdir}/new_vcf.vcf
        cat {wdir}/trimmed.vcf | grep -v '^#' | shuf -n $N  >> {wdir}/new_vcf.vcf
        rm {wdir}/trimmed.vcf
        mv {wdir}/new_vcf.vcf {wdir}/trimmed.vcf
        bgzip {wdir}/trimmed.vcf
        """


# Use vcftools to export to .ped format with the --plink switch.  Then convert the .ped to .bed using Plink.  Faststructure will read the .bed format files (there are three files for each project)
rule vcf2structure:
    """
    Convert the vcf file to bed format for FastStructure
    """
    input:
        "{wdir}/trimmed.vcf.gz"
    output:
        bed = "{wdir}/vcf2structure.bed",
        bam = "{wdir}/vcf2structure.fam",
        bim = "{wdir}/vcf2structure.bim"
    log:
        "{wdir}/logs/vcf2structure.log"
    conda:
        "../envs/vcftools.yaml"
    shell:
        """
        plink --vcf {input} --out {wdir}/vcf2structure --allow-extra-chr --vcf-half-call missing --const-fid --chr-set $(wc -l {wdir}/{dataset}.chromosomes | awk '{{ print $1 }}')
        """


rule faststructure:
    """
    FastStructure
    Detects population structure and infer genetic clusters
    Produces diagnostic plots and summary statistics
    to help selecting the best population for recombination map
    """
    input:
    	bed = "{wdir}/vcf2structure.bed",
        bam = "{wdir}/vcf2structure.fam",
        bim = "{wdir}/vcf2structure.bim"
    output:
    	"{wdir}/structure/faststructure.{k}.log",
        "{wdir}/structure/faststructure.{k}.meanP",
        "{wdir}/structure/faststructure.{k}.varP",
        "{wdir}/structure/faststructure.{k}.meanQ",
        "{wdir}/structure/faststructure.{k}.varQ"
    log:
	    "{wdir}/logs/faststructure.{k}.log"
    shell:
        """
        singularity exec --bind {wdir}:/faststructure/data faststructure.sif python /faststructure/structure.py -K {wildcards.k} --input=/faststructure/data/vcf2structure --output=/faststructure/data/structure/faststructure --full --seed=100
        """


rule distruct:
    """
    Distruct Plots
    """
    input:
        expand("{wdir}/structure/faststructure.{k}.log", wdir=wdir, k=range(1,maxk+1)),
        expand("{wdir}/structure/faststructure.{k}.meanP", wdir=wdir, k=range(1,maxk+1)),
        expand("{wdir}/structure/faststructure.{k}.varP", wdir=wdir, k=range(1,maxk+1)),
        expand("{wdir}/structure/faststructure.{k}.meanQ", wdir=wdir, k=range(1,maxk+1)),
        expand("{wdir}/structure/faststructure.{k}.varQ", wdir=wdir, k=range(1,maxk+1))
    output:
        expand("{wdir}/structure/distruct.{k}.svg", wdir=wdir, k=range(1,maxk+1))
    log:
        expand("{wdir}/logs/distruct.{k}.log", wdir=wdir, k=range(1,maxk+1))
    shell:
        """
        for k in {{1..{maxk}}}
        do
            singularity exec --bind {wdir}:/faststructure/data faststructure.sif python /faststructure/distruct.py -K $k --input=/faststructure/data/structure/faststructure --output=/faststructure/data/structure/distruct.$k.svg
        done
        """


rule choose_k:
    """
    Choose the best K value
    """
    input:
        expand("{wdir}/structure/distruct.{k}.svg", wdir=wdir, k=range(1,maxk+1))
    output:
        "{wdir}/structure/chooseK"
    log:
        "{wdir}/logs/chooseK"
    shell:
        """
        echo $(singularity exec --bind {wdir}:/faststructure/data faststructure.sif python /faststructure/chooseK.py --input=/faststructure/data/structure/faststructure) > {wdir}/structure/chooseK
        """


rule indlist:
    """
    The list of individuals
    """
    input:
        "{wdir}/structure/chooseK"
    output:
        "{wdir}/indlist"
    log:
        "{wdir}/logs/indlist.log"
    conda:
         "../envs/vcftools.yaml"
    shell:
        """
        bcftools query -l {wdir}/trimmed.vcf.gz > {wdir}/indlist
        """



rule k_statistics:
    """
    Report summary statistics for each population in each K chosen
    to help identify the best config K/population
    Pop. statistic csv file contains columns K, pop number and then columns of statistics
    Headers, statistics are:
    - number of individuals
    - missing data per locus (proportion)
    - missing data per individual (proportion)
    - He
    - Ho
    - Fis
    """
    input:
        "{wdir}/indlist"
    output:
        "{wdir}/structure/{dataset}.popstatistics.{k}.pop1.sites.pi",
        "{wdir}/structure/{dataset}.popstatistics.{k}.pop1.imiss",
        "{wdir}/structure/{dataset}.popstatistics.{k}.pop1.lmiss",
        "{wdir}/structure/{dataset}.popstatistics.{k}.pop1.het",
	    "{wdir}/{dataset}.popstatistics.{k}"
    log:
        "{wdir}/logs/{dataset}.popstatistics.{k}.log"
    conda:
        "../envs/vcftools.yaml"
    shell:
        """
        echo K pop n pi missingindv missingsites Ho He F > {wdir}/{dataset}.popstatistics.{wildcards.k}
        perl -ane '$r = 0; for my $i (1 .. $#F) {{$r = $i if $F[$i] > $F[$r];}} print $r + 1, " ";' < {wdir}/structure/faststructure.{wildcards.k}.meanQ > {wdir}/structure/clusters.{wildcards.k}.tmp
        sed -i -e "s/\\s\\+/\\n/g" {wdir}/structure/clusters.{wildcards.k}.tmp
        paste {wdir}/indlist {wdir}/structure/clusters.{wildcards.k}.tmp > {wdir}/structure/clusters2.{wildcards.k}.tmp
        paste {wdir}/structure/clusters2.{wildcards.k}.tmp {wdir}/structure/faststructure.{wildcards.k}.meanQ > {wdir}/structure/clusters.{wildcards.k}.list
        for pop in {{1..{wildcards.k}}}; do
        awk -v pop=$pop '{{if($2==pop) print $1}}' {wdir}/structure/clusters.{wildcards.k}.list > {wdir}/structure/poplist.{wildcards.k}.$pop
    	n=$(wc -l {wdir}/structure/poplist.{wildcards.k}.$pop | awk '{{ print $1 }}')
        if [ $n -gt 1 ]
		then
		vcftools --gzvcf {wdir}/trimmed.vcf.gz --keep {wdir}/structure/poplist.{wildcards.k}.$pop --site-pi --out {wdir}/structure/{dataset}.popstatistics.{wildcards.k}.pop$pop
        	pi=$(sed '1d' {wdir}/structure/{dataset}.popstatistics.{wildcards.k}.pop$pop.sites.pi | awk '{{ total += $3 }} END {{ print total/NR }}')
        	vcftools --gzvcf {wdir}/trimmed.vcf.gz --keep {wdir}/structure/poplist.{wildcards.k}.$pop --missing-indv --out {wdir}/structure/{dataset}.popstatistics.{wildcards.k}.pop$pop
        	missindv=$(sed '1d' {wdir}/structure/{dataset}.popstatistics.{wildcards.k}.pop$pop.imiss | awk '{{ total += $5 }} END {{ print total/NR }}')
		vcftools --gzvcf {wdir}/trimmed.vcf.gz --keep {wdir}/structure/poplist.{wildcards.k}.$pop --missing-site --out {wdir}/structure/{dataset}.popstatistics.{wildcards.k}.pop$pop
        	misssites=$(sed '1d' {wdir}/structure/{dataset}.popstatistics.{wildcards.k}.pop$pop.lmiss | awk '{{ total += $6 }} END {{ print total/NR }}')
		vcftools --gzvcf {wdir}/trimmed.vcf.gz --keep {wdir}/structure/poplist.{wildcards.k}.$pop --het --out {wdir}/structure/{dataset}.popstatistics.{wildcards.k}.pop$pop
        	Ho=$(sed '1d' {wdir}/structure/{dataset}.popstatistics.{wildcards.k}.pop$pop.het | awk '{{ total += $2/ $4 }} END {{ print total/NR }}')
		He=$(sed '1d' {wdir}/structure/{dataset}.popstatistics.{wildcards.k}.pop$pop.het | awk '{{ total += $3 / $4 }} END {{ print total/NR }}')
		F=$(sed '1d' {wdir}/structure/{dataset}.popstatistics.{wildcards.k}.pop$pop.het | awk '{{ total += $5 }} END {{ print total/NR }}')
		echo {wildcards.k} $pop $n $pi $missindv $misssites $Ho $He $F >> {wdir}/{dataset}.popstatistics.{wildcards.k}
       	else
		touch {wdir}/structure/{dataset}.popstatistics.{wildcards.k}.pop$pop.sites.pi
		touch {wdir}/structure/{dataset}.popstatistics.{wildcards.k}.pop$pop.imiss
		touch {wdir}/structure/{dataset}.popstatistics.{wildcards.k}.pop$pop.lmiss
		touch {wdir}/structure/{dataset}.popstatistics.{wildcards.k}.pop$pop.het
        fi
	    done
        """
