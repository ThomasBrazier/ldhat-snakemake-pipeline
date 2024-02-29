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
        "../envs/Renv.yaml"
    shell:
        """
        gzip -c {wdirpop}/ldhat/{dataset}.{chrom}.bpen{bpen}.res.txt > {wdirpop}/ldhat/{dataset}.{chrom}.bpen{bpen}.res.txt.gz
        Rscript --vanilla workflow/scripts/ldhat_MCMC.R {wdirpop} {dataset} {chrom} {bpen}
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
        "../envs/Renv.yaml"
    shell:
        """ 
        Rscript workflow/scripts/vcf_qualityreport_chrom.R {dataset} {chrom} {wdirpop} {bpen}
        # Copy the .yaml config
        cp {wdir}/config.yaml {wdirpop}/{dataset}.{chrom}.bpen{bpen}.yaml
        """
