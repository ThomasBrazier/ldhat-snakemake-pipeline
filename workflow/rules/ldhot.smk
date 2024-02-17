
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

