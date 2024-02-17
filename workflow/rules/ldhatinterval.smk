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
            "../envs/vcftools.yaml"
        log:
            "{wdirpop}/logs/{dataset}.split_dataset.{chrom}.log"
        shell:
            """
            # The first line splits up the snps into chunks of whatever size you want (-l) and then the next line loops over each file and subsets the vcf according.
            bcftools query -f'%CHROM\t%POS\n' {wdirpop}/{dataset}.chromosome.{chrom}.ldhat.vcf.gz > {wdirpop}/{dataset}.{chrom}.positions 
            python workflow/scripts/split_dataset.py {wdirpop}/{dataset}.{chrom}.positions {wdirpop}/ldhat/{dataset}.{chrom} {config[cut_size]} {config[cut_overlap]}
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
            "../envs/vcftools.yaml"
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
            "../envs/vcftools.yaml"
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
            "../envs/vcftools.yaml"
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
            echo "Loci    Mean_rho    Median    L95    U95" > {wdirpop}/ldhat/{dataset}.{chrom}/bpen{bpen}.header
            Loci="-1.000"
            MeanRho=$(cat {wdirpop}/ldhat/{dataset}.{chrom}/bpen{bpen}.res_noheader.txt | awk '{{s+=$2}} END {{printf "%.0f", s}}')
            Median=$(cat {wdirpop}/ldhat/{dataset}.{chrom}/bpen{bpen}.res_noheader.txt | awk '{{s+=$3}} END {{printf "%.0f", s}}')
            L95=$(cat {wdirpop}/ldhat/{dataset}.{chrom}/bpen{bpen}.res_noheader.txt | awk '{{s+=$4}} END {{printf "%.0f", s}}')
            U95=$(cat {wdirpop}/ldhat/{dataset}.{chrom}/bpen{bpen}.res_noheader.txt | awk '{{s+=$5}} END {{printf "%.0f", s}}')
            echo "$Loci     $MeanRho        $Median $L95    $U95"
            echo "$Loci    $MeanRho    $Median    $L95    $U95" >> {wdirpop}/ldhat/{dataset}.{chrom}/bpen{bpen}.header
            cat {wdirpop}/ldhat/{dataset}.{chrom}/bpen{bpen}.header > {wdirpop}/ldhat/{dataset}.{chrom}/bpen{bpen}.res.txt
            cat {wdirpop}/ldhat/{dataset}.{chrom}/bpen{bpen}.res_noheader.txt >> {wdirpop}/ldhat/{dataset}.{chrom}/bpen{bpen}.res.txt
            cp {wdirpop}/ldhat/{dataset}.{chrom}/bpen{bpen}.res.txt {wdirpop}/ldhat/{dataset}.{chrom}.bpen{bpen}.res.txt
            echo "Concat LDhat rates"
            bash workflow/scripts/concat_ldhat_rates.sh {wdirpop} {dataset} {chrom} {bpen} $bigchunk $smalloverlap $chunk
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
            "../envs/vcftools.yaml"
        shell:
            """
            vcftools --gzvcf {wdirpop}/{dataset}.chromosome.{chrom}.ldhat.vcf.gz --chr {chrom} --ldhat --out {wdirpop}/ldhat/{dataset}.{chrom}.{bpen}
            # TODO Test data integrity
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
            "../envs/vcftools.yaml"
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

