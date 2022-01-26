# LDRecombinationMaps_pipeline

## Installation

### Dependencies



## Usage

The pipeline can be run for one/many samples and/or populations. Thus you need to give this information:
* sample name
* population to sample (after preliminary analyses)
* chromosome to sample
Only one population can be sampled in a sample directory. For analysing more than one population, duplicate the sample directory.

The workflow take into consideration the directory where the sample data is stored. Thus working directory must be set when invoking Snakemake.

```
sample=sample
wd='/data/dir'
pop=pop
chrom=chrom
snakemake --directory $wd --config sample=$sample pop=$pop chrom=$chrom
```



## Options


## Details


## References
