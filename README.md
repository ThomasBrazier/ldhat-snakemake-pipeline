# LDRecombinationMaps_pipeline

## Installation

### Dependencies

Singularity images are required. Run in the directory

```
singularity pull docker://tombrazier/faststructure
```

To install conda environments at first run, use

```
snakemake --s data_preprocessing.smk --use-conda --use-singularity
snakemake --use-conda --use-singularity
```

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
### Files

* dataset.vcf.gz, a tabix vcf file, bgzipped
* samplelist, a one column text file with a list of individuals to keep in the original vcf


## Options


## Details


## References
