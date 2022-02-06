# LDRecombinationMaps_pipeline

## Installation

### Dependencies

Singularity images are required. Run in the directory

```
singularity pull docker://tombrazier/faststructure
```

To install conda environments at first run, use

```
snakemake -s data_preprocessing.snake --use-conda --use-singularity --cores 1
snakemake --use-conda --use-singularity
```

## Usage

The pipeline can be run for one/many samples and/or populations. Thus you need to give this information:
* sample name
* population to sample (after preliminary analyses)
* chromosome to sample
Only one population can be sampled in a sample directory. For analysing more than one population, duplicate the sample directory.

The workflow take into consideration the directory where the sample data is stored. Thus working directory must be set in 'config.yaml' when invoking Snakemake.

TO run the first step of the pipeline, invoke the 'data_preprocessing.snake' file to identify population structure in your dataset.

```
ncores=8
snakemake -s data_preprocessing.snake --use-conda --use-singularity --cores $ncores
```

Once population structure is inferred, run the main pipeline after specifying the chosen number of gneetic clusters to consider (K) and the population to sample in your config.yaml.

```
snakemake -s Snakefile --use-conda --use-singularity --cores $ncores
```

### Files

* dataset.vcf.gz, a tabix vcf file, bgzipped
* samplelist, a one column text file with a list of individuals to keep in the original vcf


## Options


## Details


## References
