# LDRecombinationMaps_pipeline

## Installation

### Dependencies

Strict dependencies have to be installed before first run:
- conda
- snakemake
- singularity (> 3.0)

To install, just clone the github directory.

```
git clone
```

Singularity images are required. Run in the directory:

```
singularity pull faststructure.sif docker://tombrazier/faststructure
singularity pull smcpp.sif docker://terhorst/smcpp
singularity pull pyrho.sif docker://tombrazier/pyrho
singularity pull ldhat.sif docker://tombrazier/ldhat
```

To install conda environments at first run, use

```
snakemake -s data_preprocessing.snake --use-conda --conda-create-envs-only
snakemake --use-conda --conda-create-envs-only
```

## Usage

Put your `<dataset>` in `data/`.

The pipeline can be run for one/many samples and/or populations. Thus you need to give this information:
* sample name
* population to sample (after preliminary analyses)
* chromosome to sample
Only one population can be sampled in a sample directory. For analysing more than one population, duplicate the sample directory.

The workflow take into consideration the directory where the sample data is stored. Thus working directory must be set in 'config.yaml' when invoking Snakemake.

TO run the first step of the pipeline, invoke the 'data_preprocessing.snake' file to identify population structure in your dataset.

```
ncores=8
snakemake -s data_preprocessing.snake --use-conda --use-singularity --cores $ncores -j $ncores --config dataset=<dataset>
```

By default, working directory is `data/`. To run in a different directory, change the value in `config.yaml` or in command line.

Once population structure is inferred, run the main pipeline after specifying the chosen number of genetic clusters to consider (K) and the population to sample in your config.yaml.

```
snakemake -s Snakefile --use-conda --use-singularity --cores $ncores --config dataset=<dataset> --K <K> --pop <pop> --chrom <chromosome>
```

### Files

* dataset.vcf.gz, a tabix vcf file, bgzipped
* samplelist, a one column text file with a list of individuals to keep in the original vcf


## Options


## Details


## References
