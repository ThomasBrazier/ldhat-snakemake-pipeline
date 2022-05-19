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

Pre-generated look up tables are necessary for LDhat.

```
wget https://github.com/auton1/LDhat/tree/master/lk_files
```


## Usage

Put your `<dataset>` in `data/`.

The pipeline can be run for one/many samples and/or populations. Thus you need to give this information:
* sample name
* population to sample (after preliminary analyses)
* chromosome to sample
Only one population can be sampled in a sample directory. For analysing more than one population, duplicate the sample directory.

The workflow take into consideration the directory where the sample data is stored. Thus working directory must be set in 'config.yaml' when invoking Snakemake.

To run the first step of the pipeline, invoke the 'data_preprocessing.snake' file to identify population structure in your dataset.

```
ncores=8
snakemake -s data_preprocessing.snake --use-conda --use-singularity --cores $ncores -j $ncores --config dataset=<dataset>
```

By default, working directory is `data/`. To run in a different directory, change the value in `config.yaml` or in command line.

Once population structure is inferred from 'popstatistics.<K>', 'structure/chooseK' and 'structure/distruct.<K>.svg', run the main pipeline after specifying the chosen <K> number of genetic clusters to consider and the <population> to sample in your config.yaml.

```
snakemake -s Snakefile --use-conda --use-singularity --cores $ncores --config dataset=<dataset> --K <K> --pop <pop> --chrom <chromosome>
```

One chromosome is run at a time. You must specify which chromosome to process in the 'config.yaml'.

### SMC++


### ShapeIt

Contig length in the vcf file header are necessary at the phasing step. You should verify this features before running analyses and annotate your vcf if necessary.


### LDhat

Theta must be specified in the `config.yaml` file. To date, only theta = 0.01 or theta = 0.001 are allowed. The look-up table will be computed from a pre-generated one with the lkgen function of LDhat. A maximum of 100 haploid individuals is allowed (50 diploids).


### Files

* <dataset>.vcf.gz, a tabix vcf file, bgzipped
* samplelist, a one column text file with a list of individuals to keep in the original vcf


## Options



## Details


## References
