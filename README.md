# A Snakemake pipeline to estimate LD-based recombination maps

**Authors: Thomas Brazier<sup>1</sup> and Sylvain Glémin<sup>1</sup>**

*Institutions: (1) UMR 6553 ECOBIO, University of Rennes*


This pipeline is a BETA version in active development (use the `main` branch as `dev` is regularly breaking things). It is fully functional but not error-prone and not optimized (e.g. better use of parallelism). Please report bugs and errors for improvement.


## Install

### Dependencies

Strict dependencies must be installed before first run:
- conda
- snakemake >= 8
- snakemake-executor-plugin-slurm
- singularity (> 3.0)

### Installation

After installing dependencies, clone the github directory where you want to perform computations. This github directory will be your working directory.

```
git clone https://github.com/ThomasBrazier/LDRecombinationMaps-pipeline.git 
```

Then, install necessary softwares with the `install.sh` script.


## Input

Required input files are:

* '<dataset>.vcf.gz' in bgzipped format
* '<dataset>.chromosomes' is a list of chromosome names (one name per row)
* a 'samplelist' file (optional) if you wish to subset a list of individuals from the original dataset



## Usage

Put your `<dataset>` directory into teh working directory `./data`.
By default, working directory is `data/`. To run in a different directory, change the value in `config.yaml` or in command line.

You can configure your pipeline settings in the `config.yaml` file and the cluster profile (e.g. slurm settings for time, memory and CPU) in `profiles/slumr/config.yaml`.


You need a first step of data preprocessing to infer the number of independent genetic populations in the sample (based on FastStructure [[1]](#1)) and output summary statistics to choose the most appropriate population in further analyses. It produces a file `structure/poplist.*.*` which contains the list of individuals to sample in the main pipeline.

Launch a job on the cluster: 

```
sbatch launch_preprocessing_structure.sh <dataset> <chromosome>
```

or alternatively on CLI:

```
ncores=<number of cores to use>
snakemake -s workflow/Snakefile --use-conda --use-singularity --cores $ncores -j $ncores --until k_statistics --config dataset=<dataset> chromosome=<chromosome>
```


Once population structure is inferred from 'popstatistics.<K>', 'structure/chooseK' and 'structure/distruct.<K>.svg', run the main pipeline for a <dataset> and a single <chromosome> after specifying the chosen <K> number of genetic clusters to consider and the <population> to sample in your config.yaml.

Launch a job on the cluster: 

```
sbatch launch.sh <dataset> <chromosome>
```

or alternatively on CLI:

```
snakemake -s workflow/Snakefile --use-conda --use-singularity --cores $ncores --config dataset=<dataset> --K <K> --pop <pop> --chrom <chromosome>
```

You can also decide that you are not interested in recombination hotspot detection

```
snakemake -s workflow/Snakefile --use-conda --use-singularity --cores $ncores --until Rmd_report --config dataset=<dataset> --K <K> --pop <pop> --chrom <chromosome>
```


Alternatively, you can use the Slurm launcher script and modify it for your custom needs.

```
./launch.sh <dataset> <chromosome> # Run in CLI
# or 
sbatch launch.sh <dataset> <chromosome> # Run in a SLURM job
```

At the current stage, you can run as many <dataset> as you want in parallel, as directories are isolated, but only one <chromosome> at a time to avoid interferences beween Snakemake parallel processes accessing the same files. This issue is on a list of future improvements.
Only one population can be sampled in a sample directory. For analysing more than one population, duplicate the sample directory.


## Details of the main pipeline


### Data trimming

VCF trimming is done by vcftools.

- maf
- missing
- maxmissing
- minQ. Filter with vcftools based on a minimum Quality value per site. Set `minQ` value to 0 if your vcf does not have Quality values.

### Phasing with ShapeIt

After subsampling the population, the genotypes are phased with ShapeIt2 [[2]](#2). 

Verify that contig length are annotated in the vcf file header as they are necessary at the phasing step.

### Making Pseudodiploids (for selfing species)

Selfing individuals can be highly homozygotes. Thus diploid genomes are exact duplicates. If you want to sample only one phased genome for those species (one haploid genome), you can set the `pseudodiploid: 1` option, otherwise set it to `pseudodiploid: 0`.

If `pseudodiploid: 1`, only one phased haploid genome per individual will be sampled and two individuals will be mixed for compatibility with `LDhat` which accepts diploid vcf only. Twice as much individuals will be sampled to gather the "pseudodiploid" vcf.

For example, if `pseudodiploid: 1`, 80 individuals are sampled and only one haploid genome of each one of them are sampled to assemble a diploid vcf with 40 diploid mixed individuals, hence false diploids.


### LDhat

Recombination rates are estimated with `LDhat` [[3]](#3).

Theta must be specified in the `config.yaml` file. The look-up table will be computed from a pre-generated one with the `lkgen` function of LDhat. A maximum of 100 haploid individuals is allowed (50 diploids). If you specify a `theta` different of 0.01 or 0.001, a complete look-up table will be computed, which require extra time. 


### The Large sample sub-pipeline for large numbers of SNPs

LDhat requires large SNP datasets but processing millions of SNPs can be computationally limiting. To circumvent this problem, you can run the pipeline with optimized steps to cut LDhat processes in chunks of thoudands of SNPs (default 2,000 SNPs overlapping by 200 SNPs). Use the option `large_sample: "yes"`.

### LDhot

Recombination hotspots are inferred from results of LDhat with LDhot [[4]]("4").


### Parallelisation

Some steps are multithreaded with GNU parallel. Do not forget to set up the number of cores properly.


## Input Files


The input files required for the pipeline are:

* `<dataset>.vcf.gz`, a tabix vcf file, bgzipped
* `samplelist`, a one column text file with a list of individuals to keep from the original vcf


## Output Files


Output files of LDhat and LDhot are placed in `data/<dataset>/K*.pop*/ldhat` and `data/<dataset>/K*.pop*/ldhot`. `ldhat/*.res.txt.gz` contains the recombination rates estimated between the SNP n-1 and the SNP n. `ldhot/*.hotspots.txt.gz` and `ldhot/*.hot_summary.txt.gz` contain the recombination hostpots inferred with LDhot.

Read the LDhat and LDhot manuals for further details on output files.


## Limitations

* Only one chromosome at a time, at least until the rule 'split_chromosome' has been successfully passed one time
* contig length in the vcf file are required for SMC++

## Known bugs





## References

<a id="1">[1]</a> 
Raj, A., Stephens, M., & Pritchard, J. K. (2014).
fastSTRUCTURE: variational inference of population structure in large SNP data sets.
Genetics, 197(2), 573-589.

<a id="2">[2]</a>
Delaneau, O., Zagury, J. F., & Marchini, J. (2013).
Improved whole-chromosome phasing for disease and population genetic studies.
Nature methods, 10(1), 5-6.

<a id="3">[3]</a>
Auton, A., Myers, S., & McVean, G. (2014).
Identifying recombination hotspots using population genetic data.
arXiv preprint arXiv:1403.4264.

<a id="4">[4]</a>
Auton, A., Myers, S., & McVean, G. (2014).
Identifying recombination hotspots using population genetic data.
arXiv preprint arXiv:1403.4264.

