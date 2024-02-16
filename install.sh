#!/bin/bash

# INSTALL SCRIPT

# LDpop must be installed locally from the github directory:
git clone https://github.com/popgenmethods/ldpop.git ldpop/

# In addition, Singularity images are required for additional softwares. Run within the working directory:
singularity pull faststructure.sif docker://tombrazier/faststructure
singularity pull ldhat.sif docker://tombrazier/ldhat:v1.0
singularity pull ldhot.sif docker://tombrazier/ldhot:v1.0
singularity pull smcpp.sif docker://terhorst/smcpp:latest

# To install conda environments at the first run of the pipeline, use
snakemake -s data_preprocessing.snake --use-conda --conda-create-envs-only
snakemake --use-conda --conda-create-envs-only