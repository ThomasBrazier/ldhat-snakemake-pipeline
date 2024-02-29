#!/usr/bin/env Rscript
# Testing convergence of the MCMC chains

args = commandArgs(trailingOnly=TRUE)

# test if there is one argument: if not, return an error
if (length(args) != 4) {
  stop("<wdirpop>, <dataset>, <chromosome> and <bpen> must be passed in argument", call.=FALSE)
}
wdirpop = args[1]
dataset = args[2]
chromosome = args[3]
bpen = args[4]
library(rmarkdown)
library(vcfR)
library(adegenet)
library(poppr)
rmarkdown::render('workflow/scripts/ldhat_MCMC.Rmd',
                  output_file = paste0(dataset, ".", chromosome, ".bpen", bpen, ".ldhat_MCMC.html"),
                  output_dir = paste0(wdirpop, "/MCMC"),
                  params = list(wdirpop = wdirpop,
                                dataset = dataset,
                                chromosome = chromosome,
                                bpen = bpen))
