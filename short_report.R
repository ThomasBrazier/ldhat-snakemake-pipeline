#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)
dataset = args[1]
chromosome = args[2]
bpen = args[3]
wdirpop = args[4]

library(rmarkdown)
library(vcfR)
library(adegenet)
library(poppr)

rmarkdown::render('short_report.Rmd',
                  params = list(dataset = dataset,
                                chromosome = chromosome,
                                bpen = bpen,
                                wdirpop = wdirpop))
