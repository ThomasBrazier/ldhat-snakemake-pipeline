#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)
dataset = args[1]
chrom = args[2]
wdirpop = as.character(args[3])
bpen = as.character(args[4])

library(rmarkdown)
library(vcfR)
library(adegenet)
library(poppr)
rmarkdown::render('workflow/scripts/vcf_qualityreport_chrom.Rmd',
                  output_file = paste0(dataset, ".", chrom, ".bpen", bpen, ".quality.html"),
                  output_dir = paste0(wdirpop),
                  params = list(dat = dataset, chrom = chrom, wdirpop = wdirpop))