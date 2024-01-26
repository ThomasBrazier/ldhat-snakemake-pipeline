#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)
vcf_file = args[1]
chrom = args[2]
library(rmarkdown)
library(vcfR)
library(adegenet)
library(poppr)
rmarkdown::render('vcf_qualityreport_chrom.Rmd',
                  output_file = paste0(dataset, ".", chromosome, ".bpen", bpen, ".quality.html"),
                  output_dir = paste0(wdirpop),
                  params = list(dat = vcf_file, chrom = chrom))