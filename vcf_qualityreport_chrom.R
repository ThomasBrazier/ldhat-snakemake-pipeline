#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)
vcf_file = args[1]
chrom = args[2]
library(rmarkdown)
library(vcfR)
library(adegenet)
library(poppr)
rmarkdown::render('vcf_qualityreport_chrom.Rmd',
                  params = list(dat = vcf_file, chrom = chrom))
