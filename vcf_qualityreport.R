#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)
vcf_file = args[1]
library(rmarkdown)
library(vcfR)
library(adegenet)
library(poppr)
rmarkdown::render('vcf_qualityreport.Rmd',
                  params = list(vcf_file = vcf_file))
