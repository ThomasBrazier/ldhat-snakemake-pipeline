#!/usr/bin/env Rscript
# args = commandArgs(trailingOnly=TRUE)
# 
# # test if there is at least one argument: if not, return an error
# if (length(args)==0) {
#   stop("At least one argument must be supplied (input file).n", call.=FALSE)
# } else if (length(args)==1) {
#   # default output file
#   args[2] = "out.txt"
# }

# vcf_file = "~/Academic/LDRecombinationMaps-pipeline/data/Oryza_sativa_McCouch2016/Oryza_sativa_McCouch2016.trimmed.vcf.gz"
# vcf_file = args[1]

# LOADING ENVIRONMENT
# Loading packages
library(ade4)
library(adegenet)
library(vcfR)

# Loading variables & objects
# vcf = read.vcfR(vcf_file)
# genlight = vcfR2genlight(vcf)
# rm(vcf)
# gc()

# save(genlight, file = "/home/tbrazier/Academic/LDRecombinationMaps-pipeline/data/Oryza_sativa_McCouch2016/genlight.Rda")
load("~/Academic/LDRecombinationMaps-pipeline/data/Oryza_sativa_McCouch2016/genlight.Rda")

genlight

# Subset will be reduced to a random subset of markers (n=10,000 at each iteration)
# Computation issues
n.loci = 10000
# Explore parameter space K = 1-5
set.seed(42)
sub.genlight = genlight[,sample(1:genlight$n.loc, n.loci)]



k.clusters = find.clusters(sub.genlight, max.n.clust = 5, choose.n.clust = FALSE)

cat("Done")

# Find the best K

# Cluster individuals

# Compute summary statistics

# Save figures and statistics 
# dapc_plot = "{wdir}/popstructure/DAPC.pdf"
# summary_stats = "{wdir}/popstructure/statistics.csv}"

