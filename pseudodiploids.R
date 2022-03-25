#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
# <wdirpop> must be passed in argument

# test if there is at least one argument: if not, return an error
if (length(args) != 2) {
 stop("<wdirpop> and <chromosome> must be passed in argument", call.=FALSE)
}
wdirpop = args[1]
chromosome = args[2]
# wdirpop = "data/Oryza_sativa_McCouch2016/K3.pop2"
require(vcfR)

vcf_file = list.files(path = wdirpop, pattern = paste(chromosome, ".phased.csv.gz$", sep = ""), full.names = TRUE)
vcf_prefix = gsub(".phased.vcf.gz", "", vcf_file)

vcf = read.vcfR(vcf_file, verbose = FALSE, convertNA = FALSE)

genotypes = as.matrix(vcf@gt)
# Extract haplotypes in a matrix
GT = genotypes[,1]
genotypes = genotypes[,-1]
# Matrix of first haplotype - phased data
matpos1 = apply(genotypes, c(1,2), function(x){strsplit(x, "|")[[1]][1]})
# Matrix of second haplotype - phased data
matpos2 = apply(genotypes, c(1,2), function(x){strsplit(x, "|")[[1]][2]})
# Make pseudodiploids by randomly pairing haplotypes
set.seed(42)
new_genotypes = paste(matpos1[sample(seq(1, nrow(matpos1)), replace = FALSE),],
                      matpos2[sample(seq(1, nrow(matpos2)), replace = FALSE),], sep = "|")
new_genotypes = matrix(new_genotypes, nrow = nrow(genotypes), ncol = ncol(genotypes),
                       dimnames = list(seq(1, nrow(genotypes)), colnames(genotypes)))
new_genotypes = cbind(GT, new_genotypes)
colnames(new_genotypes)[1] = "FORMAT"

vcf@gt = new_genotypes

# Write a new pseudodiploid vcf
out_file = paste(vcf_prefix, ".pseudodiploid.vcf", sep = "")
write.vcf(vcf, out_file)

