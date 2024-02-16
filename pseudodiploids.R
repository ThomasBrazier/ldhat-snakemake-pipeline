#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
# <wdirpop> must be passed in argument

# Reshuffling code:
# 1 - Discard one of the two haplotypes in each individual before reshuffling
# 2- Reshuffle but keep all haplotypes

# test if there are three arguments: if not, return an error
if (length(args) != 3) {
 stop("<wdirpop>, <chromosome> and the reshuffling code (1 or 2) must be passed in argument", call.=FALSE)
}
wdirpop = args[1]
chromosome = args[2]
reshuf = args[3]
# wdirpop = "data/Oryza_sativa_McCouch2016/K3.pop2"
# chromosome = "1"
require(vcfR)

vcf_file = list.files(path = wdirpop, pattern = paste("chromosome.", chromosome, ".phased.vcf.gz$", sep = ""), full.names = TRUE)
vcf_prefix = gsub(".phased.vcf.gz", "", vcf_file)

cat("Loading vcf file", vcf_file,".\n")
vcf = read.vcfR(vcf_file, verbose = FALSE, convertNA = FALSE)
print(vcf)
genotypes = as.matrix(vcf@gt)
# Extract haplotypes in a matrix
GT = genotypes[,1]
genotypes = genotypes[,-1]
cat("Random resampling of haplotypes.\n")
# Matrix of first haplotype - phased data
sep = "|"
matpos1 = apply(genotypes, c(1,2), function(x){strsplit(x, sep)[[1]][1]})
# Matrix of second haplotype - phased data
matpos2 = apply(genotypes, c(1,2), function(x){strsplit(x, sep)[[1]][3]})
# Make pseudodiploids by randomly pairing haplotypes
set.seed(42)
if (reshuf == 1) {
  # Keep only one haplotype per individual -> matpos1
  if (ncol(matpos1) %% 2 != 0) {
    matpos1 = matpos1[,-1]
  }
  idx_sample = sample(seq(1, ncol(matpos1)/2), replace = FALSE)
  # Resample individuals to make pseudodiploids from two haploid individuals
  new_genotypes = paste(matpos1[,idx_sample],
                        matpos1[,-idx_sample], sep = sep)
  new_id_sample = paste(colnames(matpos1)[idx_sample], colnames(matpos1)[-idx_sample], sep = "_")
  new_genotypes = matrix(new_genotypes, nrow = nrow(matpos1), ncol = ncol(matpos1)/2,
                         dimnames = list(seq(1, nrow(matpos1)), new_id_sample))
}
if (reshuf == 2) {
  new_genotypes = paste(matpos1[sample(seq(1, nrow(matpos1)), replace = FALSE),],
                        matpos2[sample(seq(1, nrow(matpos2)), replace = FALSE),], sep = sep)
  new_genotypes = matrix(new_genotypes, nrow = nrow(genotypes), ncol = ncol(genotypes),
                         dimnames = list(seq(1, nrow(genotypes)), colnames(genotypes)))
}
new_genotypes = cbind(GT, new_genotypes)
colnames(new_genotypes)[1] = "FORMAT"

vcf@gt = new_genotypes
print(vcf)

print(matrix(new_genotypes[1:5,1:10], nrow = 5, ncol = 10))

cat("Write new vcf file.\n")
# Write a new pseudodiploid vcf
out_file = paste(vcf_prefix, ".pseudodiploid.vcf.gz", sep = "")
write.vcf(vcf, out_file)

