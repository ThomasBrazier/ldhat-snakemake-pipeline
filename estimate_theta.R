#!/usr/bin/env Rscript
# Estimate theta on a vcf file

args = commandArgs(trailingOnly=TRUE)

# test if there is one argument: if not, return an error
if (length(args) != 1) {
  stop("<set> and <pop> must be passed in argument", call.=FALSE)
}


set = args[1]
pop = args[2]
wdir = paste("data/", set, "/", pop, "/", sep = "")
# set = "Oryza_sativa_McCouch2016"
# pop = "K1.pop1"
require(vcfR)
require(pegas)

vcf_file = list.files(path = wdir, pattern = paste(set, ".pop.vcf.gz$", sep = ""), full.names = TRUE)

cat("Loading vcf file", vcf_file,".\n")
vcf = read.vcfR(vcf_file, verbose = FALSE, convertNA = FALSE)
print(vcf)

require(adegenet)
require(hierfstat)
genind = vcfR2genind(vcf)
# Diversity
div = summary(genind)
basicstat = basic.stats(genind, diploid = TRUE, digits = 2)


vcf = readVCF(vcf_file, tid="1:40", from=1, to=10^9)  # for tid=" 1:17", I want to use SNP from all the 17 chromosomes

get.neutrality(vcf,theta=TRUE,stats=TRUE)


get.neutrality(vcf)[[1]]

# Tajima.D n.segregating.sites Rozas.R_2   Fu.Li.F   Fu.Li.D Fu.F_S Fay.Wu.H Zeng.E Strobeck.S
# 1 - 10291 0.7972621                   1        NA 0.6281656 0.3860124    NaN      NaN    NaN        NaN


# Estimate theta Watterson for each site
s = dim(vcf@gt)[1] # Number of segregating sites at a position
n = dim(vcf@gt)[2] # Number of sequences
th = theta.s(s, n, variance = TRUE)

# Compute the mean Theta and its CI
# Bootstrap


# Save results
write.table(th, file = paste(wdir, "theta.txt", sep = ""),
            quote = F, col.names = F, row.names = F, sep = "\t")





