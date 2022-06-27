#!/usr/bin/env Rscript
# Estimate theta on a vcf file

args = commandArgs(trailingOnly=TRUE)

# test if there is one argument: if not, return an error
if (length(args) != 1) {
  stop("<wdirpop> must be passed in argument", call.=FALSE)
}


wdirpop = args[1]

require(vcfR)
require(pegas)

# vcf_file = "data/Oryza_sativa_McCouch2016/Oryza_sativa_McCouch2016_test.vcf.gz"
vcf_file = list.files(path = wdirpop, pattern = paste(".pop.vcf.gz$", sep = ""), full.names = TRUE)

cat("Loading vcf file", vcf_file,".\n")
vcf = read.vcfR(vcf_file, verbose = FALSE, convertNA = FALSE)
print(vcf)

# Estimate theta Watterson per site
# Count how many segregating sites per position
geno = vcf@gt[,-1]
pos = matrix(unlist(apply(geno, c(1,2), function(x) {grepl("0/1", x) | grepl("1/0", x)})),
                nrow = dim(geno)[1], ncol = dim(geno)[2])
s = apply(pos, 2, sum)
# theta per site
# s = dim(vcf@gt)[1] # Number of segregating sites at a position
n = dim(vcf@gt)[2] # Number of sequences
th = theta.s(s, n, variance = FALSE)
th

# Compute the mean Theta and its CI
# Bootstrap
nboot = 1000
boot = numeric(nboot)
for (i in 1:nboot) {
  boot[i] = mean(sample(th, replace = TRUE), na.rm = TRUE)
}
df = data.frame(mean.theta = mean(boot, na.rm = TRUE),
                lower.theta = quantile(boot, 0.025, na.rm = TRUE),
                upper.theta = quantile(boot, 0.975, na.rm = TRUE))

# Save results
write.table(df, file = paste(wdirpop, "theta.txt", sep = ""),
            quote = F, col.names = T, row.names = F, sep = "\t")





