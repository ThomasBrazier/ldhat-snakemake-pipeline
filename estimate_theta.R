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

cat("Estimating Theta Watterson per site\n")
# Estimate theta Watterson per site
# Count how many segregating sites per position
res = apply(vcf@gt, 1, function(x) {x == "0/1" | x =="1/0" | x == "0|1" | x =="1|0"})
s = apply(res, 2, sum)
# theta per site
# s = dim(vcf@gt)[1] # Number of segregating sites at a position
n = dim(vcf@gt)[2] # Number of sequences
th = theta.s(s, n, variance = FALSE)

cat("Compute the mean Theta and 95% CI\n")
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
cat("Saving\n")
write.table(df, file = paste(wdirpop, "/theta.txt", sep = ""),
            quote = F, col.names = T, row.names = F, sep = "\t")





