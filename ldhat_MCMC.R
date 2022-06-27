#!/usr/bin/env Rscript
# Estimate theta on a vcf file

args = commandArgs(trailingOnly=TRUE)

# test if there is one argument: if not, return an error
if (length(args) != 4) {
  stop("<wdirpop>, <set>, <chromosome> and <bpen> must be passed in argument", call.=FALSE)
}

wdirpop = args[1]
set = args[2]
chromosome = args[3]
bpen = args[4]

cat("Loading MCMC sampling chain.\n")
mcmc.file = paste(wdirpop, "/ldhat/", set, ".", chromosome, ".bpen", bpen, ".rates.txt.gz", sep = "")
# mcmc.file = "data/Arabidopsis_thaliana_1001genomes.1.bpen15.rates.txt.gz"
mcmc.chain = read.table(gzfile(mcmc.file), header = F, skip = 1)
# One line per position
# Each column is a MCMC iteration
# Do the mean of each column (minus first one) and check if mean posterior is stable across iterations
mean.posterior = data.frame(iteration = 1:(dim(mcmc.chain)[2] - 1),
                            posone = as.numeric(mcmc.chain[1,-1]),
                            pos1000 = as.numeric(mcmc.chain[1000,-1]),
                            posterior = apply(mcmc.chain[,-1], 2, median),
                            totalgenetic = apply(mcmc.chain[,-1], 2, sum))

require(ggplot2)
require(ggpubr)
# Median of all positions
p1 = ggplot(mean.posterior, aes(x = iteration, y = posterior)) +
  geom_line() +
  xlab("Sample") + ylab("Median posterior")
# One position
p2 = ggplot(mean.posterior, aes(x = iteration, y = posone)) +
  geom_line() +
  xlab("Sample") + ylab("Posterior position one")

# ggplot(mean.posterior, aes(x = iteration, y = pos1000)) +
#   geom_line()

# Total genetic distance ~ sample (iteration)
p3 = ggplot(mean.posterior, aes(x = iteration, y = totalgenetic)) +
  geom_line() +
  xlab("Sample") + ylab("Total genetic distance")

p = ggarrange(p1, p2, p3, nrow = 3)
ggsave(paste(wdirpop, "/mcmc/", set, ".", chromosome, ".bpen", bpen, ".jpeg", sep = ""),
       plot = p, device = jpeg(), width = 7, height = 12, dpi = 300)

