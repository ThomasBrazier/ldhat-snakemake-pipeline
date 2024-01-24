#!/usr/bin/env Rscript
# Detect windows with low SNP density to mask
args = commandArgs(trailingOnly=TRUE)

# test if there is one argument: if not, return an error
if (length(args) != 3) {
  stop("<input file>, <bin size> and <min number of SNPs> must be passed in argument (in this order)", call.=FALSE)
}

input = args[1]
bin = as.numeric(args[2])
min = as.numeric(args[3])

snpdens = read.table(paste0(input, ".snpden"),
                 header = T)

snpden = as.numeric(min)/as.numeric(bin) * 1000

mask = snpdens[which(snpdens$VARIANTS.KB < snpden),]
mask$chromStart = mask$BIN_START - 1
mask$chromEnd = mask$BIN_START + bin - 1

colnames(mask)[1] = "chrom"

mask[,c("chrom", "chromStart", "chromEnd")]

write.table(mask[,c("chrom", "chromStart", "chromEnd")],
            paste0(input, ".bed"), sep = "\t",
            col.names = F, row.names = F, quote = F)


