#!/usr/bin/env Rscript
# Selection of best hyperparameters in Pyrho
# based on maximizing LogL2

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)!=2) {
  stop("<dataset> and <chromosome> must be arguments.", call.=FALSE)
}

dataset=args[1]
chrom=args[2]
wd=paste("data/", dataset, sep = "")

df = read.table(paste(wd, "/pyrho/", dataset, ".hyperparam.", chrom, sep = ""), header = T)

print(df)

# Keep lines maximizing LogL2
df.maxLogL2 = df[which.max(df$Log_L2),]
# If more than one, keep lowest windowsize
if (nrow(df.maxLogL2) > 1) {
    df.maxLogL2 = df[which.min(df.maxLogL2$Window_Size),]
}
# If more than one, keep lowest block penalty
if (nrow(df.maxLogL2) > 1) {
    df.maxLogL2 = df[which.min(df.maxLogL2$Block_Penalty),]
}

# Check we have finally a single pair of hyperparameters
if (nrow(df.maxLogL2) == 1) {
    windowsize = df.maxLogL2$Window_Size
    bpen = df.maxLogL2$Block_Penalty
} else {
    stop(paste("Number of results: ", nrow(df.maxLogL2), "\nOnly one pair of hyperparameters <windowsize> and <bpen> must be retain. Check your hyperparameters file."), sep = "")
}

# Save in two files
write.table(windowsize, file = paste(wd, "/pyrho/", dataset, ".windowsize.", chrom, sep = ""),
            quote = F, col.names = F, row.names = F, sep = "\t")
write.table(bpen, file = paste(wd, "/pyrho/", dataset, ".bpen.", chrom, sep = ""),
            quote = F, col.names = F, row.names = F, sep = "\t")




