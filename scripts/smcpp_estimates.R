#!/usr/bin/env Rscript
# Transform the SMC++ plot in vectors of Ne and times
args = commandArgs(trailingOnly=TRUE)

# test if there is one argument: if not, return an error
if (length(args) != 2) {
  stop("<input file> and <output dir> must be passed in argument", call.=FALSE)
}

input = args[1]
output = args[2]

df = read.table(input, header = T, sep = ",")

Ne = as.character(unique(df$y))
Ne = paste(Ne, collapse = ",")
Ne
write.table(Ne, paste0(output, "Ne.txt"), col.names = F,
            row.names = F, quote = F)

t = unlist(lapply(Ne, function(x) {min(df$x[which(df$y == x)])}))
t = t[-1]
t = paste(t, collapse = ",")
t
write.table(t, paste0(output, "times.txt"), col.names = F,
            row.names = F, quote = F)