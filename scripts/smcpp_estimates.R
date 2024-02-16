#!/usr/bin/env Rscript
# Transform the SMC++ plot in vectors of Ne and times
args = commandArgs(trailingOnly=TRUE)

# test if there is one argument: if not, return an error
if (length(args) != 4) {
  stop("<input file>, <output dir>, <mu> and <theta> must be passed in argument", call.=FALSE)
}

input = args[1]
output = args[2]
mu = as.numeric(args[3])
theta = as.numeric(args[4])

cat("mu =", mu, "\n")
cat("theta =", theta, "\n")

df = read.table(input, header = T, sep = ",")
df$x = as.numeric(df$x)
df$y = as.numeric(df$y)

df

cat("Get piecewise coalescent-scaled Ne\n")
n_ref = theta / (4. * mu)

y = unique(df$y)
Ne = paste(as.character(y/n_ref), collapse = ",")
Ne
write.table(Ne, paste0(output, "Ne.txt"), col.names = F,
            row.names = F, quote = F)

cat("Get piecewise coalescent-scaled times\n")
t = unlist(lapply(y, function(x) {min(df$x[which(df$y == as.numeric(x))])}))
t
t = t[-1]
t = t/(2*n_ref)
t = paste(as.character(t), collapse = ",")
t
write.table(t, paste0(output, "times.txt"), col.names = F,
            row.names = F, quote = F)
