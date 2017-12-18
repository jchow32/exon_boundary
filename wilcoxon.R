# R Script for p-value calculation following Kruskal wallis multiple group comparison

arguments <- commandArgs(trailingOnly = TRUE)
filename <- arguments[1]
output <- arguments[2]

df <- read.table(filename, header = FALSE)

transcript = df$V1
gene = df$V2
score1 = df$V3
score2 = df$V4

resultsPval = numeric(length(score1))
for (i in 1:length(score1)) {
	a = as.numeric(unlist(strsplit(as.character(score1[i]), ",")))
	b = as.numeric(unlist(strsplit(as.character(score2[i]),",")))
	resultsPval[i] <- wilcox.test(a, b)[3]
}

df2 <- data.frame(df, unlist(resultsPval))
write.table(format(df2, digits = 6, scientific = F), file = output, quote = F, sep = "\t", col.names = F, row.names = F)