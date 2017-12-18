# R Script for p-value calculation

arguments <- commandArgs(trailingOnly = TRUE)
filename <- arguments[1]
output <- arguments[2]

df <- read.table(filename, header = FALSE)

samSuccess <- df$V3 # missense
bgSuccess <- df$V5 # total missense
bgFailure <- df$V6 # total syn
samSize <- df$V3 + df$V4 # missense + syn

# Hypergeometrix test
resultsPval <- phyper(samSuccess,bgSuccess,bgFailure,samSize, lower.tail=TRUE)
df2 <- data.frame(df, resultsPval)

write.table(format(df2, digits = 6, scientific = F), file = output, quote = F, sep = "\t", col.names = F, row.names = F)

# Fisher option
#df <- data.frame(data[6], data[7], data[8]-data[6], data[9]-data[7])
#pvalues <- apply(df,1,function(x) fisher.test(matrix(x, nrow = 2), alt = "less")$p.value)
#df2 <- data.frame(data, df, pvalues)
#write.table(df2, file = "exacdb.nonpsych_synnorm_pvalues", quote=F, sep="\t", col.names = NA)