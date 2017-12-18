# R Script for p-value calculation

arguments <- commandArgs(trailingOnly = TRUE)
filename <- arguments[1]
output <- arguments[2]

df <- read.table(filename, header = FALSE)

scores <- df$V3
transcript_labels <- df$V6
transcript = df$V1

resultsPval = numeric(length(scores))
for (i in 1:length(scores)) {
	score = as.numeric(unlist(strsplit(as.character(scores[i]), ",")))
	label = as.factor(unlist(strsplit(as.character(transcript_labels[i]),",")))
	resultsPval[i] <- kruskal.test(score, label)[3]
	print(transcript[i])
}

df2 <- data.frame(df, unlist(resultsPval))
write.table(format(df2, digits = 6, scientific = F), file = output, quote = F, sep = "\t", col.names = F, row.names = F)