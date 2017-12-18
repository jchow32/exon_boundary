# Rscript for adjusted p-value calculation per transcript

# Allowing input file from command shell
arguments <- commandArgs(trailingOnly = TRUE)
filename <- arguments[1]
output <- arguments[2]

# Reading files
df <- read.table(filename, header = FALSE)

# Transcripts in sample
transcripts <- unique(as.vector(df$V1))

# FDR calculation for each transcript
resultsFdr <- data.frame()
for(i in 1:length(transcripts)){
  data <- df[df$V1 == transcripts[i],]
  dataFdr  <- p.adjust(data$V5, 'fdr')
  dataFdr2 <- data.frame(data, dataFdr)
  resultsFdr <- rbind(resultsFdr,dataFdr2)
}

# Writing results
write.table(format(resultsFdr, digits = 6, scientific = F), file = output, quote = F, sep = "\t", col.names = F, row.names = F)