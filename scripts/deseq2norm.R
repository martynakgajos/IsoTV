# Normalize isoform expression across all conditions

library("DESeq2")

samples <- snakemake@config["samples"]$samples
data <- read.csv(file=snakemake@input[["counts"]], header=TRUE, sep=",")
countData <- data[colnames(data)[2:(length(samples) + 1)]]
countData <- as.matrix(countData)

name <- colnames(data)[2:(length(samples) + 1)]
condition <- sapply(strsplit(name,"_"), `[`, 1)
colData <- cbind(name,condition)

dds <- DESeqDataSetFromMatrix(countData,colData,design=~condition)
dds <- estimateSizeFactors(dds)
Pcounts <- counts(dds, normalized=TRUE)

tpm <- function(x) {
  x * 1000000/sum(x)
}

Pcounts <- apply(Pcounts,2,tpm)

data[colnames(data)[2:(length(samples) + 1)]] <- Pcounts

data <- data[order(data$transcript_id),]
write.table(data,file=snakemake@output[["norm_counts"]],row.names = FALSE,quote=FALSE,sep=",")
