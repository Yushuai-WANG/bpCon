#############################   Rscript   ################################
#!/bin/R
suppressPackageStartupMessages(library(DESeq2))
args <- commandArgs(T)
data=read.table(args[1],stringsAsFactors=F,sep = "\t", header = F, row.names = 1)
sampleNames <- c("Input1_1","Input1_2","IP1_1","IP1_2")
names(data) <- sampleNames
data1 <- round(as.matrix(data[,1:4]))
condition1 <- factor(c(rep("Input1",2),rep("IP1",2)),levels = c("Input1","IP1"))
coldata1 <- data.frame(row.names = colnames(data1), condition1)
dds1 <- DESeqDataSetFromMatrix(countData=data1, colData=coldata1, design=~condition1)

# DEGs
f1 = as.numeric(args[4])
f2 = as.numeric(args[5])
f3 = as.numeric(args[6])
f4 = as.numeric(args[7])
sizeFactors(dds1) <- c(f1*2/1000000,f2*2/1000000,f3*2/1000000,f4*2/1000000)
dds1 <- estimateDispersionsGeneEst(dds1)
dispersions(dds1) <- mcols(dds1)$dispGeneEst
dds1 <- nbinomWaldTest(dds1)
res1 <- results(dds1, lfcThreshold = 0, altHypothesis = "greater")
res1 <- res1[order(res1$padj),]

write.table(res1, file = args[2],sep="\t",quote=F,col.names = F,row.names=T,eol="\n")

normalized_counts1 <- counts(dds1, normalized=TRUE)
normalized_counts1 <- as.data.frame(normalized_counts1)

write.table(normalized_counts1, file = args[3],sep="\t",quote=F,col.names = F,row.names=T,eol="\n")

############################ E    N    D ########################################

