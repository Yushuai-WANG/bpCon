#############################   Rscript   ################################
#!/bin/R
suppressPackageStartupMessages(library(DESeq2))
args <- commandArgs(T)
data=read.table(args[1],stringsAsFactors=F,sep = "\t", header = F, row.names = 1)
sampleNames <- c("Input1_1","Input1_2","IP1_1","IP1_2","Input2_1","Input2_2","IP2_1","IP2_2")
names(data) <- sampleNames
data1 <- round(as.matrix(data[,1:4]))
data2 <- round(as.matrix(data[,5:8]))
condition1 <- factor(c(rep("Input1",2),rep("IP1",2)),levels = c("Input1","IP1"))
condition2 <- factor(c(rep("Input2",2),rep("IP2",2)),levels = c("Input2","IP2"))
coldata1 <- data.frame(row.names = colnames(data1), condition1)
coldata2 <- data.frame(row.names = colnames(data2), condition2)
dds1 <- DESeqDataSetFromMatrix(countData=data1, colData=coldata1, design=~condition1)
dds2 <- DESeqDataSetFromMatrix(countData=data2, colData=coldata2, design=~condition2)

# DEGs
f1 = as.numeric(args[6])
f2 = as.numeric(args[7])
f3 = as.numeric(args[8])
f4 = as.numeric(args[9])

f5 = as.numeric(args[10])
f6 = as.numeric(args[11])
f7 = as.numeric(args[12])
f8 = as.numeric(args[13])

sizeFactors(dds1) <- c(f1*2/1000000,f2*2/1000000,f3*2/1000000,f4*2/1000000)
sizeFactors(dds2) <- c(f5*2/1000000,f6*2/1000000,f7*2/1000000,f8*2/1000000)
dds1 <- estimateDispersionsGeneEst(dds1)
dds2 <- estimateDispersionsGeneEst(dds2)
dispersions(dds1) <- mcols(dds1)$dispGeneEst
dispersions(dds2) <- mcols(dds2)$dispGeneEst
dds1 <- nbinomWaldTest(dds1)
dds2 <- nbinomWaldTest(dds2)
res1 <- results(dds1, lfcThreshold = 0, altHypothesis = "greater")
res2 <- results(dds2, lfcThreshold = 0, altHypothesis = "greater")
res1 <- res1[order(res1$padj),]
res2 <- res2[order(res2$padj),]

write.table(res1, file = args[2],sep="\t",quote=F,col.names = F,row.names=T,eol="\n")
write.table(res2, file = args[3],sep="\t",quote=F,col.names = F,row.names=T,eol="\n")

normalized_counts1 <- counts(dds1, normalized=TRUE)
normalized_counts2 <- counts(dds2, normalized=TRUE)
normalized_counts1 <- as.data.frame(normalized_counts1)
normalized_counts2 <- as.data.frame(normalized_counts2)

write.table(normalized_counts1, file = args[4],sep="\t",quote=F,col.names = F,row.names=T,eol="\n")
write.table(normalized_counts2, file = args[5],sep="\t",quote=F,col.names = F,row.names=T,eol="\n")

############################ E    N    D ########################################

