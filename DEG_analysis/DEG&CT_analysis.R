#Installing DESeq2
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
BiocManager::install("apeglm")
BiocManager::install("genefilter")
install.packages("gplots")
#library load
library(DESeq2)

#DESeq2: Differential expression analysis and Comparative transcriptomics
#Ref2011
#total RNA
load("totalRNAcountmatrix2011.RData")
head(totalRNAcountmatrix2011$counts)
#Creating colData and defining design for DESeq2
totalRNA2011counts=totalRNAcountmatrix2011$counts
#Names simplified
colnames(totalRNA2011counts)
colnames(totalRNA2011counts)<-c("sensitive1","sensitive2","sensitive3","resistant4","resistant5","resistant6")
#Conditions
condition<-c("sensitive","sensitive","sensitive","resistant","resistant","resistant")
#Seq reads types
Seqtype<-c("pairend","pairend","pairend","pairend","pairend","pairend")
#colData
coldata<-matrix(c(condition,Seqtype),nrow=6,ncol = 2)
rownames(coldata)<-colnames(totalRNA2011counts)
colnames(coldata)<-c("condition","type")
coldata=as.data.frame(coldata)
coldata$condition=as.factor(coldata$condition)
coldata$type=as.factor(coldata$type)
coldata
#creating the DESeqDataSet object
ddsTotalRNAref2011<-DESeqDataSetFromMatrix(countData = totalRNA2011counts,
                                           colData = coldata,
                                           design = ~ condition)
#Prefiltering according to tutorial
#Here we perform pre-filtering to keep only rows that have a count of at least 10 for a minimal number of samples.
#The count of 10 is a reasonable choice for bulk RNA-seq. 
#A recommendation for the minimal number of samples is to specify the smallest group size, e.g. here there are 3 treated samples.-cited from tutorial
#also for this case the smallest groupsize is 3 (3 sensitive and 3 resistant)
smallestGroupSize <- 3
keep <- rowSums(counts(ddsTotalRNAref2011) >= 10) >= smallestGroupSize
ddsTotalRNAref2011 <- ddsTotalRNAref2011[keep,]
#Releveling according to control-sensitive
ddsTotalRNAref2011$condition <- relevel(ddsTotalRNAref2011$condition, ref="sensitive")
ddsTotalRNAref2011$condition
#DEG
#Standard analysis
ddsTotalRNAref2011<-DESeq(ddsTotalRNAref2011)
resTotalRNAref2011<-results(ddsTotalRNAref2011)
resTotalRNAref2011
#LFC shrinkage
resLFCTotalRNAref2011 <- lfcShrink(ddsTotalRNAref2011, coef="condition_resistant_vs_sensitive", type="apeglm")
resNormTotalRNAref2011 <- lfcShrink(ddsTotalRNAref2011, coef="condition_resistant_vs_sensitive", type="normal")
#differential expressed genes rank
resLFCTotalRNAref2011
#MA plot
par(mfcol=c(1,2))
xlim <- c(1,1e5); ylim <- c(-3,3)
plotMA(resLFCTotalRNAref2011, main="apeglm")
plotMA(resNormTotalRNAref2011, xlim=xlim, ylim=ylim, main="normal")
#Visualize dispersion estimates:
par(mfcol=c(1,1))
plotDispEsts(ddsTotalRNAref2011)
#histograms of pvalues
hist(resNormTotalRNAref2011$pvalue,main = "Normal",breaks = 12)
hist(resLFCTotalRNAref2011$pvalue,main="apelgm shrinkage",breaks=12)
#rlog
rldTotalRNAref2011<-rlog(ddsTotalRNAref2011)
head(assay(rldTotalRNAref2011))
par( mfrow = c( 1, 2 ) )
plot( log2( 1+counts(ddsTotalRNAref2011, normalized=TRUE)[, 1:2] ), col="#00000020", pch=20, cex=0.3 )
plot( assay(rldTotalRNAref2011)[, 1:2], col="#00000020", pch=20, cex=0.3 )
#heatmap
sampleDistTotalRNAref2011<-dist( t(assay(rldTotalRNAref2011)))
sampleDistTotalRNAref2011
sampleDistMatrix <- as.matrix( sampleDistTotalRNAref2011)
sampleDistMatrix

library( "RColorBrewer" )
colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap.2( sampleDistMatrix, trace="none")
print( plotPCA(rldTotalRNAref2011 ) )

library( "genefilter" )
library("gplots")
par(mfcol=c(1,1),mar=c(1,1,1,1))
topVarGenes <- head( order( rowVars( assay(rldTotalRNAref2011) ), decreasing=TRUE ), 35 )
heatmap.2(assay(rldTotalRNAref2011)[ topVarGenes, ], scale="row",
           trace="none", dendrogram="column",
           col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),cexCol = 0.5)

#Conclusion:
"Seems doubtly, samples being of one condition does not seems to share same characteristics"
"gene cluster profiles also look abnormal"
"there are apparent outliers in the shrinkage MA plots and dispersion plots"
