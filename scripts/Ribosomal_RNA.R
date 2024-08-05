#Percentage of reads mapped to rRNA by featureCount
#FeatureCount is perfectly viable to characterize the transcriptome 
library("ggplot2")
library("ggpubr")
#Mapped to genes
#2011
#Total RNA
setwd("C:/Users/toros/OneDrive - KU Leuven/My thesis/Research/thesis/thesis-results/readmapping/ref2011/totalRNA/postAlignmentQuality/Total-RNA-seq-data--with-2011-reference-genome_multiqc_report_data/")
totalRNA2011=read.table(file="multiqc_rseqc_bam_stat.txt",header = T)
head(totalRNA2011)
#number of reads in primary alignment
totalreadprimary=totalRNA2011$mapq_lt_mapq_cut_non.unique + totalRNA2011$mapq_gte_mapq_cut_unique
totalreadprimary
setwd("C:/Users/toros/OneDrive - KU Leuven/My thesis/Research/thesis/thesis-results/count_matrices/totalRNA")
load("UniqueMapGeneTotalRNAcountmatrix2011reads.RData")
#read maps to genes
fracMapGenesTotalRNA2011=totalRNAcountmatrix2011reads$stat[1,c(5,6,7,2,3,4)]/totalreadprimary
fracMapGenesTotalRNA2011
#Poli RNA
setwd("C:/Users/toros/OneDrive - KU Leuven/My thesis/Research/thesis/thesis-results/readmapping/ref2011/poliRNA/PostAlignmentQuality/poli-A-enriched-RNA-seq-data--with-2011-reference-genome_multiqc_report_data/")
poliRNA2011=read.table(file="multiqc_rseqc_bam_stat.txt",header = T)
head(poliRNA2011)
#number of reads in primary alignment
totalreadprimary=poliRNA2011$mapq_lt_mapq_cut_non.unique + poliRNA2011$mapq_gte_mapq_cut_unique
totalreadprimary
setwd("C:/Users/toros/OneDrive - KU Leuven/My thesis/Research/thesis/thesis-results/count_matrices/poliRNA")
load("UniqueMapGenepoliRNAcountmatrix2011reads.RData")
#read maps to genes
fracMapGenesPoliRNA2011=poliRNAcountmatrix2011reads$stat[1,c(5,6,7,2,3,4)]/totalreadprimary

#2019
#Total RNA
setwd("C:/Users/toros/OneDrive - KU Leuven/My thesis/Research/thesis/thesis-results/readmapping/ref2019/totalRNA/postAlignmentQuality/Total-RNA-seq-data--with-2019-reference-genome_multiqc_report_data/")
totalRNA2019=read.table(file="multiqc_rseqc_bam_stat.txt",header = T)
head(totalRNA2019)
#number of reads in primary alignment
totalreadprimary=totalRNA2019$mapq_lt_mapq_cut_non.unique + totalRNA2019$mapq_gte_mapq_cut_unique
totalreadprimary
setwd("C:/Users/toros/OneDrive - KU Leuven/My thesis/Research/thesis/thesis-results/count_matrices/totalRNA")
load("UniqueMaptotalRNAcountmatrix2019reads.RData")
#read maps to genes
fracMapGenesTotalRNA2019=totalRNAcountmatrix2019reads$stat[1,c(5,6,7,2,3,4)]/totalreadprimary
fracMapGenesTotalRNA2019
#Poli RNA
setwd("C:/Users/toros/OneDrive - KU Leuven/My thesis/Research/thesis/thesis-results/readmapping/ref2019/poliRNA/PostAlignmentQuality/Poli-enrinched-RNA-seq-data--with-2019-reference-genome_multiqc_report_data/")
poliRNA2019=read.table(file="multiqc_rseqc_bam_stat.txt",header = T)
head(poliRNA2019)
#number of reads in primary alignment
totalreadprimary=poliRNA2019$mapq_lt_mapq_cut_non.unique + poliRNA2019$mapq_gte_mapq_cut_unique
totalreadprimary
setwd("C:/Users/toros/OneDrive - KU Leuven/My thesis/Research/thesis/thesis-results/count_matrices/poliRNA")
load("poliRNAcountmatrix2019reads.RData")
#read maps to genes
fracMapGenesPoliRNA2019=poliRNAcountmatrix2019reads$stat[1,c(5,6,7,2,3,4)]/totalreadprimary
# Plots
#MaptoGenes2011
vectortotalRNA2011<-c(fracMapGenesTotalRNA2011$`totalRNA-sensitive1-ref2011Aligned.sortedByCoord.out.bam`,fracMapGenesTotalRNA2011$`totalRNA-sensitive2-ref2011Aligned.sortedByCoord.out.bam`,fracMapGenesTotalRNA2011$`totalRNA-sensitive3-ref2011Aligned.sortedByCoord.out.bam`,fracMapGenesTotalRNA2011$`totalRNA-resistant4-ref2011Aligned.sortedByCoord.out.bam`,fracMapGenesTotalRNA2011$`totalRNA-resistant5-ref2011Aligned.sortedByCoord.out.bam`,fracMapGenesTotalRNA2011$`totalRNA-resistant6-ref2011Aligned.sortedByCoord.out.bam`)
vectortotalRNA2011=as.matrix(vectortotalRNA2011)
vectorname1=c("sens1","sens2","sens3","res4","res5","res6")
rownames(vectortotalRNA2011)=vectorname1
colnames(vectortotalRNA2011)="fraction"
vectortotalRNA2011=as.data.frame(vectortotalRNA2011)
condition=c("senstive","senstive","senstive","resistant","resistant","resistant")
cond.ordered <- factor(rownames(vectortotalRNA2011),levels=c("sens1", "sens2","sens3","res4","res5","res6"))
p1=ggplot(vectortotalRNA2011,aes(x=cond.ordered,y=fraction,fill=condition))+
  theme(axis.text.x = element_text(vjust =0.1 , hjust=0.5,size = 13),axis.title = element_text(size=13),plot.title = element_text(size=15), axis.text.y = element_text(size = 13),axis.title.x = element_text(size = 13),legend.text = element_text(size = 13),legend.title = element_text(size=13))+
  geom_bar(position="dodge",stat="identity")+
  labs(y="fraction",x="sample")+
  ggtitle("Total RNA, reference 2011")
#MaptoGenes2019
vectortotalRNA2019<-c(fracMapGenesTotalRNA2019$`totalRNA-sensitive1-ref2019Aligned.sortedByCoord.out.bam`,fracMapGenesTotalRNA2019$`totalRNA-sensitive2-ref2019Aligned.sortedByCoord.out.bam`,fracMapGenesTotalRNA2019$`totalRNA-sensitive3-ref2019Aligned.sortedByCoord.out.bam`,fracMapGenesTotalRNA2019$`totalRNA-resistant4-ref2019Aligned.sortedByCoord.out.bam`,fracMapGenesTotalRNA2019$`totalRNA-resistant5-ref2019Aligned.sortedByCoord.out.bam`,fracMapGenesTotalRNA2019$`totalRNA-resistant6-ref2019Aligned.sortedByCoord.out.bam`)
vectortotalRNA2019=as.matrix(vectortotalRNA2019)
vectorname1=c("sens1","sens2","sens3","res4","res5","res6")
rownames(vectortotalRNA2019)=vectorname1
colnames(vectortotalRNA2019)="fraction"
vectortotalRNA2019=as.data.frame(vectortotalRNA2019)

condition=c("senstive","senstive","senstive","resistant","resistant","resistant")
cond.ordered <- factor(rownames(vectortotalRNA2019), levels=c("sens1", "sens2","sens3","res4","res5","res6"))
p2=ggplot(vectortotalRNA2019, aes( x=cond.ordered,y=fraction,fill=condition))+
  theme(axis.text.x = element_text(vjust =0.1 , hjust=0.5,size = 13),axis.title = element_text(size=13),plot.title = element_text(size=15), axis.text.y = element_text(size = 13),axis.title.x = element_text(size = 13),legend.text = element_text(size = 13),legend.title = element_text(size=13))+
  geom_bar(position="dodge",stat="identity")+
  labs(y="fraction",x="sample")+
  ggtitle("Total RNA, reference 2019")

figure=ggarrange(p1,p2,ncol=1,nrow=2)
setwd("C:/Users/toros/Desktop/thesis/thesis-R/plots")
png("genes-characterization-totalRNA.png", width = 800, height = 600) 
figure
dev.off()
#Poli RNA 2011
#MaptoGenes2011
vectorpoliRNA2011<-c(fracMapGenesPoliRNA2011$`poliRNA-sensitive1-ref2011Aligned.sortedByCoord.out.bam`,fracMapGenesPoliRNA2011$`poliRNA-sensitive2-ref2011Aligned.sortedByCoord.out.bam`,fracMapGenesPoliRNA2011$`poliRNA-sensitive3-ref2011Aligned.sortedByCoord.out.bam`,fracMapGenesPoliRNA2011$`poliRNA-resistant4-ref2011Aligned.sortedByCoord.out.bam`,fracMapGenesPoliRNA2011$`poliRNA-resistant5-ref2011Aligned.sortedByCoord.out.bam`,fracMapGenesPoliRNA2011$`poliRNA-resistant6-ref2011Aligned.sortedByCoord.out.bam`)
vectorpoliRNA2011=as.matrix(vectorpoliRNA2011)
vectorname1=c("sample1","sample2","sample3","sample4","sample5","sample6")
rownames(vectorpoliRNA2011)=vectorname1
colnames(vectorpoliRNA2011)="fraction"
vectorpoliRNA2011=as.data.frame(vectorpoliRNA2011)
p3=ggplot(vectorpoliRNA2011, aes( x=rownames(vectorpoliRNA2011),y=fraction))+
  theme(axis.text.x = element_text(vjust =0.1 , hjust=0.5,size = 13),axis.title = element_text(size=13),plot.title = element_text(size=15), axis.text.y = element_text(size = 13),axis.title.x = element_text(size = 13),legend.text = element_text(size = 13),legend.title = element_text(size=13))+
  geom_bar(position="dodge",stat="identity")+
  labs(y="fraction",x="sample")+
  ggtitle(" Poli (A) enriched RNA, reference 2011")
#MapGenes2019
vectorpoliRNA2019<-c(fracMapGenesPoliRNA2019$`poliRNA-sensitive1-ref2019Aligned.sortedByCoord.out.bam`,fracMapGenesPoliRNA2019$`poliRNA-sensitive2-ref2019Aligned.sortedByCoord.out.bam`,fracMapGenesPoliRNA2019$`poliRNA-sensitive3-ref2019Aligned.sortedByCoord.out.bam`,fracMapGenesPoliRNA2019$`poliRNA-resistant4-ref2019Aligned.sortedByCoord.out.bam`,fracMapGenesPoliRNA2019$`poliRNA-resistant5-ref2019Aligned.sortedByCoord.out.bam`,fracMapGenesPoliRNA2019$`poliRNA-resistant6-ref2019Aligned.sortedByCoord.out.bam`)
vectorpoliRNA2019=as.matrix(vectorpoliRNA2019)
vectorname1=c("sample1","sample2","sample3","sample4","sample5","sample6")
rownames(vectorpoliRNA2019)=vectorname1
colnames(vectorpoliRNA2019)="fraction"
vectorpoliRNA2019=as.data.frame(vectorpoliRNA2019)
p4=ggplot(vectorpoliRNA2019, aes( x=rownames(vectorpoliRNA2019),y=fraction))+
  theme(axis.text.x = element_text(vjust =0.1 , hjust=0.5,size = 13),axis.title = element_text(size=13),plot.title = element_text(size=15), axis.text.y = element_text(size = 13),axis.title.x = element_text(size = 13),legend.text = element_text(size = 13),legend.title = element_text(size=13))+
  geom_bar(position="dodge",stat="identity")+
  labs(y="fraction",x="sample")+
  ggtitle("poli (A) enriched RNA, reference 2019")

figure=ggarrange(p3,p4,ncol=1,nrow=2)
setwd("C:/Users/toros/Desktop/thesis/thesis-R/plots")
png("genes-characterization-poliRNA.png", width = 800, height = 600) 
figure
dev.off()


#Mapped rRNA
#2011
#Total RNA

#number of reads in primary alignment
totalreadprimary=totalRNA2011$mapq_lt_mapq_cut_non.unique + totalRNA2011$mapq_gte_mapq_cut_unique
totalreadprimary
setwd("C:/Users/toros/OneDrive - KU Leuven/My thesis/Research/thesis/thesis-results/count_matrices/totalRNA")
load("UniqueMaprRNAcountmatrix2011reads.RData")
#read maps to genes
fracMaprRNATotalRNA2011=rRNAcountmatrix2011reads$stat[1,c(5,6,7,2,3,4)]/totalreadprimary
fracMaprRNATotalRNA2011
#Poli RNA
#number of reads in primary alignment
totalreadprimary=poliRNA2011$mapq_lt_mapq_cut_non.unique + poliRNA2011$mapq_gte_mapq_cut_unique
totalreadprimary
setwd("C:/Users/toros/OneDrive - KU Leuven/My thesis/Research/thesis/thesis-results/count_matrices/poliRNA")
load("UniqueMaprRNApoliRNAcountmatrix2011reads.RData")
#read maps to genes
fracMaprRNAPoliRNA2011=rRNAcountmatrix2011poliRNAreads$stat[1,c(5,6,7,2,3,4)]/totalreadprimary
fracMaprRNAPoliRNA2011
#2019
#Total RNA
#number of reads in primary alignment
totalreadprimary=totalRNA2019$mapq_lt_mapq_cut_non.unique + totalRNA2019$mapq_gte_mapq_cut_unique
totalreadprimary
setwd("C:/Users/toros/OneDrive - KU Leuven/My thesis/Research/thesis/thesis-results/count_matrices/totalRNA")
load("UniqueMaprRNAcountmatrix2019.RData")
#read maps to genes
fracMaprRNATotalRNA2019=rRNAcountmatrix2019$stat[1,c(5,6,7,2,3,4)]/totalreadprimary
fracMaprRNATotalRNA2019
#Poli RNA
#number of reads in primary alignment
totalreadprimary=poliRNA2019$mapq_lt_mapq_cut_non.unique + poliRNA2019$mapq_gte_mapq_cut_unique
totalreadprimary
setwd("C:/Users/toros/OneDrive - KU Leuven/My thesis/Research/thesis/thesis-results/count_matrices/poliRNA")
load("rRNApoliRNAcountmatrix2019.RData")
#read maps to genes
fracMaprRNAPoliRNA2019=rRNAcountmatrix2019poliRNA$stat[1,c(5,6,7,2,3,4)]/totalreadprimary
fracMaprRNAPoliRNA2019
#Plots
# Total RNA 
#MaptorRNA 2011
vectorrRNAtotalRNA2011<-c(fracMaprRNATotalRNA2011$`totalRNA-sensitive1-ref2011Aligned.sortedByCoord.out.bam`,fracMaprRNATotalRNA2011$`totalRNA-sensitive2-ref2011Aligned.sortedByCoord.out.bam`,fracMaprRNATotalRNA2011$`totalRNA-sensitive3-ref2011Aligned.sortedByCoord.out.bam`,fracMaprRNATotalRNA2011$`totalRNA-resistant4-ref2011Aligned.sortedByCoord.out.bam`,fracMaprRNATotalRNA2011$`totalRNA-resistant5-ref2011Aligned.sortedByCoord.out.bam`,fracMaprRNATotalRNA2011$`totalRNA-resistant6-ref2011Aligned.sortedByCoord.out.bam`)
vectorrRNAtotalRNA2011=as.matrix(vectorrRNAtotalRNA2011)
vectorname1=c("sens1","sens2","sens3","res4","res5","res6")
rownames(vectorrRNAtotalRNA2011)=vectorname1
colnames(vectorrRNAtotalRNA2011)="fraction"
vectorrRNAtotalRNA2011=as.data.frame(vectorrRNAtotalRNA2011)
condition=c("senstive","senstive","senstive","resistant","resistant","resistant")
cond.ordered <- factor(rownames(vectorrRNAtotalRNA2011), levels=c("sens1", "sens2","sens3","res4","res5","res6"))
p5=ggplot(vectorrRNAtotalRNA2011, aes( x=cond.ordered,y=fraction,fill=condition))+
  theme(axis.text.x = element_text(vjust =0.1 , hjust=0.5,size = 13),axis.title = element_text(size=13),plot.title = element_text(size=15), axis.text.y = element_text(size = 13),axis.title.x = element_text(size = 13),legend.text = element_text(size = 13),legend.title = element_text(size=13))+
  geom_bar(position="dodge",stat="identity")+
  labs(y="fraction",x="sample")+
  ggtitle("Total RNA, reference 2011")
#MaptorRNA 2019
vectorrRNAtotalRNA2019<-c(fracMaprRNATotalRNA2019$`totalRNA-sensitive1-ref2019Aligned.sortedByCoord.out.bam`,fracMaprRNATotalRNA2019$`totalRNA-sensitive2-ref2019Aligned.sortedByCoord.out.bam`,fracMaprRNATotalRNA2019$`totalRNA-sensitive3-ref2019Aligned.sortedByCoord.out.bam`,fracMaprRNATotalRNA2019$`totalRNA-resistant4-ref2019Aligned.sortedByCoord.out.bam`,fracMaprRNATotalRNA2019$`totalRNA-resistant5-ref2019Aligned.sortedByCoord.out.bam`,fracMaprRNATotalRNA2019$`totalRNA-resistant6-ref2019Aligned.sortedByCoord.out.bam`)
vectorrRNAtotalRNA2019=as.matrix(vectorrRNAtotalRNA2019)
vectorname1=c("sens1","sens2","sens3","res4","res5","res6")
rownames(vectorrRNAtotalRNA2019)=vectorname1
colnames(vectorrRNAtotalRNA2019)="fraction"
vectorrRNAtotalRNA2019=as.data.frame(vectorrRNAtotalRNA2019)
condition=c("senstive","senstive","senstive","resistant","resistant","resistant")
cond.ordered <- factor(rownames(vectorrRNAtotalRNA2019), levels=c("sens1", "sens2","sens3","res4","res5","res6"))
p6=ggplot(vectorrRNAtotalRNA2019, aes( x=cond.ordered,y=fraction,fill=condition))+
  theme(axis.text.x = element_text(vjust =0.1 , hjust=0.5,size = 13),axis.title = element_text(size=13),plot.title = element_text(size=15), axis.text.y = element_text(size = 13),axis.title.x = element_text(size = 13),legend.text = element_text(size = 13),legend.title = element_text(size=13))+
  geom_bar(position="dodge",stat="identity")+
  labs(y="fraction",x="sample")+
  ggtitle("Total RNA, reference 2019")

figure=ggarrange(p5,p6,ncol=1,nrow=2)
setwd("C:/Users/toros/Desktop/thesis/thesis-R/plots")
png("rRNA-characterization-totalRNA.png", width = 800, height = 600) 
figure
dev.off()
# Poli RNA 
#MaptorRNA 2011
vectorrRNApoliRNA2011<-c(fracMaprRNAPoliRNA2011$`poliRNA-sensitive1-ref2011Aligned.sortedByCoord.out.bam`,fracMaprRNAPoliRNA2011$`poliRNA-sensitive2-ref2011Aligned.sortedByCoord.out.bam`,fracMaprRNAPoliRNA2011$`poliRNA-sensitive3-ref2011Aligned.sortedByCoord.out.bam`,fracMaprRNAPoliRNA2011$`poliRNA-resistant4-ref2011Aligned.sortedByCoord.out.bam`,fracMaprRNAPoliRNA2011$`poliRNA-resistant5-ref2011Aligned.sortedByCoord.out.bam`,fracMaprRNAPoliRNA2011$`poliRNA-resistant6-ref2011Aligned.sortedByCoord.out.bam`)
vectorrRNApoliRNA2011=as.matrix(vectorrRNApoliRNA2011)
vectorname1=c("sample1","sample2","sample3","sample4","sample5","sample6")
rownames(vectorrRNApoliRNA2011)=vectorname1
colnames(vectorrRNApoliRNA2011)="fraction"
vectorrRNApoliRNA2011=as.data.frame(vectorrRNApoliRNA2011)
p7=ggplot(vectorrRNApoliRNA2011, aes( x=rownames(vectorrRNApoliRNA2011),y=fraction))+
  theme(axis.text.x = element_text(vjust =0.1 , hjust=0.5,size = 13),axis.title = element_text(size=13),plot.title = element_text(size=15), axis.text.y = element_text(size = 13),axis.title.x = element_text(size = 13),legend.text = element_text(size = 13),legend.title = element_text(size=13))+
  geom_bar(position="dodge",stat="identity")+
  labs(y="fraction",x="sample")+
  ggtitle("Poli (A) enriched RNA, reference 2011")
#MaptorRNA 2019
vectorrRNApoliRNA2019<-c(fracMaprRNAPoliRNA2019$`poliRNA-sensitive1-ref2019Aligned.sortedByCoord.out.bam`,fracMaprRNAPoliRNA2019$`poliRNA-sensitive2-ref2019Aligned.sortedByCoord.out.bam`,fracMaprRNAPoliRNA2019$`poliRNA-sensitive3-ref2019Aligned.sortedByCoord.out.bam`,fracMaprRNAPoliRNA2019$`poliRNA-resistant4-ref2019Aligned.sortedByCoord.out.bam`,fracMaprRNAPoliRNA2019$`poliRNA-resistant5-ref2019Aligned.sortedByCoord.out.bam`,fracMaprRNAPoliRNA2019$`poliRNA-resistant6-ref2019Aligned.sortedByCoord.out.bam`)
vectorrRNApoliRNA2019=as.matrix(vectorrRNApoliRNA2019)
vectorname1=c("sample1","sample2","sample3","sample4","sample5","sample6")
rownames(vectorrRNApoliRNA2019)=vectorname1
colnames(vectorrRNApoliRNA2019)="fraction"
vectorrRNApoliRNA2019=as.data.frame(vectorrRNApoliRNA2019)
p8=ggplot(vectorrRNApoliRNA2019, aes( x=rownames(vectorrRNApoliRNA2019),y=fraction))+
  theme(axis.text.x = element_text(vjust =0.1 , hjust=0.5,size = 13),axis.title = element_text(size=13),plot.title = element_text(size=15), axis.text.y = element_text(size = 13),axis.title.x = element_text(size = 13),legend.text = element_text(size = 13),legend.title = element_text(size=13))+
  geom_bar(position="dodge",stat="identity")+
  labs(y="fraction",x="sample")+
  ggtitle("poli (A) enriched RNA, reference 2019")
figure=ggarrange(p7,p8,ncol=1,nrow=2)
setwd("C:/Users/toros/Desktop/thesis/thesis-R/plots")
png("rRNA-characterizarion-poliRNA.png", width = 800, height = 600) 
figure
dev.off()

