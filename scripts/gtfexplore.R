library(reader)
library(ggplot2)
library(dplyr)
#The gtf distribution  graph
gtf2011summarylocation=find.file("gtf_features_summary2011.txt", dir = "C:/Users/toros/OneDrive - KU Leuven/My thesis/Research/thesis/thesis-results/readmapping/gtf-files-used", dirs = NULL)
gtf2019summarylocation=find.file("gtf_features_summary2019.txt", dir = "C:/Users/toros/OneDrive - KU Leuven/My thesis/Research/thesis/thesis-results/readmapping/gtf-files-used", dirs = NULL)
gtf2011summary=read.table(file=gtf2011summarylocation,sep=' ')
names(gtf2011summary)<-c("features_name","feature_counts")
gtf2019summary=read.table(file=gtf2019summarylocation,sep=' ')
names(gtf2019summary)<-c("features_name","feature_counts")
names(gtf2011summary)
gtf2011summary
row12=c("snoRNA",0)
row13=c("snRNA",0)
row14=c("three_prime_UTR",0)
row15=c("ncRNA",0)
gtf2011summary=rbind(gtf2011summary,row12)
gtf2011summary=rbind(gtf2011summary,row13)
gtf2011summary$genome=as.factor(2011)
gtf2019summary=rbind(gtf2019summary,row14)
gtf2019summary=rbind(gtf2019summary,row15)
gtf2019summary$genome=as.factor(2019)

gtfsummary=rbind(gtf2011summary,gtf2019summary)
gtfsummary$feature_counts=as.integer(gtfsummary$feature_counts)
setwd("C:/Users/toros/Desktop/thesis/thesis-R/plots")
png("gff_file_annotation_summary.png", width = 800, height = 600) 
ggplot(gtfsummary, aes(group=genome,x=feature_counts, y=features_name,fill=genome))+
  theme(axis.text.x = element_text(vjust =0.3 , hjust=1))+
  geom_bar(position="dodge",stat="identity")+
  geom_text(aes(label=feature_counts),vjust = 0.3,hjust=0,size=3,position = position_dodge(0.9))+
  labs(y="feature count",x="feature name",size=3)+
  scale_x_continuous(limits = c(0, 15000))+
  ggtitle("gff file annotation summary")
dev.off()

