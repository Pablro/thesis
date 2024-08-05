#Intsalling Feature Count
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.18")
BiocManager::install("Rsubread")
install.packages("reader")
#libraries
library(Rsubread)
library(reader)
#Ref2011
#totalRNA
#Verication of existence files:
find.file("totalRNA-sensitive1-ref2011Aligned.sortedByCoord.out.bam", dir = "C:/Users/toros/OneDrive - KU Leuven/My thesis/Research/thesis/thesis-results/readmapping/ref2011/totalRNA", dirs = NULL)
find.file("totalRNA-sensitive2-ref2011Aligned.sortedByCoord.out.bam", dir = "C:/Users/toros/OneDrive - KU Leuven/My thesis/Research/thesis/thesis-results/readmapping/ref2011/totalRNA", dirs = NULL)
find.file("totalRNA-sensitive3-ref2011Aligned.sortedByCoord.out.bam", dir = "C:/Users/toros/OneDrive - KU Leuven/My thesis/Research/thesis/thesis-results/readmapping/ref2011/totalRNA", dirs = NULL)
find.file("totalRNA-resistant4-ref2011Aligned.sortedByCoord.out.bam", dir = "C:/Users/toros/OneDrive - KU Leuven/My thesis/Research/thesis/thesis-results/readmapping/ref2011/totalRNA", dirs = NULL)
find.file("totalRNA-resistant5-ref2011Aligned.sortedByCoord.out.bam", dir = "C:/Users/toros/OneDrive - KU Leuven/My thesis/Research/thesis/thesis-results/readmapping/ref2011/totalRNA", dirs = NULL)
find.file("totalRNA-resistant6-ref2011Aligned.sortedByCoord.out.bam", dir = "C:/Users/toros/OneDrive - KU Leuven/My thesis/Research/thesis/thesis-results/readmapping/ref2011/totalRNA", dirs = NULL)
find.file("genome_2011.gff",dir="C:/Users/toros/OneDrive - KU Leuven/My thesis/Research/thesis/thesis-results/readmapping/gtf-files-used",dirs=NULL)
#files to variables:
totalsensitive1="C:/Users/toros/OneDrive - KU Leuven/My thesis/Research/thesis/thesis-results/readmapping/ref2011/totalRNA/totalRNA-sensitive1-ref2011Aligned.sortedByCoord.out.bam"
totalsensitive2="C:/Users/toros/OneDrive - KU Leuven/My thesis/Research/thesis/thesis-results/readmapping/ref2011/totalRNA/totalRNA-sensitive2-ref2011Aligned.sortedByCoord.out.bam"
totalsensitive3="C:/Users/toros/OneDrive - KU Leuven/My thesis/Research/thesis/thesis-results/readmapping/ref2011/totalRNA/totalRNA-sensitive3-ref2011Aligned.sortedByCoord.out.bam"
totalresistant4="C:/Users/toros/OneDrive - KU Leuven/My thesis/Research/thesis/thesis-results/readmapping/ref2011/totalRNA/totalRNA-resistant4-ref2011Aligned.sortedByCoord.out.bam"
totalresistant5="C:/Users/toros/OneDrive - KU Leuven/My thesis/Research/thesis/thesis-results/readmapping/ref2011/totalRNA/totalRNA-resistant5-ref2011Aligned.sortedByCoord.out.bam"
totalresistant6="C:/Users/toros/OneDrive - KU Leuven/My thesis/Research/thesis/thesis-results/readmapping/ref2011/totalRNA/totalRNA-resistant6-ref2011Aligned.sortedByCoord.out.bam"
#Vector of files: (Bam order by coordinates:)
totalRNAref2011BAMfiles<-c(totalsensitive1,totalsensitive2,totalsensitive3,totalresistant4,totalresistant5,totalresistant6)
#annotation file:
gtf2011<-"C:/Users/toros/OneDrive - KU Leuven/My thesis/Research/thesis/thesis-results/readmapping/gtf-files-used/genome_2011.gff"
#Counting reads in fragments:
#Parameters:
#Not overlapping
#unique reads-not multimappers
#both reads mapped
#Allow metafeature if exon=genes
#MQS=255-mapq score-default used by star
#exon
totalRNAcountmatrix2011=featureCounts(totalRNAref2011BAMfiles,annot.ext = gtf2011,isGTFAnnotationFile = T,GTF.featureType = "exon",isPairedEnd = T,countMultiMappingReads = F,minMQS = 255,nthreads =4,useMetaFeatures = T,primaryOnly = T,requireBothEndsMapped = T )
load(file="UniqueMapGenesTotalRNAcountmatrix2011.RData")
#Inspection
head(totalRNAcountmatrix2011$counts)
genes2011totalRNAstats=totalRNAcountmatrix2011$stat
View(genes2011totalRNAstats)
head(totalRNAcountmatrix2011$annotation)
save(totalRNAcountmatrix2011,file ="UniqueMapGenesTotalRNAcountmatrix2011.RData" )
#counting in reads
#Counting reads in fragments:
#Parameters:
#Not overlapping
#unique reads-not multimappers
#both reads mapped
#Allow metafeature if exon=genes
#MQS=255-mapq score-default used by star
#countReadPairs = F this will now count reads rather than fragments. Results will now be comparable to stats of post alignment quality tools
totalRNAcountmatrix2011reads=featureCounts(totalRNAref2011BAMfiles,annot.ext = gtf2011,isGTFAnnotationFile = T,GTF.featureType = "exon",isPairedEnd = T,countReadPairs = F,countMultiMappingReads = F,minMQS = 255,nthreads =4,useMetaFeatures = T,primaryOnly = T,requireBothEndsMapped = T )
load(file="UniqueMapGeneTotalRNAcountmatrix2011reads.RData")
#Inspection
head(totalRNAcountmatrix2011reads$counts)
genes2011totalRNAreadstats=totalRNAcountmatrix2011reads$stat
View(genes2011totalRNAreadstats)
write.table(genes2011totalRNAreadstats,file="UniqueMapGenesTotalRNAcountmatrix2011readstat.txt",quote = F,row.names = F,col.names = T)
head(totalRNAcountmatrix2011reads$annotation)
save(totalRNAcountmatrix2011reads,file ="UniqueMapGeneTotalRNAcountmatrix2011reads.RData" )
#ribosomal rna
#Parameters:
#Not overlapping
#unique reads-not multimappers
#both reads mapped
#MQS=255-mapq score-default used by star
#countReadPairs = F this will now count reads rather than fragments. Results will now be comparable to stats of post alignment quality tools
rRNAcountmatrix2011reads= featureCounts(totalRNAref2011BAMfiles,annot.ext = gtf2011,isGTFAnnotationFile = T,GTF.attrType = "ID",GTF.featureType = "rRNA",isPairedEnd = T,countReadPairs = F,countMultiMappingReads = F,minMQS = 255,nthreads =4,useMetaFeatures = F, primaryOnly = T,requireBothEndsMapped = T )
save(rRNAcountmatrix2011reads,file = "UniqueMaprRNAcountmatrix2011reads.RData")
load(file="UniqueMaprRNAcountmatrix2011reads.RData")
rRNA2011totalRNAreadstats=rRNAcountmatrix2011reads$stat
sum(rRNA2011totalRNAreadstats[,2])
write.table(rRNA2011totalRNAreadstats,file="UniqueMaprRNAcountmatrix2011totalRNAreadstat.txt",quote = F,row.names = F,col.names = T)
View(rRNAcountmatrix2011reads$stat)
#poliRNA
find.file("poliRNA-sensitive1-ref2011Aligned.sortedByCoord.out.bam", dir = "C:/Users/toros/OneDrive - KU Leuven/My thesis/Research/thesis/thesis-results/readmapping/ref2011/poliRNA", dirs = NULL)
find.file("poliRNA-sensitive2-ref2011Aligned.sortedByCoord.out.bam", dir = "C:/Users/toros/OneDrive - KU Leuven/My thesis/Research/thesis/thesis-results/readmapping/ref2011/poliRNA", dirs = NULL)
find.file("poliRNA-sensitive3-ref2011Aligned.sortedByCoord.out.bam", dir = "C:/Users/toros/OneDrive - KU Leuven/My thesis/Research/thesis/thesis-results/readmapping/ref2011/poliRNA", dirs = NULL)
find.file("poliRNA-resistant4-ref2011Aligned.sortedByCoord.out.bam", dir = "C:/Users/toros/OneDrive - KU Leuven/My thesis/Research/thesis/thesis-results/readmapping/ref2011/poliRNA", dirs = NULL)
find.file("poliRNA-resistant5-ref2011Aligned.sortedByCoord.out.bam", dir = "C:/Users/toros/OneDrive - KU Leuven/My thesis/Research/thesis/thesis-results/readmapping/ref2011/poliRNA", dirs = NULL)
find.file("poliRNA-resistant6-ref2011Aligned.sortedByCoord.out.bam", dir = "C:/Users/toros/OneDrive - KU Leuven/My thesis/Research/thesis/thesis-results/readmapping/ref2011/poliRNA", dirs = NULL)
find.file("genome_2011.gff",dir="C:/Users/toros/OneDrive - KU Leuven/My thesis/Research/thesis/thesis-results/readmapping/gtf-files-used",dirs=NULL)
#files to variables
#files to variables:
polisensitive1="C:/Users/toros/OneDrive - KU Leuven/My thesis/Research/thesis/thesis-results/readmapping/ref2011/poliRNA/poliRNA-sensitive1-ref2011Aligned.sortedByCoord.out.bam"
polisensitive2="C:/Users/toros/OneDrive - KU Leuven/My thesis/Research/thesis/thesis-results/readmapping/ref2011/poliRNA/poliRNA-sensitive2-ref2011Aligned.sortedByCoord.out.bam"
polisensitive3="C:/Users/toros/OneDrive - KU Leuven/My thesis/Research/thesis/thesis-results/readmapping/ref2011/poliRNA/poliRNA-sensitive3-ref2011Aligned.sortedByCoord.out.bam"
poliresistant4="C:/Users/toros/OneDrive - KU Leuven/My thesis/Research/thesis/thesis-results/readmapping/ref2011/poliRNA/poliRNA-resistant4-ref2011Aligned.sortedByCoord.out.bam"
poliresistant5="C:/Users/toros/OneDrive - KU Leuven/My thesis/Research/thesis/thesis-results/readmapping/ref2011/poliRNA/poliRNA-resistant5-ref2011Aligned.sortedByCoord.out.bam"
poliresistant6="C:/Users/toros/OneDrive - KU Leuven/My thesis/Research/thesis/thesis-results/readmapping/ref2011/poliRNA/poliRNA-resistant6-ref2011Aligned.sortedByCoord.out.bam"
#Vector of files: (Bam order by coordinates:)
poliRNAref2011BAMfiles<-c(polisensitive1,polisensitive2,polisensitive3,poliresistant4,poliresistant5,poliresistant6)
#annotation file:
gtf2011<-"C:/Users/toros/OneDrive - KU Leuven/My thesis/Research/thesis/thesis-results/readmapping/gtf-files-used/genome_2011.gff"
#Counting fragments:
#Parameters:
#Not overlapping
#unique reads-not multimappers
#both reads mapped
#Allow metafeature if exon=genes
#MQS=255-mapq score-default used by star
poliRNAcountmatrix2011=featureCounts(poliRNAref2011BAMfiles,annot.ext = gtf2011,isGTFAnnotationFile = T,GTF.featureType = "exon",isPairedEnd = T,countMultiMappingReads = F,minMQS = 255,nthreads =4,useMetaFeatures = T,primaryOnly = T,requireBothEndsMapped = T)
#Inspection
load(file="poliRNAcountmatrix2011.RData")
head(poliRNAcountmatrix2011$counts)
View(poliRNAcountmatrix2011$counts)
head(poliRNAcountmatrix2011$annotation)
View(poliRNAcountmatrix2011$stat)
save(poliRNAcountmatrix2011,file ="poliRNAcountmatrix2011.RData" )
Gene2011PoliRNAstats=poliRNAcountmatrix2011$stat
View(Gene2011PoliRNAstats)
sum(Gene2011PoliRNAstats[,2])
#counting reads
#Parameters:
#Not overlapping
#unique reads-not multimappers
#both reads mapped
#Allow metafeature if exon=genes
#MQS=255-mapq score-default used by star
#countReadPairs = F this will now count reads rather than fragments. Results will now be comparable to stats of post alignment quality tools
poliRNAcountmatrix2011reads=featureCounts(poliRNAref2011BAMfiles,annot.ext = gtf2011,isGTFAnnotationFile = T,GTF.featureType = "exon",isPairedEnd = T,countReadPairs = F,countMultiMappingReads = F,minMQS = 255,nthreads =4,useMetaFeatures = T,primaryOnly = T,requireBothEndsMapped = T)
#Inspection
load(file="UniqueMapGenepoliRNAcountmatrix2011reads.RData")
head(poliRNAcountmatrix2011reads$counts)
genes2011poliRNAreadstats=poliRNAcountmatrix2011reads$stat
View(genes2011poliRNAreadstats)
sum(genes2011poliRNAreadstats[,2])
View(poliRNAcountmatrix2011reads$counts)
write.table(genes2011poliRNAreadstats,file="UniqueMapGenespoliRNAcountmatrix2011readstat.txt",quote = F,row.names = F,col.names = T)
View(poliRNAcountmatrix2011reads$annotation)
sum(poliRNAcountmatrix2011reads$stat[2])
save(poliRNAcountmatrix2011reads,file ="UniqueMapGenepoliRNAcountmatrix2011reads.RData" )
#ribosomal rna
#reads
rRNAcountmatrix2011poliRNAreads= featureCounts(poliRNAref2011BAMfiles,annot.ext = gtf2011,isGTFAnnotationFile = T,GTF.attrType = "ID",GTF.featureType = "rRNA",isPairedEnd = T,countReadPairs = F,countMultiMappingReads = F,minMQS = 255,nthreads =4,useMetaFeatures = T,primaryOnly = T,requireBothEndsMapped = T )
#Inspection
head(rRNAcountmatrix2011poliRNAreads$counts)
rRNA2011poliRNAreadstats=rRNAcountmatrix2011poliRNAreads$stat
View(rRNA2011poliRNAreadstats)
sum(rRNA2011poliRNAreadstats[,2])
write.table(rRNA2011poliRNAreadstats,file="UniquerRNAcountmatrix2011poliRNAreadstat.txt",quote = F,row.names = F,col.names = T)
save(rRNAcountmatrix2011poliRNAreads,file ="UniqueMaprRNApoliRNAcountmatrix2011reads.RData" )
########################################################################################################
########################################################################################################
#Ref2019
#totalRNA
#Verication of existence files:
find.file("totalRNA-sensitive1-ref2019Aligned.sortedByCoord.out.bam", dir = "C:/Users/toros/OneDrive - KU Leuven/My thesis/Research/thesis/thesis-results/readmapping/ref2019/totalRNA", dirs = NULL)
find.file("totalRNA-sensitive2-ref2019Aligned.sortedByCoord.out.bam", dir = "C:/Users/toros/OneDrive - KU Leuven/My thesis/Research/thesis/thesis-results/readmapping/ref2019/totalRNA", dirs = NULL)
find.file("totalRNA-sensitive3-ref2019Aligned.sortedByCoord.out.bam", dir = "C:/Users/toros/OneDrive - KU Leuven/My thesis/Research/thesis/thesis-results/readmapping/ref2019/totalRNA", dirs = NULL)
find.file("totalRNA-resistant4-ref2019Aligned.sortedByCoord.out.bam", dir = "C:/Users/toros/OneDrive - KU Leuven/My thesis/Research/thesis/thesis-results/readmapping/ref2019/totalRNA", dirs = NULL)
find.file("totalRNA-resistant5-ref2019Aligned.sortedByCoord.out.bam", dir = "C:/Users/toros/OneDrive - KU Leuven/My thesis/Research/thesis/thesis-results/readmapping/ref2019/totalRNA", dirs = NULL)
find.file("totalRNA-resistant6-ref2019Aligned.sortedByCoord.out.bam", dir = "C:/Users/toros/OneDrive - KU Leuven/My thesis/Research/thesis/thesis-results/readmapping/ref2019/totalRNA", dirs = NULL)
find.file("genome_2011.gff",dir="C:/Users/toros/OneDrive - KU Leuven/My thesis/Research/thesis/thesis-results/readmapping/gtf-files-used",dirs=NULL)
#files to variables:
totalsensitive1="C:/Users/toros/OneDrive - KU Leuven/My thesis/Research/thesis/thesis-results/readmapping/ref2019/totalRNA/totalRNA-sensitive1-ref2019Aligned.sortedByCoord.out.bam"
totalsensitive2="C:/Users/toros/OneDrive - KU Leuven/My thesis/Research/thesis/thesis-results/readmapping/ref2019/totalRNA/totalRNA-sensitive2-ref2019Aligned.sortedByCoord.out.bam"
totalsensitive3="C:/Users/toros/OneDrive - KU Leuven/My thesis/Research/thesis/thesis-results/readmapping/ref2019/totalRNA/totalRNA-sensitive3-ref2019Aligned.sortedByCoord.out.bam"
totalresistant4="C:/Users/toros/OneDrive - KU Leuven/My thesis/Research/thesis/thesis-results/readmapping/ref2019/totalRNA/totalRNA-resistant4-ref2019Aligned.sortedByCoord.out.bam"
totalresistant5="C:/Users/toros/OneDrive - KU Leuven/My thesis/Research/thesis/thesis-results/readmapping/ref2019/totalRNA/totalRNA-resistant5-ref2019Aligned.sortedByCoord.out.bam"
totalresistant6="C:/Users/toros/OneDrive - KU Leuven/My thesis/Research/thesis/thesis-results/readmapping/ref2019/totalRNA/totalRNA-resistant6-ref2019Aligned.sortedByCoord.out.bam"
#Vector of files: (Bam order by coordinates:)
totalRNAref2019BAMfiles<-c(totalsensitive1,totalsensitive2,totalsensitive3,totalresistant4,totalresistant5,totalresistant6)
#annotation file:
gtf2019<-"C:/Users/toros/OneDrive - KU Leuven/My thesis/Research/thesis/thesis-results/readmapping/gtf-files-used/genome_2019.gff"
#Parameters:
#Parameters:
#Not overlapping
#unique reads-not multimappers
#both reads mapped
#Allow metafeature if exon=genes
#MQS=255-mapq score-default used by star
totalRNAcountmatrix2019=featureCounts(totalRNAref2019BAMfiles,annot.ext = gtf2019,isGTFAnnotationFile = T,GTF.featureType = "exon",countMultiMappingReads = F,isPairedEnd = T,nthreads =4,useMetaFeatures = T,primaryOnly = T,requireBothEndsMapped = T )
save(totalRNAcountmatrix2019,file = "UniqueMaptotalRNAcountmatrix2019.RData")
load(file="UniqueMaptotalRNAcountmatrix2019.RData")
Gene2019totalRNAstats=totalRNAcountmatrix2019$stat
#Inspection
View(totalRNAcountmatrix2019$counts)
head(totalRNAcountmatrix2019$annotation)
head(totalRNAcountmatrix2019$stat)
#Counting reads:
#Parameters:
#Not overlapping
#unique reads-not multimappers
#both reads mapped
#Allow metafeature if exon=genes
#MQS=255-mapq score-default used by star
#countReadPairs = F this will now count reads rather than fragments. Results will now be comparable to stats of post alignment quality tools
totalRNAcountmatrix2019reads=featureCounts(totalRNAref2019BAMfiles,annot.ext = gtf2019,isGTFAnnotationFile = T,GTF.featureType = "exon",isPairedEnd = T,countReadPairs = F,countMultiMappingReads = F,minMQS = 255,nthreads =4,useMetaFeatures = T,primaryOnly = T,requireBothEndsMapped = T )
save(totalRNAcountmatrix2019reads,file = "UniqueMaptotalRNAcountmatrix2019reads.RData")
load(file="UniqueMaptotalRNAcountmatrix2019reads.RData")
Gene2019totalRNAreadstats=totalRNAcountmatrix2019reads$stat
write.table(Gene2019totalRNAreadstats,file="UniqueMapGeneTotalRNAcountmatrix2019readstat.txt",row.names = F,col.names = T,quote = F)
#Inspection
View(Gene2019totalRNAreadstats)
sum(Gene2019totalRNAreadstats[,2])
head(totalRNAcountmatrix2019reads$counts)
head(totalRNAcountmatrix2019reads$annotation)
#ribosomal rna
#Parameters:
#Not overlapping
#unique reads-not multimappers
#both reads mapped
#MQS=255-mapq score-default used by star
#countReadPairs = F this will now count reads rather than fragments. Results will now be comparable to stats of post alignment quality tools
rRNAcountmatrix2019= featureCounts(totalRNAref2019BAMfiles,annot.ext = gtf2019,isGTFAnnotationFile = T,GTF.attrType = "ID",GTF.featureType = "rRNA",isPairedEnd = T,countReadPairs = F,countMultiMappingReads = F,minMQS = 255,nthreads =4,primaryOnly = T,useMetaFeatures = F,requireBothEndsMapped = T )
save(rRNAcountmatrix2019,file = "UniqueMaprRNAcountmatrix2019.RData")
load(file="UniqueMaprRNAcountmatrix2019.RData")
rRNA2019totalRNAstats=rRNAcountmatrix2019$stat
write.table(rRNA2019totalRNAstats,file="UniqueMaprRNAcountmatrix2019stat.txt",quote = F,row.names = F,col.names = T)
View(rRNAcountmatrix2019$counts)
rRNAcountmatrix2019$counts
sum(rRNAcountmatrix2019$counts[,2])
View(rRNA2019totalRNAstats)
length(rRNAcountmatrix2019$counts[,1])
#Inspection
View(rRNA2019totalRNAstats)
head(totalRNAcountmatrix2019reads$counts)
#poliRNA
find.file("poliRNA-sensitive1-ref2019Aligned.sortedByCoord.out.bam", dir = "C:/Users/toros/OneDrive - KU Leuven/My thesis/Research/thesis/thesis-results/readmapping/ref2019/poliRNA", dirs = NULL)
find.file("poliRNA-sensitive2-ref2019Aligned.sortedByCoord.out.bam", dir = "C:/Users/toros/OneDrive - KU Leuven/My thesis/Research/thesis/thesis-results/readmapping/ref2019/poliRNA", dirs = NULL)
find.file("poliRNA-sensitive3-ref2019Aligned.sortedByCoord.out.bam", dir = "C:/Users/toros/OneDrive - KU Leuven/My thesis/Research/thesis/thesis-results/readmapping/ref2019/poliRNA", dirs = NULL)
find.file("poliRNA-resistant4-ref2019Aligned.sortedByCoord.out.bam", dir = "C:/Users/toros/OneDrive - KU Leuven/My thesis/Research/thesis/thesis-results/readmapping/ref2019/poliRNA", dirs = NULL)
find.file("poliRNA-resistant5-ref2019Aligned.sortedByCoord.out.bam", dir = "C:/Users/toros/OneDrive - KU Leuven/My thesis/Research/thesis/thesis-results/readmapping/ref2019/poliRNA", dirs = NULL)
find.file("poliRNA-resistant6-ref2019Aligned.sortedByCoord.out.bam", dir = "C:/Users/toros/OneDrive - KU Leuven/My thesis/Research/thesis/thesis-results/readmapping/ref2019/poliRNA", dirs = NULL)
find.file("genome_2019.gff",dir="C:/Users/toros/OneDrive - KU Leuven/My thesis/Research/thesis/thesis-results/readmapping/gtf-files-used",dirs=NULL)

#files to variables:
poliRNA2019sensitive1="C:/Users/toros/OneDrive - KU Leuven/My thesis/Research/thesis/thesis-results/readmapping/ref2019/poliRNA/poliRNA-sensitive1-ref2019Aligned.sortedByCoord.out.bam"
poliRNA2019sensitive2="C:/Users/toros/OneDrive - KU Leuven/My thesis/Research/thesis/thesis-results/readmapping/ref2019/poliRNA/poliRNA-sensitive2-ref2019Aligned.sortedByCoord.out.bam"
poliRNA2019sensitive3="C:/Users/toros/OneDrive - KU Leuven/My thesis/Research/thesis/thesis-results/readmapping/ref2019/poliRNA/poliRNA-sensitive3-ref2019Aligned.sortedByCoord.out.bam"
poliRNA2019resistant4="C:/Users/toros/OneDrive - KU Leuven/My thesis/Research/thesis/thesis-results/readmapping/ref2019/poliRNA/poliRNA-resistant4-ref2019Aligned.sortedByCoord.out.bam"
poliRNA2019resistant5="C:/Users/toros/OneDrive - KU Leuven/My thesis/Research/thesis/thesis-results/readmapping/ref2019/poliRNA/poliRNA-resistant5-ref2019Aligned.sortedByCoord.out.bam"
poliRNA2019resistant6="C:/Users/toros/OneDrive - KU Leuven/My thesis/Research/thesis/thesis-results/readmapping/ref2019/poliRNA/poliRNA-resistant6-ref2019Aligned.sortedByCoord.out.bam"
#Vector of files: (Bam order by coordinates:)
poliRNAref2019BAMfiles<-c(poliRNA2019sensitive1,poliRNA2019sensitive2,poliRNA2019sensitive3,poliRNA2019resistant4,poliRNA2019resistant5,poliRNA2019resistant6)
#annotation file:
gtf2019<-"C:/Users/toros/OneDrive - KU Leuven/My thesis/Research/thesis/thesis-results/readmapping/gtf-files-used/genome_2019.gff"
#counting fragments
#Parameters:
#Not overlapping
#unique reads-not multimappers
#both reads mapped
#MQS=255-mapq score-default used by star
poliRNAcountmatrix2019=featureCounts(poliRNAref2019BAMfiles,annot.ext = gtf2019,isGTFAnnotationFile = T,GTF.featureType = "exon",isPairedEnd = T,countMultiMappingReads = F,minMQS = 255,nthreads =4,primaryOnly = T,useMetaFeatures = F,requireBothEndsMapped = T  )
sum(poliRNAcountmatrix2019$stat[,2])
#Inspection
load(file ="poliRNAcountmatrix2019.RData")
head(poliRNAcountmatrix2019$counts)
head(poliRNAcountmatrix2019$annotation)
head(poliRNAcountmatrix2019$stat)
save(poliRNAcountmatrix2019,file ="poliRNAcountmatrix2019.RData" )
#Counting reads:
#Parameters:
#Not overlapping
#unique reads-not multimappers
#both reads mapped
#MQS=255-mapq score-default used by star
#countReadPairs = F this will now count reads rather than fragments. Results will now be comparable to stats of post alignment quality tools
poliRNAcountmatrix2019reads=featureCounts(poliRNAref2019BAMfiles,annot.ext = gtf2019,isGTFAnnotationFile = T,GTF.featureType = "exon",isPairedEnd = T,countReadPairs = F,countMultiMappingReads = F,minMQS = 255,nthreads =4,primaryOnly = T,useMetaFeatures = F,requireBothEndsMapped = T  )
sum(poliRNAcountmatrix2019reads$stat[,2])
poliRNA2019readstats=poliRNAcountmatrix2019reads$stat
#Inspection
load(file="poliRNAcountmatrix2019reads.RData")
head(poliRNAcountmatrix2019reads$counts)
head(poliRNAcountmatrix2019reads$annotation)
View(poliRNAcountmatrix2019reads$stat)
save(poliRNAcountmatrix2019reads,file ="poliRNAcountmatrix2019reads.RData" )
write.table(poliRNA2019readstats,file="UniqueMapGenePoliRNAcountmatrix2019stat.txt",row.names = F,col.names = T,quote = F)

#ribosomal rna

rRNAcountmatrix2019poliRNA= featureCounts(poliRNAref2019BAMfiles,annot.ext = gtf2019,isGTFAnnotationFile = T,GTF.attrType = "ID",GTF.featureType = "rRNA",isPairedEnd = T,countReadPairs = F,countMultiMappingReads = F,minMQS = 255,nthreads =4,primaryOnly = T,useMetaFeatures = F,requireBothEndsMapped = T )
save(rRNAcountmatrix2019poliRNA,file ="rRNApoliRNAcountmatrix2019.RData" )
load(file ="rRNApoliRNAcountmatrix2019.RData")
rRNA2019PoliRNAstats=rRNAcountmatrix2019poliRNA$stat
#Inspection
load(file="poliRNAcountmatrix2019reads.RData")
head(rRNAcountmatrix2019poliRNA$counts)
head(rRNAcountmatrix2019poliRNA$annotation)
View(rRNAcountmatrix2019poliRNA$stat)
write.table(rRNA2019PoliRNAstats,file="UniqueMaprRNAPoliRNAcountmatrix2019stat.txt",row.names = F,col.names = T,quote = F)

#mapq evaluation
find.file("totalRNA-sensitive2-ref2019Aligned.sortedByCoord.out.bam", dir = "C:/Users/toros/OneDrive - KU Leuven/My thesis/Research/thesis/thesis-results/readmapping/ref2019/totalRNA", dirs = NULL)
find.file("totalRNA-sensitive3-ref2019Aligned.sortedByCoord.out.bam", dir = "C:/Users/toros/OneDrive - KU Leuven/My thesis/Research/thesis/thesis-results/readmapping/ref2019/totalRNA", dirs = NULL)
find.file("totalRNA-sensitive2-ref2019-mapq0Aligned.sortedByCoord.out.bam", dir = "C:/Users/toros/OneDrive - KU Leuven/My thesis/Research/thesis/thesis-results/readmapping/ref2019/totalRNA/mapq0-50", dirs = NULL)
find.file("totalRNA-sensitive3-ref2019-mapq0Aligned.sortedByCoord.out.bam", dir = "C:/Users/toros/OneDrive - KU Leuven/My thesis/Research/thesis/thesis-results/readmapping/ref2019/totalRNA/mapq0-50", dirs = NULL)
find.file("totalRNA-sensitive2-ref2019-mapq50Aligned.sortedByCoord.out.bam", dir = "C:/Users/toros/OneDrive - KU Leuven/My thesis/Research/thesis/thesis-results/readmapping/ref2019/totalRNA/mapq0-50", dirs = NULL)
find.file("totalRNA-sensitive3-ref2019-mapq50Aligned.sortedByCoord.out.bam", dir = "C:/Users/toros/OneDrive - KU Leuven/My thesis/Research/thesis/thesis-results/readmapping/ref2019/totalRNA/mapq0-50", dirs = NULL)
find.file("genome_2019.gff",dir="C:/Users/toros/OneDrive - KU Leuven/My thesis/Research/thesis/thesis-results/readmapping/gtf-files-used",dirs=NULL)

totalsensitive2="C:/Users/toros/OneDrive - KU Leuven/My thesis/Research/thesis/thesis-results/readmapping/ref2019/totalRNA/totalRNA-sensitive2-ref2019Aligned.sortedByCoord.out.bam"
totalsensitive3="C:/Users/toros/OneDrive - KU Leuven/My thesis/Research/thesis/thesis-results/readmapping/ref2019/totalRNA/totalRNA-sensitive3-ref2019Aligned.sortedByCoord.out.bam"
totalsensitive2mapq0="C:/Users/toros/OneDrive - KU Leuven/My thesis/Research/thesis/thesis-results/readmapping/ref2019/totalRNA/mapq0-50/totalRNA-sensitive2-ref2019-mapq0Aligned.sortedByCoord.out.bam"
totalsensitive3mapq0="C:/Users/toros/OneDrive - KU Leuven/My thesis/Research/thesis/thesis-results/readmapping/ref2019/totalRNA/mapq0-50/totalRNA-sensitive3-ref2019-mapq0Aligned.sortedByCoord.out.bam"
totalsensitive2mapq50="C:/Users/toros/OneDrive - KU Leuven/My thesis/Research/thesis/thesis-results/readmapping/ref2019/totalRNA/mapq0-50/totalRNA-sensitive2-ref2019-mapq50Aligned.sortedByCoord.out.bam"
totalsensitive3mapq50="C:/Users/toros/OneDrive - KU Leuven/My thesis/Research/thesis/thesis-results/readmapping/ref2019/totalRNA/mapq0-50/totalRNA-sensitive3-ref2019-mapq50Aligned.sortedByCoord.out.bam"
#mapq255
totalRNAref2019mapq255BAMfiles<-c(totalsensitive2,totalsensitive3)
#annotation file:
gtf2019<-"C:/Users/toros/OneDrive - KU Leuven/My thesis/Research/thesis/thesis-results/readmapping/gtf-files-used/genome_2019.gff"

#Counting reads:
#Parameters:
#Not overlapping
#unique reads-not multimappers
#both reads mapped
#Allow metafeature if exon=genes
#MQS=255-mapq score-default used by star
#countReadPairs = F this will now count reads rather than fragments. Results will now be comparable to stats of post alignment quality tools
totalRNAmapq255countmatrix2019reads=featureCounts(totalRNAref2019mapq255BAMfiles,annot.ext = gtf2019,isGTFAnnotationFile = T,GTF.featureType = "exon",isPairedEnd = T,countReadPairs = F,countMultiMappingReads = F,minMQS = 255,nthreads =4,useMetaFeatures = T,primaryOnly = T,requireBothEndsMapped = T )
save(totalRNAmapq255countmatrix2019reads,file = "UniqueMaptotalRNAcountmatrix2019readsmapq255.RData")
load(file="UniqueMaptotalRNAcountmatrix2019readsmapq255.RData")
Gene2019totalRNAreadmapq255stats=totalRNAmapq255countmatrix2019reads$stat
write.table(Gene2019totalRNAreadmapq255stats,file="UniqueMapGeneTotalRNAcountmatrix2019readmapq255stat.txt",row.names = F,col.names = T,quote = F)
#Inspection
View(Gene2019totalRNAreadmapq255stats)
sum(Gene2019totalRNAreadmapq255stats[,2])
head(totalRNAmapq255countmatrix2019reads$counts)
head(totalRNAmapq255countmatrix2019reads$annotation)
#mapq 0
totalRNAref2019mapq0BAMfiles<-c(totalsensitive2mapq0,totalsensitive3mapq0)
#annotation file:
gtf2019<-"C:/Users/toros/OneDrive - KU Leuven/My thesis/Research/thesis/thesis-results/readmapping/gtf-files-used/genome_2019.gff"

#Counting reads:
#Parameters:
#Not overlapping
#unique reads-not multimappers
#both reads mapped
#Allow metafeature if exon=genes
#MQS=0-mapq score-default used by star
#countReadPairs = F this will now count reads rather than fragments. Results will now be comparable to stats of post alignment quality tools
totalRNAmapq0countmatrix2019reads=featureCounts(totalRNAref2019mapq0BAMfiles,annot.ext = gtf2019,isGTFAnnotationFile = T,GTF.featureType = "exon",isPairedEnd = T,countReadPairs = F,countMultiMappingReads = F,minMQS = 0,nthreads =4,useMetaFeatures = T,primaryOnly = T,requireBothEndsMapped = T )
save(totalRNAmapq0countmatrix2019reads,file = "UniqueMaptotalRNAcountmatrix2019readsmapq0.RData")
load(file="UniqueMaptotalRNAcountmatrix2019readsmapq0.RData")
Gene2019totalRNAreadmapq0stats=totalRNAmapq0countmatrix2019reads$stat
write.table(Gene2019totalRNAreadmapq0stats,file="UniqueMapGeneTotalRNAcountmatrix2019readmapq0stat.txt",row.names = F,col.names = T,quote = F)
#Inspection
View(Gene2019totalRNAreadmapq0stats)
sum(Gene2019totalRNAreadmapq0stats[,2])
head(totalRNAmapq0countmatrix2019reads$counts)
head(totalRNAmapq0countmatrix2019reads$annotation)

#mapq50
totalRNAref2019mapq50BAMfiles<-c(totalsensitive2mapq50,totalsensitive3mapq50)
#annotation file:
gtf2019<-"C:/Users/toros/OneDrive - KU Leuven/My thesis/Research/thesis/thesis-results/readmapping/gtf-files-used/genome_2019.gff"

#Counting reads:
#Parameters:
#Not overlapping
#unique reads-not multimappers
#both reads mapped
#Allow metafeature if exon=genes
#MQS=50-mapq score-default used by star
#countReadPairs = F this will now count reads rather than fragments. Results will now be comparable to stats of post alignment quality tools
totalRNAmapq50countmatrix2019reads=featureCounts(totalRNAref2019mapq50BAMfiles,annot.ext = gtf2019,isGTFAnnotationFile = T,GTF.featureType = "exon",isPairedEnd = T,countReadPairs = F,countMultiMappingReads = F,minMQS = 50,nthreads =4,useMetaFeatures = T,primaryOnly = T,requireBothEndsMapped = T )
save(totalRNAmapq50countmatrix2019reads,file = "UniqueMaptotalRNAcountmatrix2019readsmapq50.RData")
load(file="UniqueMaptotalRNAcountmatrix2019readsmapq50.RData")
Gene2019totalRNAreadmapq50stats=totalRNAmapq50countmatrix2019reads$stat
write.table(Gene2019totalRNAreadmapq50stats,file="UniqueMapGeneTotalRNAcountmatrix2019readmapq50stat.txt",row.names = F,col.names = T,quote = F)
#Inspection
View(Gene2019totalRNAreadmapq50stats)
sum(Gene2019totalRNAreadmapq50stats[,2])
head(totalRNAmapq50countmatrix2019reads$counts)
head(totalRNAmapq50countmatrix2019reads$annotation)