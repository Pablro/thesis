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
#Counting reads:
#Parameters:
#Defaults in True: Multireads, no overlapping, duplicates not ignore
#Define:
#IsPairedEnd-True

totalRNAcountmatrix2011=featureCounts(totalRNAref2011BAMfiles,annot.ext = gtf2011,isGTFAnnotationFile = T,GTF.featureType = "exon",isPairedEnd = T,nthreads =4 )
#Inspection
head(totalRNAcountmatrix2011$counts)
head(totalRNAcountmatrix2011$annotation)
head(totalRNAcountmatrix2011$stat)
save(totalRNAcountmatrix2011,file ="totalRNAcountmatrix2011.RData" )
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
#Counting reads:
#Parameters:
#Defaults in True: Multireads, no overlapping, duplicates not ignore
#Define:
#IsPairedEnd-True

poliRNAcountmatrix2011=featureCounts(poliRNAref2011BAMfiles,annot.ext = gtf2011,isGTFAnnotationFile = T,GTF.featureType = "exon",isPairedEnd = T,nthreads =4 )
#Inspection
head(poliRNAcountmatrix2011$counts)
View(poliRNAcountmatrix2011$counts)
View(poliRNAcountmatrix2011$annotation)
head(poliRNAcountmatrix2011$stat)
save(poliRNAcountmatrix2011,file ="poliRNAcountmatrix2011.RData" )
