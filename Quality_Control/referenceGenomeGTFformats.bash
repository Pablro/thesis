#! /bin/bash
#The following scripts contain all the commands used to obtain the equivalent files to the gff format  and derivative files, e.g. ribosomal 
#RNA coordinates and the gene prediction file for Picard RNA metrics.
#Note: I run this comands in the terminal. Use this script as reference to adapted it for your needs.
#The gff file was extracted from TriTrypDB
#Requirements: eg. bedops tool installation for gff to bed file conversion
#reference genome 2011
$VSC_DATA/workflows/thesis/tools/bedops/bin/gff2bed $VSC_DATA/workflows/thesis/reference-genome/ref2011/genome_2011.gff $VSC_DATA/workflows/thesis/reference-genome/ref2011/genome_2011.bed
#reference genome 2019
$VSC_DATA/workflows/thesis/tools/bedops/bin/gff2bed $VSC_DATA/workflows/thesis/reference-genome/ref2019/genome_2019.gff $VSC_DATA/workflows/thesis/reference-genome/ref2019/genome_2019.bed
#Summary of features of gtf files
#2011
cut -f8 $VSC_DATA/workflows/thesis/reference-genome/ref2011/genome_2011.bed|sort|uniq -c|awk '{print $2 " " $1}' > $VSC_DATA/workflows/thesis/reference-genome/ref2011/gtf_features_summary.txt
#2019
cut -f8 $VSC_DATA/workflows/thesis/reference-genome/ref2011/genome_2019.bed|sort|uniq -c|awk '{print $2 " " $1}' > $VSC_DATA/workflows/thesis/reference-genome/ref2019/gtf_features_summary.txt
#Building redflat file for Picard
#Requirements: gff3toGenePred tool from UCSC
# raw gene prediction file from the tool 
#ref2011
$VSC_DATA/workflows/thesis/tools/gff3ToGenePred $VSC_DATA/workflows/thesis/reference-genome/ref2011/genome_2011.gff $VSC_DATA/workflows/thesis/reference-genome/ref2011/genome_2011.genePred
#Reflat format according to: 
#https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=1613988045_1uIvC8LkQDyzJIXRJVAXR9o9amv0&hgta_doSchemaDb=hg38&hgta_doSchemaTable=refFlat
awk 'BEGIN { OFS="\t"} {print $12, $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' $VSC_DATA/workflows/thesis/reference-genome/ref2011/genome_2011.genePred > $VSC_DATA/workflows/thesis/reference-genome/ref2011/genome_2011_ref_flat.txt
#ref2019
$VSC_DATA/workflows/thesis/tools/gff3ToGenePred $VSC_DATA/workflows/thesis/reference-genome/ref2019/genome_2019.gff $VSC_DATA/workflows/thesis/reference-genome/ref2019/genome_2019.genePred
awk 'BEGIN { OFS="\t"} {print $12, $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' $VSC_DATA/workflows/thesis/reference-genome/ref2019/genome_2019.genePred > $VSC_DATA/workflows/thesis/reference-genome/ref2019/genome_2019_ref_flat.txt

#Building ribosomal intervals
#Extracting headers-required by ribosomal interval file for Picard
#ref2011
samtools view -H $VSC_DATA/workflows/thesis/results/readmapping/totalRNAsensitive2011/totalRNA-sensitive1-ref2011Aligned.sortedByCoord.out.bam > $VSC_DATA/workflows/thesis/reference-genome/ref2011/headers2011.txt
#ref2019
samtools view -H $VSC_DATA/workflows/thesis/results/readmapping/totalRNAsensitive2019/totalRNA-sensitive1-ref2019Aligned.sortedByCoord.out.bam > $VSC_DATA/workflows/thesis/reference-genome/ref2019/headers2019.txt
#ribosomal coordinates
#ref2011
grep rRNA $VSC_DATA/workflows/thesis/reference-genome/ref2011/genome_2011.gff | awk 'BEGIN{OFS="\t"}{print $9 }'| cut -d';' -f1|cut -d'=' -f2>$VSC_DATA/workflows/thesis/reference-genome/ref2011/intervalname.txt
grep rRNA $VSC_DATA/workflows/thesis/reference-genome/ref2011/genome_2011.gff | awk 'BEGIN{OFS="\t"}{print $1, $4-1, $5}' > $VSC_DATA/workflows/thesis/reference-genome/ref2011/rRNA.bed
paste $VSC_DATA/workflows/thesis/reference-genome/ref2011/rRNA.bed $VSC_DATA/workflows/thesis/reference-genome/ref2011/intervalname.txt > $VSC_DATA/workflows/thesis/reference-genome/ref2011/tmp.txt
mv $VSC_DATA/workflows/thesis/reference-genome/ref2011/tmp.txt $VSC_DATA/workflows/thesis/reference-genome/ref2011/rRNA.bed

#include headers
cat $VSC_DATA/workflows/thesis/reference-genome/ref2011/headers2011.txt $VSC_DATA/workflows/thesis/reference-genome/ref2011/rRNA.bed > $VSC_DATA/workflows/thesis/reference-genome/ref2011/rRNA_interval_list.txt
#ref2019
grep rRNA $VSC_DATA/workflows/thesis/reference-genome/ref2019/genome_2019.gff | awk 'BEGIN{OFS="\t"}{print $9 }'| cut -d';' -f1|cut -d'=' -f2>$VSC_DATA/workflows/thesis/reference-genome/ref2019/intervalname.txt
grep rRNA $VSC_DATA/workflows/thesis/reference-genome/ref2019/genome_2019.gff | awk 'BEGIN{OFS="\t"}{print $1, $4-1, $5}' > $VSC_DATA/workflows/thesis/reference-genome/ref2019/rRNA.bed
paste $VSC_DATA/workflows/thesis/reference-genome/ref2019/rRNA.bed $VSC_DATA/workflows/thesis/reference-genome/ref2019/intervalname.txt > $VSC_DATA/workflows/thesis/reference-genome/ref2019/tmp.txt
mv $VSC_DATA/workflows/thesis/reference-genome/ref2019/tmp.txt $VSC_DATA/workflows/thesis/reference-genome/ref2019/rRNA.bed
#include headers
cat $VSC_DATA/workflows/thesis/reference-genome/ref2019/headers2019.txt $VSC_DATA/workflows/thesis/reference-genome/ref2019/rRNA.bed > $VSC_DATA/workflows/thesis/reference-genome/ref2019/rRNA_interval_list.txt
