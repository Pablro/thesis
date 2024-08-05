#! /bin/bash
#
#SBATCH --job-name=postAlignmentQC_2019
#SBATCH --ntasks=24 --cpus-per-task=3 --mem=3g
#SBATCH --time=10:00:00

# Build the environment
module --force purge
module load calcua/all
module load picard/2.25.1-Java-11
module load Python/3.11.3-GCCcore-12.3.0
module load SAMtools/1.18-GCC-12.3.0
module load R/4.3.2-gfbf-2023a
#example
#Required installation of RSEQC which provide several post alignment measures. 
#Yet I just perform basic stats which are perfectly achieved by Samtools.
#The reason I was hoping to achieve what Picard achieved but with this package. However the  read_distribution.py did not work with the bedfile.
#So I took a longer procedure with Picard 
export PYTHONPATH="./scripts/python_lib/lib/python3.11/site-packages/:${PYTHONPATH}"
export PATH="./scripts/python_lib/bin/:${PATH}"
#TotalRNAref2011
#Sensitive 1 
#Basic
bam_stat.py -q 255 -i ./results/totalRNA-sensitive1-ref2019Aligned.sortedByCoord.out.bam > ./results/totalRNA-sensitive1_2019_AligmentStats.txt
#Clipping
clipping_profile.py -q 255 -i ./results/totalRNA-sensitive1-ref2019Aligned.sortedByCoord.out.bam -s "PE" -o ./results/totalRNA-sensitive1-ref2019
#Deletion Profile
deletion_profile.py -q 255 -i ./results/totalRNA-sensitive1-ref2019Aligned.sortedByCoord.out.bam -l 151 -o ./results/totalRNA-sensitive1-ref2019
#Insertion profile:
insertion_profile.py -q 255 -s "PE" -i ./results/totalRNA-sensitive1-ref2019Aligned.sortedByCoord.out.bam -o ./results/totalRNA-sensitive1-ref2019
#Picard to evaluated the fraction of reads that map to ribosome.
java -jar $EBROOTPICARD/picard.jar CollectRnaSeqMetrics -I ./results/totalRNA-sensitive1-ref2019Aligned.sortedByCoord.out.bam -O ./results/totalRNA-sensitive1-ref2019.RNA_Metrics --REF_FLAT ./reference-genome/ref2019/genome_2019_ref_flat.txt --RIBOSOMAL_INTERVALS ./reference-genome/ref2019/rRNA_interval_list.txt -STRAND NONE
#Samtools stats:
samtools stats -@3  ./results/totalRNA-sensitive1-ref2019Aligned.sortedByCoord.out.bam > ./results/totalRNA-sensitive1_2019_BamStats.txt
#sensitive 2
#Basic
bam_stat.py -q 255 -i ./results/totalRNA-sensitive2-ref2019Aligned.sortedByCoord.out.bam > ./results/totalRNA-sensitive2_2019_AligmentStats.txt
#Clipping
clipping_profile.py -q 255 -i ./results/totalRNA-sensitive2-ref2019Aligned.sortedByCoord.out.bam -s "PE" -o ./results/totalRNA-sensitive2-ref2019
#Deletion Profile
deletion_profile.py -q 255 -i ./results/totalRNA-sensitive2-ref2019Aligned.sortedByCoord.out.bam -l 151 -o ./results/totalRNA-sensitive2-ref2019
#Insertion profile:
insertion_profile.py -q 255 -s "PE" -i ./results/totalRNA-sensitive2-ref2019Aligned.sortedByCoord.out.bam -o ./results/totalRNA-sensitive2-ref2019
#Picard to evaluated the fraction of reads that map to ribosome.
java -jar $EBROOTPICARD/picard.jar CollectRnaSeqMetrics -I ./results/totalRNA-sensitive2-ref2019Aligned.sortedByCoord.out.bam -O ./results/totalRNA-sensitive2-ref2019.RNA_Metrics --REF_FLAT ./reference-genome/ref2019/genome_2019_ref_flat.txt --RIBOSOMAL_INTERVALS ./reference-genome/ref2019/rRNA_interval_list.txt -STRAND NONE 
#Samtools stats:
samtools stats -@3  ./results/totalRNA-sensitive2-ref2019Aligned.sortedByCoord.out.bam > ./results/totalRNA-sensitive2_2019_BamStats.txt
#sensitive 3
#Basic
bam_stat.py -q 255 -i ./results/totalRNA-sensitive3-ref2019Aligned.sortedByCoord.out.bam > ./results/totalRNA-sensitive3_2019_AligmentStats.txt
#Clipping
clipping_profile.py -q 255 -i ./results/totalRNA-sensitive3-ref2019Aligned.sortedByCoord.out.bam -s "PE" -o ./results/totalRNA-sensitive3-ref2019
#Deletion Profile
deletion_profile.py -q 255 -i ./results/totalRNA-sensitive3-ref2019Aligned.sortedByCoord.out.bam -l 151 -o ./results/totalRNA-sensitive3-ref2019
#Insertion profile:
insertion_profile.py -q 255 -s "PE" -i ./results/totalRNA-sensitive3-ref2019Aligned.sortedByCoord.out.bam -o ./results/totalRNA-sensitive3-ref2019
#Picard to evaluated the fraction of reads that map to ribosome.
java -jar $EBROOTPICARD/picard.jar CollectRnaSeqMetrics -I ./results/totalRNA-sensitive3-ref2019Aligned.sortedByCoord.out.bam -O ./results/totalRNA-sensitive3-ref2019.RNA_Metrics --REF_FLAT ./reference-genome/ref2019/genome_2019_ref_flat.txt --RIBOSOMAL_INTERVALS ./reference-genome/ref2019/rRNA_interval_list.txt -STRAND NONE
#Samtools stats:
samtools stats -@3 ./results/totalRNA-sensitive3-ref2019Aligned.sortedByCoord.out.bam > ./results/totalRNA-sensitive3_2019_BamStats.txt

#Resistant4
#Basic
bam_stat.py -q 255 -i ./results/totalRNA-resistant4-ref2019Aligned.sortedByCoord.out.bam > ./results/totalRNA-resistant4-ref2019_AligmentStats.txt
#Clipping
clipping_profile.py -q 255 -i ./results/totalRNA-resistant4-ref2019Aligned.sortedByCoord.out.bam -s "PE" -o ./results/totalRNA-resistant4-ref2019
#Deletion Profile
deletion_profile.py -q 255 -i ./results/totalRNA-resistant4-ref2019Aligned.sortedByCoord.out.bam -l 151 -o ./results/totalRNA-resistant4-ref2019
#Insertion profile:
insertion_profile.py -q 255 -s "PE" -i ./results/totalRNA-resistant4-ref2019Aligned.sortedByCoord.out.bam -o ./results/totalRNA-resistant4-ref2019
#Picard to evaluated the fraction of reads that map to ribosome.
java -jar $EBROOTPICARD/picard.jar CollectRnaSeqMetrics -I ./results/totalRNA-resistant4-ref2019Aligned.sortedByCoord.out.bam -O ./results/totalRNA-resistant4-ref2019.RNA_Metrics --REF_FLAT ./reference-genome/ref2019/genome_2019_ref_flat.txt --RIBOSOMAL_INTERVALS ./reference-genome/ref2019/rRNA_interval_list.txt -STRAND NONE
#Samtools stats:
samtools stats -@3  ./results/totalRNA-resistant4-ref2019Aligned.sortedByCoord.out.bam > ./results/totalRNA-resistant4-ref2019_BamStats.txt
#Resistant5
#Basic
bam_stat.py -q 255 -i ./results/totalRNA-resistant5-ref2019Aligned.sortedByCoord.out.bam > ./results/totalRNA-resistant5-ref2019_AligmentStats.txt
#Clipping
clipping_profile.py -q 255 -i ./results/totalRNA-resistant5-ref2019Aligned.sortedByCoord.out.bam -s "PE" -o ./results/totalRNA-resistant5-ref2019
#Deletion Profile
deletion_profile.py -q 255 -i ./results/totalRNA-resistant5-ref2019Aligned.sortedByCoord.out.bam -l 151 -o ./results/totalRNA-resistant5-ref2019
#Insertion profile:
insertion_profile.py -q 255 -s "PE" -i ./results/totalRNA-resistant5-ref2019Aligned.sortedByCoord.out.bam -o ./results/totalRNA-resistant5-ref2019
#Picard to evaluated the fraction of reads that map to ribosome.
java -jar $EBROOTPICARD/picard.jar CollectRnaSeqMetrics -I ./results/totalRNA-resistant5-ref2011Aligned.sortedByCoord.out.bam -O ./results/totalRNA-resistant5-ref2019.RNA_Metrics --REF_FLAT ./reference-genome/ref2019/genome_2019_ref_flat.txt --RIBOSOMAL_INTERVALS ./reference-genome/ref2019/rRNA_interval_list.txt -STRAND NONE
#Samtools stats:
samtools stats -@3  ./results/totalRNA-resistant5-ref2019Aligned.sortedByCoord.out.bam > ./results/totalRNA-resistant5-ref2019_BamStats.txt
#Resistant 6
#Basic
bam_stat.py -q 255 -i ./results/totalRNA-resistant6-ref2019Aligned.sortedByCoord.out.bam > ./results/totalRNA-resistant6-ref2019_AligmentStats.txt
#Clipping
clipping_profile.py -q 255 -i ./results/totalRNA-resistant6-ref2019Aligned.sortedByCoord.out.bam -s "PE" -o ./results/totalRNA-resistant6-ref2019
#Deletion Profile
deletion_profile.py -q 255 -i ./results/totalRNA-resistant6-ref2019Aligned.sortedByCoord.out.bam -l 151 -o ./results/totalRNA-resistant6-ref2019
#Insertion profile:
insertion_profile.py -q 255 -s "PE" -i ./results/totalRNA-resistant6-ref2019Aligned.sortedByCoord.out.bam -o ./results/totalRNA-resistant6-ref2019
#Picard to evaluated the fraction of reads that map to ribosome.
java -jar $EBROOTPICARD/picard.jar CollectRnaSeqMetrics -I ./results/totalRNA-resistant6-ref2019Aligned.sortedByCoord.out.bam -O ./results/totalRNA-resistant6-ref2019.RNA_Metrics --REF_FLAT ./reference-genome/ref2019/genome_2019_ref_flat.txt --RIBOSOMAL_INTERVALS ./reference-genome/ref2019/rRNA_interval_list.txt -STRAND NONE
#Samtools stats:
samtools stats -@3  ./results/totalRNA-resistant6-ref2019Aligned.sortedByCoord.out.bam > ./results/totalRNA-resistant6-ref2019_BamStats.txt


