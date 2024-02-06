#example
#Required installation of RSEQC which provide several post alignment measures. 
#Yet I just perform basic stats which are perfectly achieved by Samtools.
#The reason I was hoping to achieve what Picard achieved but with this package. However the  read_distribution.py did not work with the bedfile.
#So I took a longer procedure with Picard 
export PYTHONPATH="${VSC_DATA}/python_lib/lib/python3.11/site-packages/:${PYTHONPATH}"
export PATH="${VSC_DATA}/python_lib/bin/:${PATH}"
#TotalRNAref2011
#Sensitive
#Basic
bam_stat.py -i totalRNA-sensitive1-ref2011Aligned.sortedByCoord.out.bam > totalRNA-sensitive1_AligmentStats.txt
#Picard to evaluated the fraction of reads that map to ribosome.
java -jar $EBROOTPICARD/picard.jar CollectRnaSeqMetrics -I totalRNA-sensitive1-ref2011Aligned.sortedByCoord.out.bam -O totalRNA-sensitive1-ref2011.RNA_Metrics --REF_FLAT $VSC_DATA/workflows/thesis/reference-genome/ref2011/genome_2011_ref_flat.txt --RIBOSOMAL_INTERVALS $VSC_DATA/workflows/thesis/reference-genome/ref2011/rRNA_interval_list.txt -STRAND NONE -CHART totalRNA-sensitive1-ref2011_Picard.pdf
#Resistant
#Basic
bam_stat.py -i totalRNA-sensitive1-ref2011Aligned.sortedByCoord.out.bam > totalRNA-sensitive1_AligmentStats.txt
#Picard to evaluated the fraction of reads that map to ribosome.
java -jar $EBROOTPICARD/picard.jar CollectRnaSeqMetrics -I totalRNA-sensitive1-ref2011Aligned.sortedByCoord.out.bam -O totalRNA-sensitive1-ref2011.RNA_Metrics --REF_FLAT $VSC_DATA/workflows/thesis/reference-genome/ref2011/genome_2011_ref_flat.txt --RIBOSOMAL_INTERVALS $VSC_DATA/workflows/thesis/reference-genome/ref2011/rRNA_interval_list.txt -STRAND NONE -CHART totalRNA-sensitive1-ref2011_Picard.pdf
#PoliRNAref2011
#Sensitive
#Basic
bam_stat.py -i totalRNA-sensitive1-ref2011Aligned.sortedByCoord.out.bam > totalRNA-sensitive1_AligmentStats.txt
#Picard to evaluated the fraction of reads that map to ribosome.
java -jar $EBROOTPICARD/picard.jar CollectRnaSeqMetrics -I totalRNA-sensitive1-ref2011Aligned.sortedByCoord.out.bam -O totalRNA-sensitive1-ref2011.RNA_Metrics --REF_FLAT $VSC_DATA/workflows/thesis/reference-genome/ref2011/genome_2011_ref_flat.txt --RIBOSOMAL_INTERVALS $VSC_DATA/workflows/thesis/reference-genome/ref2011/rRNA_interval_list.txt -STRAND NONE -CHART totalRNA-sensitive1-ref2011_Picard.pdf
#Resistant
#Basic
bam_stat.py -i totalRNA-sensitive1-ref2011Aligned.sortedByCoord.out.bam > totalRNA-sensitive1_AligmentStats.txt
#Picard to evaluated the fraction of reads that map to ribosome.
java -jar $EBROOTPICARD/picard.jar CollectRnaSeqMetrics -I totalRNA-sensitive1-ref2011Aligned.sortedByCoord.out.bam -O totalRNA-sensitive1-ref2011.RNA_Metrics --REF_FLAT $VSC_DATA/workflows/thesis/reference-genome/ref2011/genome_2011_ref_flat.txt --RIBOSOMAL_INTERVALS $VSC_DATA/workflows/thesis/reference-genome/ref2011/rRNA_interval_list.txt -STRAND NONE -CHART totalRNA-sensitive1-ref2011_Picard.pdf
#TotalRNAref2019
#Sensitive
#Basic
bam_stat.py -i totalRNA-sensitive1-ref2011Aligned.sortedByCoord.out.bam > totalRNA-sensitive1_AligmentStats.txt
#Picard to evaluated the fraction of reads that map to ribosome.
java -jar $EBROOTPICARD/picard.jar CollectRnaSeqMetrics -I totalRNA-sensitive1-ref2011Aligned.sortedByCoord.out.bam -O totalRNA-sensitive1-ref2011.RNA_Metrics --REF_FLAT $VSC_DATA/workflows/thesis/reference-genome/ref2011/genome_2011_ref_flat.txt --RIBOSOMAL_INTERVALS $VSC_DATA/workflows/thesis/reference-genome/ref2011/rRNA_interval_list.txt -STRAND NONE -CHART totalRNA-sensitive1-ref2011_Picard.pdf
#Resistant
#Basic
bam_stat.py -i totalRNA-sensitive1-ref2011Aligned.sortedByCoord.out.bam > totalRNA-sensitive1_AligmentStats.txt
#Picard to evaluated the fraction of reads that map to ribosome.
java -jar $EBROOTPICARD/picard.jar CollectRnaSeqMetrics -I totalRNA-sensitive1-ref2011Aligned.sortedByCoord.out.bam -O totalRNA-sensitive1-ref2011.RNA_Metrics --REF_FLAT $VSC_DATA/workflows/thesis/reference-genome/ref2011/genome_2011_ref_flat.txt --RIBOSOMAL_INTERVALS $VSC_DATA/workflows/thesis/reference-genome/ref2011/rRNA_interval_list.txt -STRAND NONE -CHART totalRNA-sensitive1-ref2011_Picard.pdf
#PoliRNAref2019
#Sensitive
#Basic
bam_stat.py -i totalRNA-sensitive1-ref2011Aligned.sortedByCoord.out.bam > totalRNA-sensitive1_AligmentStats.txt
#Picard to evaluated the fraction of reads that map to ribosome.
java -jar $EBROOTPICARD/picard.jar CollectRnaSeqMetrics -I totalRNA-sensitive1-ref2011Aligned.sortedByCoord.out.bam -O totalRNA-sensitive1-ref2011.RNA_Metrics --REF_FLAT $VSC_DATA/workflows/thesis/reference-genome/ref2011/genome_2011_ref_flat.txt --RIBOSOMAL_INTERVALS $VSC_DATA/workflows/thesis/reference-genome/ref2011/rRNA_interval_list.txt -STRAND NONE -CHART totalRNA-sensitive1-ref2011_Picard.pdf
#Resistant
#Basic
bam_stat.py -i totalRNA-sensitive1-ref2011Aligned.sortedByCoord.out.bam > totalRNA-sensitive1_AligmentStats.txt
#Picard to evaluated the fraction of reads that map to ribosome.
java -jar $EBROOTPICARD/picard.jar CollectRnaSeqMetrics -I totalRNA-sensitive1-ref2011Aligned.sortedByCoord.out.bam -O totalRNA-sensitive1-ref2011.RNA_Metrics --REF_FLAT $VSC_DATA/workflows/thesis/reference-genome/ref2011/genome_2011_ref_flat.txt --RIBOSOMAL_INTERVALS $VSC_DATA/workflows/thesis/reference-genome/ref2011/rRNA_interval_list.txt -STRAND NONE -CHART totalRNA-sensitive1-ref2011_Picard.pdf
