#! /bin/bash
#
#SBATCH --job-name=read-mapping-ref2011
#SBATCH --ntasks=1 --cpus-per-task=6 --mem=4g
#SBATCH --time=10:00:00

# Build the environment
module --force purge
module load calcua/2022a
module load STAR/2.7.10b-GCC-11.3.0


# Star alignment
# reference genome 2011: total rna seq
# Bam file per sample to further do count matrix
STAR --runThreadN 5  --genomeDir ./reference-genome/ref2019/STAR --outFileNamePrefix ./results/totalRNA-sensitive1-ref2019 --outSAMtype BAM Unsorted --readFilesIn ./raw_data/totalRNA/sample1R1.fastq ./raw_data/totalRNA/sample1R2.fastq
#STAR --runThreadN 5  --genomeDir ./reference-genome/ref2019/STAR --outFileNamePrefix ./results/totalRNA-sensitive2-ref2019 --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 6000000000 --outBAMsortingBinsN 200 --readFilesCommand zcat --readFilesIn ./raw_data/totalRNA/sample2R1.fastq.gz ./raw_data/totalRNA/sample2R2.fastq.gz
#STAR --runThreadN 5  --genomeDir ./reference-genome/ref2019/STAR --outFileNamePrefix ./results/totalRNA-sensitive3-ref2019 --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 6000000000 --outBAMsortingBinsN 200 --readFilesCommand zcat --readFilesIn ./raw_data/totalRNA/sample3R1.fastq.gz ./raw_data/totalRNA/sample3R2.fastq.gz
#STAR --runThreadN 5  --genomeDir ./reference-genome/ref2019/STAR --outFileNamePrefix ./results/totalRNA-resistant4-ref2019 --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 6000000000 --outBAMsortingBinsN 200 --readFilesCommand zcat --readFilesIn ./raw_data/totalRNA/sample4R1.fastq.gz ./raw_data/totalRNA/sample4R2.fastq.gz
#STAR --runThreadN 5  --genomeDir ./reference-genome/ref2019/STAR --outFileNamePrefix ./results/totalRNA-resistant5-ref2019 --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 6000000000 --outBAMsortingBinsN 200 --readFilesCommand zcat --readFilesIn ./raw_data/totalRNA/sample5R1.fastq.gz ./raw_data/totalRNA/sample5R2.fastq.gz
#STAR --runThreadN 5  --genomeDir ./reference-genome/ref2019/STAR --outFileNamePrefix ./results/totalRNA-resistant6-ref2019 --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 6000000000 --outBAMsortingBinsN 200 --readFilesCommand zcat --readFilesIn ./raw_data/totalRNA/sample6R1.fastq.gz ./raw_data/totalRNA/sample6R2.fastq.gz
#reference genome 2011: polienrinched rna seq
# Bam file per sample to further do count matrix
#STAR --runThreadN 5 --genomeDir ./reference-genome/ref2019/STAR --outFileNamePrefix ./results/poliRNA-sensitive1-ref2019 --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 6000000000 --outBAMsortingBinsN 200 --readFilesCommand zcat --readFilesIn ./raw_data/polienrichedRNA/sample1enrichedR1.fastq.gz ./raw_data/polienrichedRNA/sample1enrichedR2.fastq.gz
#STAR --runThreadN 5 --genomeDir ./reference-genome/ref2019/STAR --outFileNamePrefix ./results/poliRNA-sensitive2-ref2019 --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 6000000000 --outBAMsortingBinsN 200 --readFilesCommand zcat --readFilesIn ./raw_data/polienrichedRNA/sample2enrichedR1.fastq.gz ./raw_data/polienrichedRNA/sample2enrichedR2.fastq.gz
#STAR --runThreadN 5 --genomeDir ./reference-genome/ref2019/STAR --outFileNamePrefix ./results/poliRNA-sensitive3-ref2019 --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 6000000000 --outBAMsortingBinsN 200 --readFilesCommand zcat --readFilesIn ./raw_data/polienrichedRNA/sample3enrichedR1.fastq.gz ./raw_data/polienrichedRNA/sample3enrichedR2.fastq.gz
#STAR --runThreadN 5 --genomeDir ./reference-genome/ref2019/STAR --outFileNamePrefix ./results/poliRNA-resistant4-ref2019 --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 6000000000 --outBAMsortingBinsN 200 --readFilesCommand zcat --readFilesIn ./raw_data/polienrichedRNA/sample4enrichedR1.fastq.gz ./raw_data/polienrichedRNA/sample4enrichedR2.fastq.gz
#STAR --runThreadN 5 --genomeDir ./reference-genome/ref2019/STAR --outFileNamePrefix ./results/poliRNA-resistant5-ref2019 --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 6000000000 --outBAMsortingBinsN 200 --readFilesCommand zcat --readFilesIn ./raw_data/polienrichedRNA/sample5enrichedR1.fastq.gz ./raw_data/polienrichedRNA/sample5enrichedR2.fastq.gz
#STAR --runThreadN 5 --genomeDir ./reference-genome/ref2019/STAR --outFileNamePrefix ./results/poliRNA-resistant6-ref2019 --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 6000000000 --outBAMsortingBinsN 200 --readFilesCommand zcat --readFilesIn ./raw_data/polienrichedRNA/sample6enrichedR1.fastq.gz ./raw_data/polienrichedRNA/sample6enrichedR2.fastq.gz





