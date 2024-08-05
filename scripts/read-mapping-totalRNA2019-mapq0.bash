#! /bin/bash
#
#SBATCH --job-name=read-mapping-ref2019-mapq50
#SBATCH --ntasks=1 --cpus-per-task=6 --mem=4g
#SBATCH --time=24:00:00

# Build the environment
module --force purge
module load calcua/2022a
module load STAR/2.7.10b-GCC-11.3.0


# Star alignment
# reference genome 2019: total rna seq
# Bam file per sample
STAR --runThreadN 5  --genomeDir ./reference-genome/ref2019/STAR --outFileNamePrefix ./results/totalRNA-sensitive1-ref2019-mapq50 --outSAMtype BAM SortedByCoordinate --outSAMmapqUnique 50 --readFilesIn ./raw_data/totalRNA/sample1R1.fastq ./raw_data/totalRNA/sample1R2.fastq
STAR --runThreadN 5  --genomeDir ./reference-genome/ref2019/STAR --outFileNamePrefix ./results/totalRNA-sensitive2-ref2019-mapq50 --outSAMtype BAM SortedByCoordinate --outSAMmapqUnique 50 --limitBAMsortRAM 6000000000 --outBAMsortingBinsN 200 --readFilesCommand zcat --readFilesIn ./raw_data/totalRNA/sample2R1.fastq.gz ./raw_data/totalRNA/sample2R2.fastq.gz
STAR --runThreadN 5  --genomeDir ./reference-genome/ref2019/STAR --outFileNamePrefix ./results/totalRNA-sensitive3-ref2019-mapq50 --outSAMtype BAM SortedByCoordinate --outSAMmapqUnique 50 --limitBAMsortRAM 6000000000 --outBAMsortingBinsN 200 --readFilesCommand zcat --readFilesIn ./raw_data/totalRNA/sample3R1.fastq.gz ./raw_data/totalRNA/sample3R2.fastq.gz
STAR --runThreadN 5  --genomeDir ./reference-genome/ref2019/STAR --outFileNamePrefix ./results/totalRNA-resistant4-ref2019-mapq50 --outSAMtype BAM SortedByCoordinate --outSAMmapqUnique 50 --limitBAMsortRAM 6000000000 --outBAMsortingBinsN 200 --readFilesCommand zcat --readFilesIn ./raw_data/totalRNA/sample4R1.fastq.gz ./raw_data/totalRNA/sample4R2.fastq.gz
STAR --runThreadN 5  --genomeDir ./reference-genome/ref2019/STAR --outFileNamePrefix ./results/totalRNA-resistant5-ref2019-mapq50 --outSAMtype BAM SortedByCoordinate --outSAMmapqUnique 50 --limitBAMsortRAM 6000000000 --outBAMsortingBinsN 200 --readFilesCommand zcat --readFilesIn ./raw_data/totalRNA/sample5R1.fastq.gz ./raw_data/totalRNA/sample5R2.fastq.gz
STAR --runThreadN 5  --genomeDir ./reference-genome/ref2019/STAR --outFileNamePrefix ./results/totalRNA-resistant6-ref2019-mapq50 --outSAMtype BAM SortedByCoordinate --outSAMmapqUnique 50 --limitBAMsortRAM 6000000000 --outBAMsortingBinsN 200 --readFilesCommand zcat --readFilesIn ./raw_data/totalRNA/sample6R1.fastq.gz ./raw_data/totalRNA/sample6R2.fastq.gz





