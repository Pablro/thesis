#!/bin/bash
#
#SBATCH --job-name=genome-index
#SBATCH --ntasks=2 --cpus-per-task=5 --mem=4g
#SBATCH --time=24:00:00

# Build the environment
module --force purge
module load calcua/2022a
module load STAR/2.7.10b-GCC-11.3.0



# Run the program
STAR --runThreadN 4 --runMode genomeGenerate --genomeSAindexNbases 11 --genomeDir ./reference-genome/ref2011/STAR --genomeFastaFiles ./reference-genome/ref2011/STAR/2011_Genome.fasta --sjdbGTFfile ./reference-genome/ref2011/STAR/genome_2011.gff
STAR --runThreadN 4 --runMode genomeGenerate --genomeSAindexNbases 11 --genomeDir ./reference-genome/ref2019/STAR --genomeFastaFiles ./reference-genome/ref2019/STAR/2019_Genome.fasta --sjdbGTFfile ./reference-genome/ref2019/STAR/genome_2019.gff
