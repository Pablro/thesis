#! /bin/bash
#SBATCH --job-name=fastqc_analysis
#SBATCH --ntasks=1 --cpus-per-task=6 --mem=2g
#SBATCH --time=40:00
 
#Building environment

#Running fastqc
./tools/fastqc -q ./raw_data/*fastq.gz -o ./results -t 2
