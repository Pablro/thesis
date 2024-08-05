#! /bin/bash
#
#SBATCH --job-name=postAlignmenQCMapq_2019
#SBATCH --ntasks=24 --cpus-per-task=3 --mem=3g
#SBATCH --time=10:00:00

# Build the environment
module --force purge
module load calcua/all
module load Python/3.11.3-GCCcore-12.3.0
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
bam_stat.py -q 50 -i ./results/totalRNA-sensitive2-ref2019-mapq50Aligned.sortedByCoord.out.bam > ./results/totalRNA-sensitive2-ref2019-mapq50_AligmentStats.txt


bam_stat.py -q 0 -i ./results/totalRNA-sensitive2-ref2019-mapq0Aligned.sortedByCoord.out.bam > ./results/totalRNA-sensitive2-ref2019-mapq0_AligmentStats.txt


bam_stat.py -q 50 -i ./results/totalRNA-sensitive3-ref2019-mapq50Aligned.sortedByCoord.out.bam > ./results/totalRNA-sensitive3-ref2019-mapq50_AligmentStats.txt

bam_stat.py -q 0 -i ./results/totalRNA-sensitive3-ref2019-mapq0Aligned.sortedByCoord.out.bam > ./results/totalRNA-sensitive3-ref2019-mapq0_AligmentStats.txt
