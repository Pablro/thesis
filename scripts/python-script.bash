#! /bin/bash
#
#SBATCH --job-name=transcriptional_landscape
#SBATCH --ntasks=1 --cpus-per-task=1 --mem=2g
#SBATCH --time=24:00:00

# Build the environment
module --force purge
module load calcua/all
module load Python/3.10.4-GCCcore-11.3.0
module load Biopython/1.79-foss-2022a
module load SAMtools/1.18-GCC-12.3.0
#export library
export PYTHONPATH="./scripts/python_lib/lib/python3.10/site-packages/:${PYTHONPATH}"

#script
python ./scripts/transcriptome_exploration.py
