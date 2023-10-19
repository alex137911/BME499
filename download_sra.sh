#!/bin/bash
#SBATCH --job-name=download_sra
#SBATCH --output=download_sra.out
#SBATCH --error=download_sra.err
#SBATCH --time=12:00:00
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=8
#SBATCH --mem-per-cpu=8G  # Adjust this as necessary

# Load necessary modules
module load StdEnv/2020
module load gcc/9.3.0
module load sra-toolkit/2.10.8

# Path to the .sra file
SRA_PATH="/home/a252chan/projects/def-acdoxey/a252chan/SRR15197167/SRR15197167.sra"

fasterq-dump $SRA_PATH -e 16 -O /home/a252chan/projects/def-acdoxey/a252chan/SRR15197167/fastq