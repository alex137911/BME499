#!/bin/bash
#SBATCH --job-name=download_sra
#SBATCH --output=download_sra.out
#SBATCH --error=download_sra.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G

module load sra-toolkit

# Path to the .sra file
SRA_PATH="/home/a252chan/projects/def-acdoxey/a252chan/SRR15197167/SRR15197167.sra"

fasterq-dump $SRA_PATH -e 32 -O /home/a252chan/projects/def-acdoxey/a252chan/SRR15197167/fastq