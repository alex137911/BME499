#!/bin/bash

#SBATCH --job-name=bowtie2_alignment
#SBATCH --output=bowtie2_alignment.out
#SBATCH --error=bowtie2_alignment.err
#SBATCH --time=04:00:00
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=8
#SBATCH --mem-per-cpu=8G

# Load necessary modules
module load StdEnv/2020
module load bowtie2/2.5.1
module load samtools/1.17

# Define the output directory
OUTPUT_DIR="/home/a252chan/projects/def-acdoxey/a252chan/SRR15197167"

# Align reads using Bowtie2
bowtie2 -x /home/a252chan/projects/def-acdoxey/a252chan/closDifficile/FN545816 \
-1 /home/a252chan/projects/def-acdoxey/a252chan/SRR15197167/fastq/trimmed_SRR15197167_1.fastq \
-2 /home/a252chan/projects/def-acdoxey/a252chan/SRR15197167/fastq/trimmed_SRR15197167_2.fastq \
| samtools view -bS - > /home/a252chan/projects/def-acdoxey/a252chan/SRR15197167/aligned_reads.bam