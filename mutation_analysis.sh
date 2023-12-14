#!/bin/bash

#SBATCH --job-name=mutation_analysis
#SBATCH --output=mutation_analysis.out
#SBATCH --error=mutation_analysis.err
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem-per-cpu=8G

# Load necessary modules
module load StdEnv/2020
module load samtools/1.17
module load bcftools/1.11  # Default version

# Define the reference genome path
REF_GENOME="/home/a252chan/scratch/closDifficile.fasta"

# Define the output directory for the specific SRA ID
OUTPUT_DIR="/home/a252chan/scratch/ERR2835282/"

# SRA ID
SRA_ID="ERR2835282"

echo "Starting variant calling for $SRA_ID"

# Run bcftools commands
bcftools mpileup -f "$REF_GENOME" "${OUTPUT_DIR}aligned_${SRA_ID}_reads_sorted.bam" | \
bcftools call -mv -Ov -o "${OUTPUT_DIR}${SRA_ID}_variants.vcf"

echo "Variant calling completed for $SRA_ID"