#!/bin/bash

#SBATCH --job-name=extracted_alignment
#SBATCH --output=extracted_alignment.out
#SBATCH --error=extracted_alignment.err
#SBATCH --time=24:00:00
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=8
#SBATCH --mem-per-cpu=8G

# Load necessary modules
module load StdEnv/2020
module load samtools/1.17

# Define the output directory
OUTPUT_DIR="/home/a252chan/scratch/SRR13622898"

# Check integrity of BAM file
samtools quickcheck $OUTPUT_DIR/aligned_SRR13622898reads.bam

# Sort the BAM file
samtools sort $OUTPUT_DIR/aligned_SRR13622898reads.bam -o $OUTPUT_DIR/aligned_SRR13622898reads_sorted.bam

# Index the sorted BAM file
samtools index $OUTPUT_DIR/aligned_SRR13622898reads_sorted.bam

# Retrieve alignments specific to the reference genome "FN545816.1"
samtools view -b $OUTPUT_DIR/aligned_SRR13622898reads_sorted.bam "FN545816.1" > $OUTPUT_DIR/extracted_SRR13622898reads.bam

# Index the extracted BAM file (to view in IGV)
samtools index $OUTPUT_DIR/extracted_SRR13622898reads.bam

# # Define the output directory
# OUTPUT_DIR="/home/a252chan/projects/def-acdoxey/a252chan/SRR13622898"