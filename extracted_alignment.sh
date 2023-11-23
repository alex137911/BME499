#!/bin/bash

#SBATCH --job-name=extracted_alignment
#SBATCH --output=extracted_alignment.out
#SBATCH --error=extracted_alignment.err
#SBATCH --time=48:00:00
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=8
#SBATCH --mem-per-cpu=8G

# Load necessary modules
module load StdEnv/2020
module load samtools/1.17
module load bcftools/1.11   # Default Version

# Base directory for SRA data
SRA_DIR="/scratch/a252chan"

# Path to reference genome
REF_GENOME="/home/a252chan/scratch/closDifficile.fasta"

# Output CSV file
CSV_FILE="${SRA_DIR}/mutation_analysis.csv"
echo "SRA_ID,gyrA_mutation,gyrA_coverage,nimB_mutation,nimB_coverage" > "$CSV_FILE"

# Loop through each SRA ID directory
for SRA_ID_DIR in $SRA_DIR/*/; do
    # Extract the SRA ID from the directory name
    SRA_ID=$(basename $SRA_ID_DIR)

    # Define the output directory for the current SRA ID
    OUTPUT_DIR="${SRA_ID_DIR}"

    # Check if the BAM file exists
    if [[ -f "${OUTPUT_DIR}aligned_${SRA_ID}_reads.bam" ]]; then
        # Sort and index the BAM file
        samtools sort "${OUTPUT_DIR}aligned_${SRA_ID}_reads.bam" -o "${OUTPUT_DIR}aligned_${SRA_ID}_reads_sorted.bam"
        samtools index "${OUTPUT_DIR}aligned_${SRA_ID}_reads_sorted.bam"

        # Call variants using bcftools
        bcftools mpileup -f "$REF_GENOME" "${OUTPUT_DIR}aligned_${SRA_ID}_reads_sorted.bam" | bcftools call -mv -Ov -o "${OUTPUT_DIR}${SRA_ID}_variants.vcf"

        # Check coverage and extract variants
        coverage_at_gyrA=$(samtools depth -r FN545816.1:6310 "${OUTPUT_DIR}aligned_${SRA_ID}_reads_sorted.bam" | awk '{print $3}')
        coverage_at_nimB=$(samtools depth -r FN545816.1:1547984 "${OUTPUT_DIR}aligned_${SRA_ID}_reads_sorted.bam" | awk '{print $3}')

        # Checks if there is no coverage at all or if it is 0
        if [[ -z "$coverage_at_gyrA" || "$coverage_at_gyrA" -eq 0 ]]; then
            gyrA_base="Not Covered"
        else
            gyrA_base=$(bcftools query -f '%ALT' -r FN545816.1:6310 "${OUTPUT_DIR}${SRA_ID}_variants.vcf")
        fi

        if [[ -z "$coverage_at_nimB" || "$coverage_at_nimB" -eq 0 ]]; then
            nimB_base="Not Covered"
        else
            nimB_base=$(bcftools query -f '%ALT' -r FN545816.1:1547984 "${OUTPUT_DIR}${SRA_ID}_variants.vcf")
        fi

        # Write to CSV with coverage check
        echo "${SRA_ID},${gyrA_base:-N},${coverage_at_gyrA:-0},${nimB_base:-N},${coverage_at_nimB:-0}" >> "$CSV_FILE"
    else
        echo "BAM file for $SRA_ID not found, skipping."
    fi
done