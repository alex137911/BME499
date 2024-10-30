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
module load htslib/1.17     # Default version (for bgzip and tabix)

# Base directory for SRA data
SRA_DIR="/scratch/a252chan"

# Path to reference genome
REF_GENOME="/home/a252chan/scratch/closDifficile.fasta"

# Output CSV file
CSV_FILE="${SRA_DIR}/mutation_analysis.csv"

# Check if CSV file exists and add header if it doesn't
# gyrA mutation/epidemic strain (6310), rpoB mutation/fidaxomicin (97336), 
# nimB mutation/metronidazole (1547984), vanR mutation/vancomycin (1794733)
if [ ! -f "$CSV_FILE" ]; then
    echo "SRA_ID,gyrA_mutation,gyrA_coverage,rpoB_mutation,rpoB_coverage,nimB_mutation,nimB_coverage,vanR_mutation,vanR_coverage" > "$CSV_FILE"
fi

# Loop through each SRA ID directory
for SRA_ID_DIR in $SRA_DIR/*/; do
    # Extract the SRA ID from the directory name
    SRA_ID=$(basename $SRA_ID_DIR)

    # Define the output directory for the current SRA ID
    OUTPUT_DIR="${SRA_ID_DIR}"

    # Check if the BAM file exists
    if [[ -f "${OUTPUT_DIR}aligned_${SRA_ID}_reads.bam" ]]; then

        # Remove any existing temporary files from incomplete runs
        rm -f "${OUTPUT_DIR}aligned_${SRA_ID}_reads_sorted.bam.tmp.*"

        # Check for EOF marker in BAM file (indicative of a complete file)
        if ! samtools quickcheck "${OUTPUT_DIR}aligned_${SRA_ID}_reads.bam"; then
            echo "BAM file for $SRA_ID is corrupted or incomplete, skipping." >> "${SRA_DIR}/error_log.txt"
            continue
        fi

        # Check if the sorted BAM file and compressed VCF file exist (from completed runs, if yes then skip)
        if [[ -f "${OUTPUT_DIR}aligned_${SRA_ID}_reads_sorted.bam" && -f "${OUTPUT_DIR}${SRA_ID}_variants.vcf.gz" ]]; then
            echo "Sample $SRA_ID has already been processed, skipping."
            continue
        fi

        # Sort and index the BAM file
        samtools sort "${OUTPUT_DIR}aligned_${SRA_ID}_reads.bam" -o "${OUTPUT_DIR}aligned_${SRA_ID}_reads_sorted.bam"
        samtools index "${OUTPUT_DIR}aligned_${SRA_ID}_reads_sorted.bam"

        # Call variants using bcftools
        bcftools mpileup -f "$REF_GENOME" "${OUTPUT_DIR}aligned_${SRA_ID}_reads_sorted.bam" | bcftools call -mv -Ov -o "${OUTPUT_DIR}${SRA_ID}_variants.vcf"

        # Compress and index the VCF file
        bgzip -c "${OUTPUT_DIR}${SRA_ID}_variants.vcf" > "${OUTPUT_DIR}${SRA_ID}_variants.vcf.gz"
        tabix -p vcf "${OUTPUT_DIR}${SRA_ID}_variants.vcf.gz"

        # # Check average average depth at relevant positions and extract variants
        # coverage_at_gyrA=$(samtools depth -r FN545816.1:6310 "${OUTPUT_DIR}aligned_${SRA_ID}_reads_sorted.bam" | awk '{sum+=$3; cnt++} END {if (cnt>0) print sum/cnt; else print 0}')
        # coverage_at_rpoB=$(samtools depth -r FN545816.1:97336 "${OUTPUT_DIR}aligned_${SRA_ID}_reads_sorted.bam" | awk '{sum+=$3; cnt++} END {if (cnt>0) print sum/cnt; else print 0}')
        # coverage_at_nimB=$(samtools depth -r FN545816.1:1547984 "${OUTPUT_DIR}aligned_${SRA_ID}_reads_sorted.bam" | awk '{sum+=$3; cnt++} END {if (cnt>0) print sum/cnt; else print 0}')
        # coverage_at_vanR=$(samtools depth -r FN545816.1:1794733 "${OUTPUT_DIR}aligned_${SRA_ID}_reads_sorted.bam" | awk '{sum+=$3; cnt++} END {if (cnt>0) print sum/cnt; else print 0}')

        # Check average average depth 10 bases before and after relevant positions and extract variants
        coverage_at_gyrA=$(samtools depth -r FN545816.1:6300-6320 "${OUTPUT_DIR}aligned_${SRA_ID}_reads_sorted.bam" | awk '{sum+=$3; cnt++} END {if (cnt>0) print sum/cnt; else print 0}')
        coverage_at_rpoB=$(samtools depth -r FN545816.1:97326-97346 "${OUTPUT_DIR}aligned_${SRA_ID}_reads_sorted.bam" | awk '{sum+=$3; cnt++} END {if (cnt>0) print sum/cnt; else print 0}')
        coverage_at_nimB=$(samtools depth -r FN545816.1:1547974-1547994 "${OUTPUT_DIR}aligned_${SRA_ID}_reads_sorted.bam" | awk '{sum+=$3; cnt++} END {if (cnt>0) print sum/cnt; else print 0}')
        coverage_at_vanR=$(samtools depth -r FN545816.1:1794723-1794743 "${OUTPUT_DIR}aligned_${SRA_ID}_reads_sorted.bam" | awk '{sum+=$3; cnt++} END {if (cnt>0) print sum/cnt; else print 0}')

        # Checks if there is no coverage at all or if it is 0
        if [[ -z "$coverage_at_gyrA" || "$coverage_at_gyrA" == "0" ]]; then
            gyrA_base="Not Covered"
        else
            gyrA_base=$(bcftools query -f '%ALT' -r FN545816.1:6310 "${OUTPUT_DIR}${SRA_ID}_variants.vcf.gz")
        fi

        if [[ -z "$coverage_at_rpoB" || "$coverage_at_rpoB" == "0" ]]; then
            rpoB_base="Not Covered"
        else
            rpoB_base=$(bcftools query -f '%ALT' -r FN545816.1:97336 "${OUTPUT_DIR}${SRA_ID}_variants.vcf.gz")
        fi

        if [[ -z "$coverage_at_nimB" || "$coverage_at_nimB" == "0" ]]; then
            nimB_base="Not Covered"
        else
            nimB_base=$(bcftools query -f '%ALT' -r FN545816.1:1547984 "${OUTPUT_DIR}${SRA_ID}_variants.vcf.gz")
        fi

        if [[ -z "$coverage_at_vanR" || "$coverage_at_vanR" == "0" ]]; then
            vanR_base="Not Covered"
        else
            vanR_base=$(bcftools query -f '%ALT' -r FN545816.1:1794733 "${OUTPUT_DIR}${SRA_ID}_variants.vcf.gz")
        fi

        # Write to CSV with coverage check
        echo "${SRA_ID},${gyrA_base:-N},${coverage_at_gyrA:-0},${rpoB_base:-N},${coverage_at_rpoB:-0},${nimB_base:-N},${coverage_at_nimB:-0},${vanR_base:-N},${coverage_at_vanR:-0}" >> "$CSV_FILE"
    else
        echo "BAM file for $SRA_ID not found, skipping."
    fi
done