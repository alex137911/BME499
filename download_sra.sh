#!/bin/bash

#SBATCH --job-name=download_sra
#SBATCH --output=download_sra.out
#SBATCH --error=download_sra.err
#SBATCH --time=48:00:00
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=8
#SBATCH --mem-per-cpu=8G  # Adjust this as necessary

# Load necessary modules
module load StdEnv/2020
module load gcc/9.3.0
module load sra-toolkit/3.0.0
# module load StdEnv/2023
# module load gcc/12.3
# module load sra-toolkit/3.0.9

CSV_PATH="/home/a252chan/scratch/c_diff-human_gut_metagenomes.csv"
SRA_DIR="/home/a252chan/scratch"

# Create an array of all SRA accession IDs (excluding the header)
mapfile -t SRA_IDS < <(tail -n +2 "$CSV_PATH" | cut -d ',' -f1)

# Loop through each SRA accession ID
for SRA_ID in "${SRA_IDS[@]}"; do
    # Define the path for the trimmed FASTQ files
    TRIMMED_FASTQ_PATH_1="/scratch/a252chan/${SRA_ID}/fastq/trimmed_${SRA_ID}_1.fastq"
    TRIMMED_FASTQ_PATH_2="/scratch/a252chan/${SRA_ID}/fastq/trimmed_${SRA_ID}_2.fastq"

    # Check if the trimmed FASTQ files exist
    if [[ -f "$TRIMMED_FASTQ_PATH_1" && -f "$TRIMMED_FASTQ_PATH_2" ]]; then
        echo "Trimmed FASTQ files for $SRA_ID already exist, skipping."
        continue
    fi

    # Create a directory for each SRA ID
    SRA_ID_DIR="$SRA_DIR/$SRA_ID"
    mkdir -p "$SRA_ID_DIR"

    # Create a directory for the fastq files
    FASTQ_DIR="/scratch/a252chan/$SRA_ID/fastq"
    mkdir -p "$FASTQ_DIR"
    
    # Run fasterq-dump on the SRA ID to convert to FASTQ format
    if ! fasterq-dump "$SRA_ID" -e 16 -O "$FASTQ_DIR"; then
        echo "Failed to run fasterq-dump on $SRA_ID" >&2
    fi
done