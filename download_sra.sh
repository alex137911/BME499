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

CSV_PATH="/home/a252chan/scratch/c_diff-human_gut_metagenomes.csv"
SRA_DIR="/home/a252chan/scratch"

# Create an array of the first 10 SRA accession IDs (excluding the header)
# mapfile -t SRA_IDS < <(tail -n +2 "$CSV_PATH" | cut -d ',' -f1 | head -10)

# Create an array of all SRA accession IDs (excluding the header)
mapfile -t SRA_IDS < <(tail -n +2 "$CSV_PATH" | cut -d ',' -f1)

# Loop through each SRA accession ID
for SRA_ID in "${SRA_IDS[@]}"; do
    # Prefetch .sra files
    # if ! prefetch "$SRA_ID" --output-directory "$SRA_DIR"; then
    #     echo "Failed to prefetch $SRA_ID" >&2
    #     continue  # Skip to the next SRA ID
    # fi
    
    # Define the path to the .sra file
    SRA_PATH="$SRA_DIR/$SRA_ID/$SRA_ID.sra"
    
    # Create a directory for the fastq files
    FASTQ_DIR="/scratch/a252chan/$SRA_ID/fastq"
    mkdir -p "$FASTQ_DIR"
    
    # Run fasterq-dump on the .sra file
    if ! fasterq-dump "$SRA_ID" -e 16 -O "$FASTQ_DIR"; then
        echo "Failed to run fasterq-dump on $SRA_ID" >&2
    fi
done