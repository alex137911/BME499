#!/bin/bash

#SBATCH --job-name=fastp_qc
#SBATCH --output=fastp_qc.out
#SBATCH --error=fastp_qc.err
#SBATCH --time=48:00:00
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=8
#SBATCH --mem-per-cpu=8G  # Adjust this as necessary

# Load the fastp module
module load fastp/0.23.4

# Loop through each SRA ID
for SRA_ID in /scratch/a252chan/*/; do
    SRA_ID=$(basename $SRA_ID)
    
    # Delete .SRA files
    rm -r /scratch/a252chan/${SRA_ID}/${SRA_ID}.sra

    # Specify the path to the FASTQ files
    FASTQ_PATH="/scratch/a252chan/${SRA_ID}/fastq"

    # Define the path for the trimmed FASTQ files
    TRIMMED_FASTQ_PATH_1="/scratch/a252chan/${SRA_ID}/fastq/trimmed_${SRA_ID}_1.fastq"
    TRIMMED_FASTQ_PATH_2="/scratch/a252chan/${SRA_ID}/fastq/trimmed_${SRA_ID}_2.fastq"

    # Check if the trimmed FASTQ files exist
    if [[ -f "$TRIMMED_FASTQ_PATH_1" && -f "$TRIMMED_FASTQ_PATH_2" ]]; then
        echo "Trimmed FASTQ files for $SRA_ID already exist, skipping."
        continue
    fi

    # Run fastp for quality control and trimming
    fastp -i "${FASTQ_PATH}/${SRA_ID}_1.fastq" -o "${FASTQ_PATH}/trimmed_${SRA_ID}_1.fastq" \
          -I "${FASTQ_PATH}/${SRA_ID}_2.fastq" -O "${FASTQ_PATH}/trimmed_${SRA_ID}_2.fastq" \
          --html "${FASTQ_PATH}/fastp_report_${SRA_ID}.html" --json "${FASTQ_PATH}/fastp_report_${SRA_ID}.json"
done