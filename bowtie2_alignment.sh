#!/bin/bash

#SBATCH --job-name=bowtie2_alignment
#SBATCH --output=bowtie2_alignment.out
#SBATCH --error=bowtie2_alignment.err
#SBATCH --time=48:00:00
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=8
#SBATCH --mem-per-cpu=8G

# Load necessary modules
module load StdEnv/2020
module load bowtie2/2.5.1
module load samtools/1.17

# Base directory for SRA data
SRA_DIR="/scratch/a252chan"

# Loop through each SRA ID directory
for SRA_ID_DIR in $SRA_DIR/*/; do
    # Extract the SRA ID from the directory name
    SRA_ID=$(basename $SRA_ID_DIR)

    # Specify the path to the trimmed FASTQ files
    FASTQ_PATH="${SRA_ID_DIR}fastq"

    # Check if trimmed FASTQ files exist
    if [[ -f "${FASTQ_PATH}/trimmed_${SRA_ID}_1.fastq" && -f "${FASTQ_PATH}/trimmed_${SRA_ID}_2.fastq" ]]; then
        # Align reads using Bowtie2
        bowtie2 -x /scratch/a252chan/closDifficile/FN545816 \
        -1 "${FASTQ_PATH}/trimmed_${SRA_ID}_1.fastq" \
        -2 "${FASTQ_PATH}/trimmed_${SRA_ID}_2.fastq" \
        | samtools view -bS - > "${SRA_ID_DIR}aligned_${SRA_ID}_reads.bam"
        
        # Remove untrimmed FASTQ files
        find "${FASTQ_PATH}" -type f -name '*.fastq' ! -name '*trimmed*' -exec rm {} +
    
    else
        echo "Trimmed FASTQ files for $SRA_ID not found."
    fi
done


# # Define the output directory
# OUTPUT_DIR="/scratch/a252chan/SRR13622898"

# # Align reads using Bowtie2
# bowtie2 -x /scratch/a252chan/closDifficile/FN545816 \
# -1 /scratch/a252chan/SRR13622898/fastq/trimmed_SRR13622898_1.fastq \
# -2 /scratch/a252chan/SRR13622898/fastq/trimmed_SRR13622898_2.fastq \
# | samtools view -bS - > $OUTPUT_DIR/aligned_SRR13622898reads.bam