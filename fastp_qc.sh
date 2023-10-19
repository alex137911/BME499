#!/bin/bash

#SBATCH --job-name=fastp_qc
#SBATCH --output=fastp_qc.out
#SBATCH --error=fastp_qc.err
#SBATCH --time=04:00:00
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=4
#SBATCH --mem-per-cpu=4G  # Adjust this as necessary

# Load the fastp module
module load fastp/0.23.4

# Specify the path to the FASTQ files
FASTQ_PATH="/home/a252chan/projects/def-acdoxey/a252chan/SRR15197167/fastq"

# Run fastp for quality control and trimming
fastp -i ${FASTQ_PATH}/SRR15197167.sra_1.fastq -o ${FASTQ_PATH}/trimmed_SRR15197167_1.fastq \
      -I ${FASTQ_PATH}/SRR15197167.sra_2.fastq -O ${FASTQ_PATH}/trimmed_SRR15197167_2.fastq \
      --html ${FASTQ_PATH}/fastp_report.html --json ${FASTQ_PATH}/fastp_report.json