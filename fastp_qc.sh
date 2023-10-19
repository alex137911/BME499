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

# Run fastp for quality control and trimming
fastp -i SRR15197167.sra_1.fastq -o trimmed_SRR15197167_1.fastq \
      -I SRR15197167.sra_2.fastq -O trimmed_SRR15197167_2.fastq \
      --html fastp_report.html --json fastp_report.json