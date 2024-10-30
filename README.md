# Investigating the prevalence of antibiotic-resistant Clostridiodes difficile in human gut metagenomes

Metagenome study to assess the global prevalence of antibiotic-resistant Clostridioides difficile – a major contributor to hospital-acquired infections – in human gut metagenomes. Developed a single nucleotide polymorphism-based detection pipeline for antibiotic resistance in *C. difficile*.

### Pipeline Overview
1. **_download_sra.sh_**: Downloads and processes SRA files containing hits for both "_C. difficile_ and "human gut metagenome" (i.e., NCBI Taxonomy ID 408170), using accession IDs from a .CSV file. Converts SRA data to FASTQ format using fasterq-dump for each ID listed.

2. **_fastp_qc.sh_**: Performs quality control and trimming on FASTQ files using fastp for each SRA accession ID. Generates HTML and JSON reports for each ID’s quality control.

3. **_bowtie2_alignment.sh_**: Aligns trimmed FASTQ files to a reference genome (FN545816.1) using Bowtie2 and generates BAM files for each SRA accession ID.

4. **_extracted_alignment.sh_**: Processes BAM files to analyze SNPs in specific genes for each SRA accession ID. Verifies the existence and completeness of BAM files, skipping those with valid output. Sorts and indexes BAM files, then calls variants using bcftools. Checks coverage at key mutation sites (within genes _gyrA_, _rpoB_, _nimB_, _vanR_) and extracts variants if covered. Writes mutation and coverage data to a CSV file, logging any missing or corrupted BAM files.

5. **_variantProcessing.py_**: Processes variant data from _C. difficile_ samples.
