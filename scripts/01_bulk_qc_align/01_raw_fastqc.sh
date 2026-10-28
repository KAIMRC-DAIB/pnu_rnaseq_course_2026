#!/bin/bash

# Stop the script if a command fails, so later steps don't run on incomplete results.
set -e

# Create a folder for the raw-read QC reports.
mkdir -p results/fastqc_raw

# Check each bulk FASTQ file.
fastqc data/raw/bulk_fastq/*.fastq.gz -o results/fastqc_raw -t 8

# Combine the individual FastQC reports.
multiqc results/fastqc_raw -o results/fastqc_raw
