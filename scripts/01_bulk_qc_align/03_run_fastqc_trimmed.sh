#!/bin/bash

# Stop the script if a command fails, so later steps don't run on incomplete results.
set -e

# Work from the project folder so the paths below are consistent.
cd "$HOME/islet-t2d-project"

trimmed="data/processed/bulk_fastq_trimmed"
fastqc_reports="results/fastqc_trimmed"
multiqc_reports="results/multiqc_trimmed"

# Create folders for the new QC reports.
mkdir -p "$fastqc_reports" "$multiqc_reports"

runs=(ERR1629982 ERR1629983 ERR1629984 ERR1629985 ERR1629986 ERR1629987 ERR1629988)

# Check both trimmed FASTQs for each run and collect their filenames.
files=()
for run in "${runs[@]}"; do
    for mate in 1 2; do
        file="$trimmed/${run}_${mate}.trimmed.fastq.gz"
        if [[ ! -s "$file" ]]; then
            echo "Missing or empty trimmed FASTQ: $file" >&2
            exit 1
        fi
        files+=("$file")
    done
done

echo "Running fastQC on ${#files[@]} trimmed FASTQ files"
fastqc --threads 4 --outdir "$fastqc_reports" "${files[@]}"

echo "Combining the trimmed-read FastQC reports"
multiqc "$fastqc_reports" --outdir "$multiqc_reports" --filename multiqc_trimmed.html

echo "Finished post-trimming FastQC and MultiQC."
