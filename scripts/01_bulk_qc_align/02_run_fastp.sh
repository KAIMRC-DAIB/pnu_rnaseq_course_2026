#!/bin/bash

# Script stops on errors
set -euo pipefail

# Paths are relative to the project folder
cd "$HOME/islet-t2d-project"

# Set short names for variables
raw="data/raw/bulk_fastq"
processed="data/processed/bulk_fastq_trimmed"
reports="results/fastp_reports"

mkdir -p "$processed" "$reports"

# These are the seven bulk runs verified in the SDRF
runs=(ERR1629982 ERR1629983 ERR1629984 ERR1629985 ERR1629986 ERR1629987 ERR1629988)

# Check that both mates exist for every run before processing
for run in "${runs[@]}"; do
  if [[ ! -s "$raw/${run}_1.fastq.gz" || ! -s "$raw/${run}_2.fastq.gz" ]]; then
    echo "Missing or empty FASTQ for $run" >&2
    exit 1
  fi
done

# Process each read-1/red-2 pair together
for run in "${runs[@]}"; do
  echo "Processing $run"
  fastp \
    -i "$raw/${run}_1.fastq.gz" \
    -I "$raw/${run}_2.fastq.gz" \
    -o "$processed/${run}_1.trimmed.fastq.gz" \
    -O "$processed/${run}_2.trimmed.fastq.gz" \
    --detect_adapter_for_pe \
    --qualified_quality_phred 20 \
    --length_required 36 \
    --thread 8 \
    --json "$reports/${run}.json" \
    --html "$reports/${run}.html" \
    2> "$reports/${run}.log"
done

echo "Finished processing all seven paired-ends runs."
