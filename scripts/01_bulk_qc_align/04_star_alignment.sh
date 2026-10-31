#!/usr/bin/env bash
# Reproduce bulk STAR alignment settings documented for E-MTAB-5060.
# Run in Bash on a suitably provisioned Linux workstation or cluster.
# Usage: bash 04_star_alignment.sh MODE
# Local Ubuntu paths use the islet-t2d-work bulk directory and ~/reference.
# Modes: index, align, check, summary. No analysis runs without an explicit mode.
set -euo pipefail
usage() {
    printf '%s\n' 'Usage: bash 04_star_alignment.sh {index|align|check|summary}'
    printf '%s\n' 'Optional overrides: PROJECT_ROOT, TRIM_DIR, REF_DIR, INDEX_DIR, BAM_DIR, LOG_DIR, ALIGN_THREADS, INDEX_THREADS'
}
if [[ $# -eq 1 && "$1" == '--help' ]]; then usage; exit 0; fi
if [[ $# -ne 1 ]]; then usage >&2; exit 2; fi
PROJECT_ROOT=${PROJECT_ROOT:-"$HOME/islet-t2d-work"}
MODE=$1
case "$MODE" in index|align|check|summary) ;; *) usage >&2; exit 2 ;; esac
TRIM_DIR=${TRIM_DIR:-"$PROJECT_ROOT/data/processed/bulk_fastq_trimmed"}
REF_DIR=${REF_DIR:-"$HOME/reference"}
INDEX_DIR=${INDEX_DIR:-"$REF_DIR/STAR_index_GRCh38"}
BAM_DIR=${BAM_DIR:-"$PROJECT_ROOT/bulk_bam"}
LOG_DIR=${LOG_DIR:-"$PROJECT_ROOT/results/alignment_logs"}
ALIGN_THREADS=${ALIGN_THREADS:-8}
INDEX_THREADS=${INDEX_THREADS:-16}
FASTA="$REF_DIR/GRCh38.primary_assembly.genome.fa"
GTF="$REF_DIR/gencode.v44.primary_assembly.annotation.gtf"
SAMPLES=(ERR1629982 ERR1629983 ERR1629984 ERR1629985 ERR1629986 ERR1629987 ERR1629988)
fail() { printf 'ERROR: %s\n' "$*" >&2; exit 1; }
require_file() { [[ -s "$1" ]] || fail "Missing or empty file: $1"; }
require_tool() { command -v "$1" >/dev/null 2>&1 || fail "Required command not found: $1"; }
for threads in "$ALIGN_THREADS" "$INDEX_THREADS"; do
    [[ "$threads" =~ ^[1-9][0-9]*$ ]] || fail 'Thread counts must be positive integers.'
done
if [[ "$MODE" == index || "$MODE" == align ]]; then
    require_tool STAR
    [[ $(STAR --version) == '2.7.11b' ]] || fail 'STAR 2.7.11b is required to reproduce this workflow.'
fi
if [[ "$MODE" == align || "$MODE" == check ]]; then require_tool samtools; fi
if [[ "$MODE" == index ]]; then
    # Build a reference index once. Never overwrite an existing nonempty index.
    require_file "$FASTA"
    require_file "$GTF"
    if [[ -d "$INDEX_DIR" && -n $(find "$INDEX_DIR" -mindepth 1 -maxdepth 1 -print -quit) ]]; then
        fail "Index directory is not empty; reuse it with align mode: $INDEX_DIR"
    fi
    mkdir -p "$INDEX_DIR"
    STAR --runMode genomeGenerate \
        --runThreadN "$INDEX_THREADS" \
        --genomeDir "$INDEX_DIR" \
        --genomeFastaFiles "$FASTA" \
        --sjdbGTFfile "$GTF" \
        --sjdbOverhang 100 \
        --genomeSAindexNbases 14 \
        --genomeChrBinNbits 18 \
        --outFileNamePrefix "$INDEX_DIR/index_build_"
    exit 0
fi
if [[ "$MODE" == align ]]; then
    require_tool gzip
    require_tool zcat
    # Validate every pair before starting any sample.
    for file in Genome SA SAindex genomeParameters.txt chrName.txt chrLength.txt; do
        require_file "$INDEX_DIR/$file"
    done
    for sample in "${SAMPLES[@]}"; do
        require_file "$TRIM_DIR/${sample}_1.trimmed.fastq.gz"
        require_file "$TRIM_DIR/${sample}_2.trimmed.fastq.gz"
        # Refuse existing sample outputs, including partial outputs from failed runs.
        if [[ -d "$BAM_DIR" && -n $(find "$BAM_DIR" -maxdepth 1 -name "${sample}_*" -print -quit) ]]; then
            fail "Existing outputs for $sample; choose a new BAM_DIR to rerun."
        fi
    done
    for sample in "${SAMPLES[@]}"; do
        gzip -t "$TRIM_DIR/${sample}_1.trimmed.fastq.gz" "$TRIM_DIR/${sample}_2.trimmed.fastq.gz"
    done
    mkdir -p "$BAM_DIR" "$LOG_DIR"
    session_log=$(mktemp "$LOG_DIR/alignment_run.XXXXXX")
    exec > >(tee "$session_log") 2>&1
    printf 'Project root: %s\nSTAR: %s\n' "$PROJECT_ROOT" "$(STAR --version)"
    samtools --version
    # Process samples sequentially to avoid loading multiple human indexes at once.
    for sample in "${SAMPLES[@]}"; do
        printf '\nAligning %s\n' "$sample"
        STAR --runMode alignReads \
            --runThreadN "$ALIGN_THREADS" \
            --genomeDir "$INDEX_DIR" \
            --readFilesIn "$TRIM_DIR/${sample}_1.trimmed.fastq.gz" "$TRIM_DIR/${sample}_2.trimmed.fastq.gz" \
            --readFilesCommand zcat \
            --outSAMtype BAM SortedByCoordinate \
            --outFilterMultimapNmax 20 \
            --alignSJoverhangMin 8 \
            --alignSJDBoverhangMin 1 \
            --outFilterMismatchNmax 999 \
            --outFilterMismatchNoverReadLmax 0.04 \
            --alignIntronMin 20 \
            --alignIntronMax 1000000 \
            --outFileNamePrefix "$BAM_DIR/${sample}_" \
            --quantMode GeneCounts
        bam="$BAM_DIR/${sample}_Aligned.sortedByCoord.out.bam"
        samtools quickcheck -v "$bam"
        samtools index "$bam"
    done
    printf '\nAlignment and BAM indexing completed. Execution log: %s\n' "$session_log"
    exit 0
fi
if [[ "$MODE" == check ]]; then
    # Basic structural check only; does not read every alignment record.
    bams=()
    for sample in "${SAMPLES[@]}"; do
        bam="$BAM_DIR/${sample}_Aligned.sortedByCoord.out.bam"
        require_file "$bam"
        bams+=("$bam")
    done
    samtools quickcheck -v "${bams[@]}"
    echo 'All seven BAM files passed quickcheck.'
    exit 0
fi
# Summarize existing logs without rerunning alignment or modifying the BAMs.
require_tool python
python - "$BAM_DIR" <<'PYTHON'
import csv
import sys
from pathlib import Path
folder = Path(sys.argv[1])
keys = ['Number of input reads', 'Uniquely mapped reads number',
        'Uniquely mapped reads %', '% of reads mapped to multiple loci',
        '% of reads unmapped: too short', 'Mismatch rate per base, %']
rows = []
for i in range(2, 9):
    sample = f'ERR162998{i}'
    path = folder / f'{sample}_Log.final.out'
    stats = {}
    for line in path.read_text().splitlines():
        if '|' in line:
            key, value = line.split('|', 1)
            stats[key.strip()] = value.strip()
    rows.append([sample, *[stats[key] for key in keys]])
writer = csv.writer(sys.stdout, delimiter='\t', lineterminator='\n')
writer.writerow(['sample_id', *keys])
writer.writerows(rows)
PYTHON
