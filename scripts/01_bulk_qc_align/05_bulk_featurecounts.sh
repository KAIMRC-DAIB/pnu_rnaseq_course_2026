#!/usr/bin/env bash

# Stop on command failures, undefined variables, and pipeline errors.
set -euo pipefail

# Activate the bulk RNA-seq environment.
source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate bulk-rna

# Define the original local input paths.
# Override these variables when running on another computer.
BAM_DIR="${BAM_DIR:-$HOME/islet-t2d-work/bulk_bam}"
GTF="${GTF:-$HOME/reference/gencode.v44.primary_assembly.annotation.gtf}"
DONOR_METADATA="${DONOR_METADATA:-$HOME/Downloads/donor_metadata_clean.csv.csv}"
SDRF="${SDRF:-$HOME/Downloads/E-MTAB-5060.sdrf.txt}"

export DONOR_METADATA SDRF

# Enter the existing BAM directory and define the sample order.
cd "$BAM_DIR"
BAMS=(ERR162998{2..8}_Aligned.sortedByCoord.out.bam)

# Verify all required inputs before running the analysis.
for file in "$GTF" "$DONOR_METADATA" "$SDRF" "${BAMS[@]}"; do
    if [[ ! -s "$file" ]]; then
        printf 'Missing or empty input: %s\n' "$file" >&2
        exit 1
    fi
done

# Display the installed version and BAM file sizes.
featureCounts -v
ls -lh "${BAMS[@]}"

# Check BAM headers and end-of-file markers.
# This is a basic integrity check, not a full record-by-record check.
samtools quickcheck -v "${BAMS[@]}"
echo "BAM quickcheck passed"

# Display STAR gene-count totals under each strandedness assumption.
# Exclude the first four summary rows.
for sample in ERR162998{2..8}; do
    file="${sample}_ReadsPerGene.out.tab"
    if [[ ! -s "$file" ]]; then
        printf 'Missing STAR gene-count file: %s\n' "$file" >&2
        exit 1
    fi

    awk '
        NR > 4 {
            unstranded += $2
            forward += $3
            reverse += $4
        }
        END {
            print FILENAME
            print "  Unstranded:", unstranded
            print "  Forward:", forward
            print "  Reverse:", reverse
        }
    ' "$file"
done

# Reverse-stranded counting was selected after reviewing these totals.
# Count paired-end fragments across exons grouped by gene_id.
# Rerunning this script replaces the count and derived output files.
featureCounts \
    -T 4 \
    -p \
    --countReadPairs \
    -s 2 \
    -t exon \
    -g gene_id \
    -a "$GTF" \
    -o bulk_all_gene_counts.txt \
    "${BAMS[@]}"

# Summarize assignment, validate the matrix, and match donor metadata.
python - <<'PY'
import csv
import os
from pathlib import Path

expected = [f"ERR162998{i}" for i in range(2, 9)]

# Read assignment statistics.
with Path("bulk_all_gene_counts.txt.summary").open() as f:
    summary = list(csv.reader(f, delimiter="\t"))

summary_samples = [
    Path(name).name.split("_")[0] for name in summary[0][1:]
]
assert summary_samples == expected, "Unexpected summary sample order"
assert all(len(row) == 8 for row in summary), "Invalid summary columns"
status = {row[0]: [int(v) for v in row[1:]] for row in summary[1:]}

for i, sample in enumerate(expected):
    total = sum(values[i] for values in status.values())
    assert total > 0, f"No counting records: {sample}"
    assigned = status["Assigned"][i]
    print(f"\n{sample}")
    print(f"  Assigned: {assigned:,} ({100 * assigned / total:.1f}%)")
    for category in (
        "Unassigned_MultiMapping",
        "Unassigned_NoFeatures",
        "Unassigned_Ambiguity",
    ):
        print(f"  {category}: {100 * status[category][i] / total:.1f}%")

# Retain raw counts and versioned Ensembl gene identifiers.
data = []
seen = set()
totals = [0] * len(expected)

with Path("bulk_all_gene_counts.txt").open() as f:
    reader = csv.reader(
        (line for line in f if not line.startswith("#")),
        delimiter="\t",
    )
    header = next(reader)
    samples = [Path(name).name.split("_")[0] for name in header[6:]]
    assert samples == expected, f"Unexpected samples: {samples}"

    for row in reader:
        assert len(row) == 13, "Unexpected count-table column count"
        gene_id = row[0]
        assert gene_id and gene_id not in seen, f"Invalid gene ID: {gene_id}"
        seen.add(gene_id)
        assert all(v.isdigit() for v in row[6:]), (
            f"Invalid nonnegative integer count: {gene_id}"
        )
        counts = [int(v) for v in row[6:]]
        totals = [a + b for a, b in zip(totals, counts)]
        data.append([gene_id, *counts])

assert len(data) == 62754, f"Unexpected gene count: {len(data)}"
assert totals == status["Assigned"], "Matrix totals do not match Assigned"

# Read cleaned donor metadata.
with Path(os.environ["DONOR_METADATA"]).open(encoding="utf-8-sig") as f:
    donors = [
        {k.strip(): (v or "").strip() for k, v in row.items()}
        for row in csv.DictReader(f)
    ]

by_run = {}
for donor in donors:
    run = donor["Bulk ENA run ID"]
    if run:
        assert run not in by_run, f"Duplicate donor mapping: {run}"
        by_run[run] = donor

# Read study metadata and check repeated paired-read records.
fields = [
    "Characteristics[individual]",
    "Characteristics[disease]",
    "Characteristics[sex]",
    "Characteristics[age]",
    "Characteristics[body mass index]",
]
study = {}
with Path(os.environ["SDRF"]).open(encoding="utf-8-sig") as f:
    for row in csv.DictReader(f, delimiter="\t"):
        run = row["Comment[ENA_RUN]"].strip()
        values = tuple(row[key].strip() for key in fields)
        if run in study:
            assert study[run] == values, f"Conflicting SDRF metadata: {run}"
        study[run] = values

conditions = {
    "normal": "Non-diabetic",
    "type II diabetes mellitus": "T2D",
}

# Match metadata to the matrix sample order.
metadata_rows = []
for run in samples:
    assert run in by_run and run in study, f"Missing metadata: {run}"
    donor = by_run[run]
    source_id, disease, sex, age, bmi = study[run]
    condition = conditions[disease]

    assert donor["Bulk"].lower() == "yes", f"Not marked as bulk: {run}"
    assert donor["Disease"] == condition, f"Disease mismatch: {run}"
    assert donor["Sex"].lower() == sex.lower(), f"Sex mismatch: {run}"
    assert float(donor["Age"]) == float(age), f"Age mismatch: {run}"

    if donor["Source Donor ID"]:
        assert donor["Source Donor ID"] == source_id, (
            f"Source donor ID mismatch: {run}"
        )

    # BMI comes from the SDRF; HbA1c comes from the cleaned donor table.
    metadata_rows.append([
        run, donor["Donor"], source_id, condition,
        sex.lower(), age, bmi, donor["HbA1c"],
    ])

assert len({row[1] for row in metadata_rows}) == 7, "Duplicate donor IDs"
assert len({row[2] for row in metadata_rows}) == 7, "Duplicate source IDs"

# Write derived files only after the validation checks pass.
with Path("bulk_count_matrix.tsv").open("w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerow(["gene_id", *samples])
    writer.writerows(data)

with Path("bulk_sample_metadata.tsv").open("w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerow([
        "sample_id", "donor_id", "source_donor_id", "condition",
        "sex", "age", "bmi", "hba1c",
    ])
    writer.writerows(metadata_rows)

print(f"\nGenes: {len(data):,} | Samples: {len(samples)}")
print("Count-matrix and metadata validation passed")
for row in metadata_rows:
    print(f"{row[0]} -> {row[1]} | {row[3]}")
PY

# Display the final output files.
ls -lh \
    bulk_all_gene_counts.txt \
    bulk_all_gene_counts.txt.summary \
    bulk_count_matrix.tsv \
    bulk_sample_metadata.tsv

echo "Bulk counting and matrix preparation completed"
