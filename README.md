# Human Pancreatic Islets in Type 2 Diabetes

## Project Overview:

This project investigates gene expression differences between human pancreatic islets from donors with type 2 diabetes (T2D) and non-diabetic donors using publicly available bulk RNA-seq data.

The workflow includes read quality assessment, preprocessing, alignment, gene-level quantification, normalization, exploratory analysis, differential expression analysis, and gene set enrichment analysis (GSEA).

The analysis supports comparison with relevant findings from Segerstolpe et al. (2016). The findings are interpreted in relation to donor coverage and study limitations, taking differences in reference files and analytical methods into account.


## Data Sources:

### RNA-seq Dataset:

The bulk RNA-seq dataset is available under accession **E-MTAB-5060** and is associated with Segerstolpe et al. (2016).

The analysis includes seven whole-islet samples from seven donors: four with T2D and three without diabetes.

| Run accession | Donor | Condition |
|---|---|---|
| ERR1629982 | T2D1 | T2D |
| ERR1629983 | H6 | Non-diabetic |
| ERR1629984 | H3 | Non-diabetic |
| ERR1629985 | H4 | Non-diabetic |
| ERR1629986 | T2D2 | T2D |
| ERR1629987 | T2D3 | T2D |
| ERR1629988 | T2D4 | T2D |

Dataset: https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-5060

FASTQ download links and source metadata are retained in `data/metadata/`.

### Reference Files:

- **Genome:** GRCh38 primary assembly.
- **Annotation:** GENCODE release 44.
- **Alignment index:** STAR index generated from the reference genome and annotation.

The genome, annotation and complete STAR index are stored locally in `reference/` and excluded from GitHub.

The small STAR-generated `geneInfo.tab` file is retained because the GSEA notebook uses it to map gene identifiers to gene symbols.

GENCODE release 44: https://www.gencodegenes.org/human/release_44.html

## Project Organization:

The repository is named `pnu_rnaseq_course_2026`. Its top-level directory is the project root and contains `README.md`, `.gitignore` and the shared configuration.

Paths in `config.yaml` are relative to the project root. The local repository folder can have a different name without changing those relative paths.

| Location | Contents |
|---|---|
| `data/metadata/` | Source metadata, cleaned donor metadata, FASTQ download links and run-to-donor mapping |
| `data/raw/` | Raw sequencing reads stored locally |
| `data/processed/bulk_fastq_trimmed/` | Preprocessed FASTQ files stored locally |
| `data/processed/bulk_bam/` | Alignment files stored locally |
| `data/processed/bulk_counts/Count-matrix/` | Raw gene count matrix and matching sample metadata |
| `data/processed/bulk_counts/featureCounts/` | Gene counts and read-assignment summaries |
| `data/processed/bulk_counts/PyDESeq/Normalization/` | Normalized count matrices and size factors |
| `envs/` | Conda environment specification |
| `reference/` | Local reference resources and the gene mapping file |
| `results/fastqc_raw/` | Raw-read FastQC reports and the raw MultiQC summary |
| `results/fastp_reports/` | Preprocessing reports |
| `results/fastqc_trimmed/` | FastQC reports after preprocessing |
| `results/multiqc_trimmed/` | Combined QC reports after preprocessing |
| `results/figures/` | Figures organized by analysis stage |
| `results/tables/` | Tables organized by analysis stage |
| `results/gsea_native_output/` | Native GSEApy outputs |
| `scripts/01_bulk_qc_align/` | Upstream QC, preprocessing and quantification scripts |
| `scripts/02_bulk_deseq2/` | Normalization, differential expression and GSEA notebooks |
| `config.yaml` | Shared downstream paths and analysis settings |
| `README.md` | Project documentation |
| `.gitignore` | Repository exclusion rules |

Figures and tables use analysis-specific subfolders, including:

- `QC/`
- `normalization/`
- `PCA/`
- `sample_distances/`
- `differential_expression/`
- `gsea/`

Additional exploratory outputs are stored according to the shared configuration.

Folder names are case-sensitive on Linux. Preserve the spelling and capitalization used by the code and configuration.

Excluded local folders may not appear in a fresh GitHub clone. Create or populate them when obtaining data or generating the relevant outputs.

## Environment Setup:

### Requirements:

- Linux, or Ubuntu through WSL2 on Windows.
- Conda.
- A Jupyter-compatible environment for running the downstream notebooks.

### Create and Activate the Environment:

Run these commands from the project root:

```bash
# Create the analysis environment.
conda env create -f envs/bulk-rna-environment.yml

# Activate the environment.
conda activate bulk-rna
```

Create the environment once. In subsequent sessions, activate the existing environment:

```bash
conda activate bulk-rna
```
### Software Versions:

The following versions were recorded in the Linux bulk RNA-seq environment:

| Software | Version | Purpose |
|---|---|---|
| Python | 3.10 | Downstream analysis |
| FastQC | 0.12.1 | Read quality assessment |
| fastp | 0.23.4 | Adapter trimming and read filtering |
| MultiQC | 1.21 | Combined QC reports |
| STAR | 2.7.11b | Read alignment |
| SAMtools | 1.19.2 | Alignment-file processing |
| Subread / featureCounts | 2.0.6 | Gene-level quantification |
| SRA Toolkit | 3.1.1 | Sequencing-data retrieval tools |
| PyDESeq2 | 0.4.10 | Normalization and differential expression |
| GSEApy | 1.1.3 | Gene set enrichment analysis |
| NumPy | 2.2.6 | Numerical operations |
| pandas | 2.3.3 | Data handling |
| SciPy | 1.15.3 | Scientific computing |
| scikit-learn | 1.4.2 | PCA and sample exploration |
| Matplotlib | 3.8.4 | Plotting |
| seaborn | 0.13.2 | Statistical visualization |
| PyYAML | 6.0.3 | Reading the shared configuration |

The environment YAML records the complete dependency list, including
supporting packages. Use the tested specification to recreate the
environment.

### Open the Notebooks:

The notebooks can be opened using Jupyter Notebook or VS Code.

**Jupyter Notebook:**

```bash
# Run from the project root.
jupyter notebook
```

**VS Code — optional:**

```bash
# Open the project from the activated environment.
code .
```

On Windows, open the project through WSL and enable the Python and Jupyter extensions in the WSL environment.

Select the `bulk-rna` environment as the notebook kernel. Activating the environment in the terminal does not automatically select it as the notebook kernel.

To check the interpreter used by a notebook:

```python
import sys

print(sys.executable)
print(sys.version)
```

## Configuration:

The shared `config.yaml` defines downstream input paths, output locations and analysis settings.

Project-relative paths allow team members to keep their repositories in different local locations without inserting personal usernames into the analysis code.

| Section | Purpose |
|---|---|
| `inputs` | Count matrix, sample metadata, gene annotation and differential expression results |
| `normalization` | Normalization output locations |
| `quality_control` | QC output locations |
| `exploration` | Exploratory analysis output locations |
| `pca` | PCA output locations and number of components |
| `sample_distances` | Sample-distance output locations |
| `differential_expression` | Output locations, significance threshold and design |
| `gsea` | Output locations, gene set collection and enrichment settings |

Keep `config.yaml` at the project root.

The upstream Bash scripts define their paths and parameters separately. Review these settings before running the upstream workflow on another computer.

## Donor Metadata:

The files in `data/metadata/` document donor characteristics, disease status and the relationship between sequencing runs and donors.

| File | Purpose |
|---|---|
| `E-MTAB-5060.sdrf.txt` | Original study sample and sequencing information |
| `ERP017123_fastq_links.tsv` | Sequencing-read download links |
| `metadata_bulk.csv.xlsx` | Supplied bulk metadata workbook |
| `donor_metadata_clean.csv` | Cleaned donor characteristics |
| `Run-to-Donor-bulk.csv` | Sequencing run-to-donor mapping |

The downstream sample metadata file is:

```text
data/processed/bulk_counts/Count-matrix/bulk_sample_metadata.tsv
```

Sample identifiers must match the count-matrix columns. Disease labels and sample order must be checked before model fitting.

Each donor represents an independent biological replicate.

## Bulk RNA-seq Workflow:

Run the upstream stages in order. Run the downstream notebooks from top to bottom using the `bulk-rna` kernel.

### 1. Raw Read Quality Assessment:

**Input:** Seven paired-end bulk RNA-seq runs, comprising fourteen FASTQ files, in `data/raw/bulk_fastq/`.

FastQC 0.12.1 assesses each read file. MultiQC 1.21 combines the reports to support comparison across samples.

Run from the project root:

```bash
bash scripts/01_bulk_qc_align/01_raw_fastqc.sh
```

**Outputs:**

- Individual FastQC HTML reports: `results/fastqc_raw/`.
- Combined report: `results/fastqc_raw/multiqc_report.html`.
- Supporting summary data: `results/fastqc_raw/multiqc_data/`.

These reports support assessment of read quality, adapter content and other sequencing characteristics before preprocessing.

### 2. Paired-End Read Preprocessing:

**Input:** Read-1/read-2 FASTQ pairs in `data/raw/bulk_fastq/`.

fastp 0.23.4 processes paired reads together, performs adapter trimming and filters reads using the settings below. The original FASTQ files are preserved.

```bash
bash scripts/01_bulk_qc_align/02_run_fastp.sh
```

| Setting | Value |
|---|---|
| Paired-end adapter detection | Enabled |
| Qualified-base quality threshold | Phred 20 |
| Minimum retained read length | 36 bp |
| Threads | 8 |
| Deduplication | Not enabled |

The Phred threshold defines a qualified base for quality filtering; it does not mean that every base in every retained read has a quality score of at least 20.

**Outputs:**

- Fourteen trimmed paired-end FASTQ files: `data/processed/bulk_fastq_trimmed/`.
- Seven HTML reports, seven JSON reports and execution logs: `results/fastp_reports/`.

The JSON reports record read-retention and filtering statistics. Trimmed FASTQ files remain outside GitHub; the selected reports are included in the repository.

### 3. Quality Assessment After Preprocessing:

**Input:** Trimmed FASTQ files in `data/processed/bulk_fastq_trimmed/`.

FastQC is repeated after preprocessing, followed by MultiQC. The reports support comparison of read quality and adapter content before and after preprocessing.

```bash
bash scripts/01_bulk_qc_align/03_run_fastqc_trimmed.sh
```

**Outputs:**

- Individual FastQC HTML reports: `results/fastqc_trimmed/`.
- Combined report: `results/multiqc_trimmed/multiqc_trimmed.html`.
- Supporting summary data: `results/multiqc_trimmed/multiqc_trimmed_data/`.

FastQC ZIP archives are excluded from the repository. The HTML reports and supporting MultiQC data are retained.

The upstream paths and parameters are defined in the Bash scripts.

### 4. Reference Indexing and Alignment:

**Script:** [04_star_alignment.sh](scripts/01_bulk_qc_align/04_star_alignment.sh)

**Requirements:** the `bulk-rna` environment, STAR v2.7.11b,
samtools, Python, and gzip utilities. Human-genome indexing and
alignment require a workstation or cluster with sufficient RAM and storage.

**Reference files:** place these matching GENCODE v44 / GRCh38
primary-assembly files in `$HOME/reference`:
- `GRCh38.primary_assembly.genome.fa`
- `gencode.v44.primary_assembly.annotation.gtf`

**Default input and output paths:**
- Trimmed FASTQ directory: `$HOME/islet-t2d-work/data/processed/bulk_fastq_trimmed`
- STAR index: `$HOME/reference/STAR_index_GRCh38`
- BAM and STAR outputs: `$HOME/islet-t2d-work/bulk_bam`
- Execution logs: `$HOME/islet-t2d-work/results/alignment_logs`

For each sample from ERR1629982 to ERR1629988, the trimmed
paired-end files must be named:
- `<sample>_1.trimmed.fastq.gz`
- `<sample>_2.trimmed.fastq.gz`

Run the following commands from the repository root:

```bash
# Activate the analysis environment
conda activate bulk-rna

# Build the reference index only if a compatible index is unavailable
bash scripts/01_bulk_qc_align/04_star_alignment.sh index

# Align the seven samples and index their coordinate-sorted BAM files
bash scripts/01_bulk_qc_align/04_star_alignment.sh align

# Check the existing BAM files without rerunning alignment
bash scripts/01_bulk_qc_align/04_star_alignment.sh check

# Display mapping statistics from existing STAR logs
bash scripts/01_bulk_qc_align/04_star_alignment.sh summary
```

Skip `index` when reusing a compatible existing index.
The script refuses to overwrite a nonempty index or existing
sample alignment outputs.

The index uses 16 threads, `sjdbOverhang=100`,
`genomeSAindexNbases=14`, and `genomeChrBinNbits=18`.
Alignment uses eight threads, coordinate-sorted BAM output,
and STAR `GeneCounts`. Other alignment parameters are recorded
in the script.

Paths and thread counts can be overridden through environment
variables. For example:

```bash
# Use trimmed FASTQ files stored outside the default directory
TRIM_DIR="/path/to/trimmed_fastq" \
  bash scripts/01_bulk_qc_align/04_star_alignment.sh align
```

**Outputs:** `<sample>_Aligned.sortedByCoord.out.bam`, its BAM
index, STAR logs, `<sample>_SJ.out.tab`, and
`<sample>_ReadsPerGene.out.tab`. Large reference and alignment
files are not stored in this repository.

### 5. Gene Quantification and Count Matrix Preparation:

**Script:** [05_bulk_featurecounts.sh](scripts/01_bulk_qc_align/05_bulk_featurecounts.sh)

**Requirements:** the `bulk-rna` environment with featureCounts
v2.0.6, samtools, and Python.

**Inputs:**
- Seven coordinate-sorted BAM files from section 4.
- The corresponding STAR `ReadsPerGene.out.tab` files.
- `$HOME/reference/gencode.v44.primary_assembly.annotation.gtf`.
- The cleaned donor table: `$HOME/Downloads/donor_metadata_clean.csv.csv`.
- The original study metadata: `$HOME/Downloads/E-MTAB-5060.sdrf.txt`.

The original donor table and SDRF are required because the script
validates donor metadata as well as preparing the count matrix.
The archived `bulk_sample_metadata.tsv` is an output, not a replacement
for these two inputs.

Run from the repository root:

```bash
# Run counting, matrix validation, and donor-metadata matching
bash scripts/01_bulk_qc_align/05_bulk_featurecounts.sh
```

The script activates `bulk-rna` through
`$HOME/miniforge3/etc/profile.d/conda.sh`. Adjust this activation
path if Conda is installed elsewhere.

The default BAM directory is `$HOME/islet-t2d-work/bulk_bam`.
Input paths can be overridden using `BAM_DIR`, `GTF`,
`DONOR_METADATA`, and `SDRF`.

Counting uses `-T 4 -p --countReadPairs -s 2 -t exon -g gene_id`:
paired-end fragments are counted with reverse-stranded settings,
and exon features are grouped by gene ID. The selected strandedness
reflects the STAR gene-count totals reviewed for this dataset.

The script validates sample order, unique gene identifiers,
nonnegative integer counts, and agreement between matrix column
totals and featureCounts Assigned totals. It also matches donor
metadata to the matrix sample order. Counts remain unnormalized
and unfiltered.

**Generated files** are written to the selected BAM directory.
Rerunning the script replaces files with the same output names.
The repository copies are:
- [Full featureCounts table](data/processed/bulk_counts/featureCounts/bulk_all_gene_counts.txt)
- [Assignment summary](data/processed/bulk_counts/featureCounts/bulk_all_gene_counts.txt.summary)
- [Raw gene-by-donor count matrix](data/processed/bulk_counts/Count-matrix/bulk_count_matrix.tsv)
- [Matched donor metadata](data/processed/bulk_counts/Count-matrix/bulk_sample_metadata.tsv)

For downstream analysis without rerunning alignment or counting,
use the archived count matrix and matched metadata, then continue
with section 6.

### 6. QC, Normalization, and Exploratory Analysis:

> Notebook: [01_bulk_rna_qc_normalization_pca.ipynb](scripts/02_bulk_deseq2/01_bulk_rna_qc_normalization_pca.ipynb)


**Inputs**: Gene quantification and count matrix results (`bulk_count_matrix.tsv`) and (`bulk_sample_metadate.tsv`).

This phase assesses the integrity and uniformity of bulk RNA-seq count data  before downstream analysis

  **Quality control (QC):** Examine raw library sizes, identified genes, zero-count genes, and evaluate the count matrix quality across samples

  **Normalization:** To account for variations in sequencing depth, PyDESeq2 is used to estimate sample-specific size factors and produce normalized gene-expression counts.

  **Log Transformation:** Modify the expression values for exploratory analysis by applying a log2 transformation to normalized counts using **log2(normalized counts + 1)**.

  **Exploratory Filtering:** Retain genes with at least  10 in two samples for exploratory analysis

  **Principal Component Analysis (PCA):** Examines differences in sample associations and gene expression profiles between T2D and non-diabetic groups.

  **Sample Distance Analysis:** Produces a Euclidean distance heatmap to examine the similarities and differences between samples

**Outputs**
- Tables:

QC: `results/tables/QC/qc_metrics.csv`,`results/tables/QC/QC_summary.csv`

Normalization:`results/tables/normalization/library_totals_before_and_after.csv`

PCA: `results/tables/PCA/PCA_coordinates.csv`, `results/tables/PCA/ERR1629984_PC2_gene_contributions.csv`, `results/tables/PCA/PC2_gene_loadings.csv`

Sample Distance: `results/tables/sample_distances/sample_distance_matrix.tsv`

- Figures:

QC: `results/figures/QC/raw_library_size_per_sample.png`

Normalization: `figures/normalization/library_totals.png`

PCA: `results/figures/PCA/PCA_bulk_samples.png`, `results/figures/PCA/PC1 vs hba1c.png`, `results/figures/PCA/PC1 vs bmi.png`, `results/figures/PCA/PC1 vs age.png`, `results/figures/PCA/PC1 by condition.png`, `results/figures/PCA/PC2 by condition.png`, `results/figures/PCA/correlation between Principal Components and Clinical Parameters`

Sample Distance: `results/figures/sample_distances/sample_distance_heatmap.png`

### 7. Differential Expression Analysis:

**Notebook:** [02_bulk_rna_DEseq.ipynb](scripts/02_bulk_deseq2/02_bulk_rna_DEseq.ipynb)

**Requirements:** the `bulk-rna` environment with PyDESeq2 v0.4.10,
Python, pandas, NumPy, Matplotlib, PyYAML, and Jupyter Notebook.

**Inputs:**
- Raw gene count matrix: `data/processed/bulk_counts/Count-matrix/bulk_count_matrix.tsv`
- Sample metadata: `data/processed/bulk_counts/Count-matrix/bulk_sample_metadata.tsv`

This stage performs differential expression analysis to compare
gene expression between T2D and non-diabetic pancreatic islet samples.

Before model fitting, sample identifiers in the count matrix
and metadata are checked and aligned. Count values are checked
to ensure they are numeric, nonnegative, and integer-valued.
Genes with zero counts across all samples are removed before
analysis.

PyDESeq2 is used to fit a single-factor model with `condition`
as the design factor. The analysis tests T2D against
non-diabetic samples using an explicit contrast.

**Analysis Parameters:**

| Parameter | Setting |
|---|---|
| Design factor | `condition` |
| Contrast | T2D vs Non-diabetic |
| Statistical test | Wald test |
| Significance threshold | `0.05` |
| Multiple testing correction | Benjamini-Hochberg (BH) |
| Independent filtering | Enabled (PyDESeq2 default) |
| Cook's distance filtering | Enabled (PyDESeq2 default) |
| Gene filtering | Remove genes with zero counts across all samples |
| Model-fitting CPUs | 2 |

**Outputs:**

- Full differential expression results:
  `results/tables/differential_expression/differential_expression_results.csv`
- Summary table:
  `results/tables/differential_expression/differential_expression_summary.csv`
- Significant genes:
  `results/tables/differential_expression/significant_genes.csv`
- Upregulated genes:
  `results/tables/differential_expression/upregulated_genes.csv`
- Downregulated genes:
  `results/tables/differential_expression/downregulated_genes.csv`

**Figures:**
- Volcano plot: `results/figures/differential_expression/volcano_plot.png`
- MA plot: `results/figures/differential_expression/MA_plot.png`

### 8. Gene Set Enrichment Analysis:


## Reproducing the Analysis:

### From Raw Reads:

1. Obtain the sequencing reads using the dataset accession and recorded download links.
2. Obtain the specified reference genome and annotation.
3. Create and activate the analysis environment.
4. Check the upstream paths and parameters.
5. Run QC, preprocessing, alignment and quantification in order.
6. Generate the combined raw count matrix and matching metadata.
7. Run the downstream notebooks.

Full reproduction from raw reads also requires the alignment/index-generation and count-matrix preparation instructions described above.

### From the Count Matrix:

1. Obtain the raw count matrix and matching sample metadata.
2. Create and activate the environment.
3. Confirm that `geneInfo.tab` and the required gene set file are available.
4. Check `config.yaml` and its integration with the notebooks.
5. Select the `bulk-rna` kernel.
6. Run each notebook from a fresh session, following the workflow order.
7. Confirm that outputs are saved in the intended locations.

Software versions, reference releases, input data, donor selection and analysis parameters must be preserved when reproducing the reported results.

## Files Excluded from GitHub:

The repository excludes large input and intermediate files, including:

- Raw FASTQ files.
- Preprocessed FASTQ files.
- FastQC ZIP archives.
- BAM files and associated alignment indexes.
- The reference genome and annotation.
- STAR index files, except `geneInfo.tab`.
- .DS_Store

Excluded data must be obtained or generated locally. files excluded from the repository are listed in `.gitignore`.

## Contributors:

| Team member | Main analysis contribution |
|---|---|
| Jood Almaymouni | Donor metadata preparation, read quality assessment and preprocessing |
| Aljadel Alazmi | Reference preparation, alignment, gene quantification and count matrix preparation |
| Anar Aldhwaela | QC, normalization, exploratory analysis and PCA |
| Haya Aloshaywan | Differential expression analysis |
| Areej Alownizy | Gene set enrichment analysis |

All team members contribute to workflow review, interpretation, project documentation and report preparation.

## References:

Segerstolpe, Å., Palasantza, A., Eliasson, P., Andersson, E.-M., Andréasson, A.-C., Sun, X., Picelli, S., Sabirsh, A., Clausen, M., Bjursell, M. K., Smith, D. M., Kasper, M., Ämmälä, C., & Sandberg, R. (2016). Single-cell transcriptome profiling of human pancreatic islets in health and type 2 diabetes. *Cell Metabolism, 24*(4), 593–607. https://doi.org/10.1016/j.cmet.2016.08.020

### Software Documentation:

- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [fastp](https://github.com/OpenGene/fastp)
- [MultiQC](https://multiqc.info/)
- [STAR](https://github.com/alexdobin/STAR)
- [Subread / featureCounts](https://subread.sourceforge.net/)
- [SAMtools](https://www.htslib.org/)
- [PyDESeq2](https://pydeseq2.readthedocs.io/)
- [GSEApy](https://gseapy.readthedocs.io/)
