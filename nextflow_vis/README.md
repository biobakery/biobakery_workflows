# bioBakery Visualization Workflow — Nextflow

A **Nextflow DSL2** port of the [bioBakery workflows](https://github.com/biobakery/biobakery_workflows) AnADAMA2 visualization pipeline (`vis.py`).

Generates publication-quality **PDF or HTML reports** from whole-metagenome shotgun sequencing (WGX) or 16S rRNA amplicon (16S) data processed by bioBakery workflows. All original visualizations are preserved exactly.

---

## Table of Contents

1. [Features](#features)
2. [Requirements](#requirements)
3. [Installation](#installation)
4. [Quick Start](#quick-start)
5. [Sample Run — HMP2 IBD Cohort](#sample-run--hmp2-ibd-cohort)
6. [Usage](#usage)
7. [Parameters](#parameters)
8. [Execution Environments](#execution-environments)
9. [Output Description](#output-description)
10. [Monitoring (Web Dashboard)](#monitoring-web-dashboard)
11. [Multi-Project Runs](#multi-project-runs)
12. [Troubleshooting](#troubleshooting)

---

## Features

| Feature | AnADAMA2 (`vis.py`) | Nextflow (`main.nf`) |
|---------|---------------------|----------------------|
| WGX (metagenomics) report | ✓ | ✓ |
| 16S amplicon report | ✓ | ✓ |
| Alpha diversity (R/vegan) | ✓ | ✓ |
| PCoA ordination plots | ✓ | ✓ |
| Taxonomy heatmaps (Spearman, Bray-Curtis, Z-score, Filtered) | ✓ | ✓ |
| Stacked / grouped / average barplots | ✓ | ✓ |
| Functional profiling (pathways, ECs) | ✓ | ✓ |
| Quality control read count plots | ✓ | ✓ |
| Metadata stratification | ✓ | ✓ |
| EC name annotation | ✓ | ✓ |
| PDF or HTML output | ✓ | ✓ |
| ZIP archive of all outputs | ✓ | ✓ |
| Docker / Singularity containers | — | ✓ |
| Conda environment | — | ✓ |
| SLURM / PBS / LSF / SGE | — | ✓ |
| AWS Batch | — | ✓ |
| Google Cloud Batch | — | ✓ |
| Multi-project parallelism | — | ✓ |
| Automatic retry on failure | — | ✓ |
| HTML timeline / trace reports | — | ✓ |
| **Web monitoring dashboard** | — | ✓ |
| **Simple primitives (floor/ceiling)** | — | ✓ |
| **Live SLURM job tracking** | — | ✓ |

---

## Requirements

### Software

| Tool | Version | Notes |
|------|---------|-------|
| [Nextflow](https://www.nextflow.io/) | ≥ 23.04 | `curl -s https://get.nextflow.io \| bash` |
| [biobakery_workflows](https://github.com/biobakery/biobakery_workflows) | ≥ 3.0 | Pre-installed in the Docker image |
| [anadama2](https://github.com/biobakery/anadama2) | ≥ 0.10.0 | Dependency of biobakery_workflows |
| Python | ≥ 3.8 | |
| R | ≥ 4.0 | |
| Pandoc | ≥ 2.11 | For PDF/HTML conversion |
| LaTeX (`texlive`) | — | For PDF output only |

**With Docker or Singularity, no manual installation is needed.** The container includes everything.

### Input Data

The `--input` directory must be the output folder of a bioBakery data-processing workflow:

**WGX (Whole Metagenome Shotgun)** — produced by `biobakery_workflows wmgx`:

| File (pattern) | Required | Description |
|----------------|----------|-------------|
| `*_taxonomic_profiles.tsv` | ✓ | MetaPhlAn merged taxonomic profiles |
| `*_pathabundance*.tsv` | Optional | HUMAnN merged pathway abundances |
| `*_ecs*.tsv` | Optional | HUMAnN merged EC abundances |
| `*kneaddata_read_count_table.tsv` | Optional | KneadData QC read counts |
| `*humann_read_and_species_count_table.tsv` | Optional | HUMAnN alignment counts |
| `*humann_feature_counts.tsv` | Optional | HUMAnN feature detection counts |

**16S rRNA Amplicon** — produced by `biobakery_workflows 16s`:

| File (pattern) | Required | Description |
|----------------|----------|-------------|
| `*_otu_table.biom` or `*_otu_table.tsv` | ✓ | OTU table (QIIME1/USEARCH) |
| `*_asv_table.tsv` | ✓ | ASV table (DADA2) |

---

## Installation

### 1. Install Nextflow

```bash
curl -s https://get.nextflow.io | bash
chmod +x nextflow && sudo mv nextflow /usr/local/bin/
nextflow -version   # verify
```

### 2. Clone this repository

```bash
git clone https://github.com/biobakery/biobakery_workflows.git
cd biobakery_workflows/nextflow_vis
```

### 3. Choose your environment (pick one)

**Option A — Docker** (recommended for reproducibility):
```bash
docker pull biobakery/workflows:latest
# OR build the extended image (extra LaTeX/Pandoc packages)
docker build -t biobakery/workflows-vis:latest environments/
```

**Option B — Singularity** (recommended for HPC):
```bash
singularity pull biobakery_vis.sif docker://biobakery/workflows:latest
```

**Option C — Conda**:
```bash
conda env create -f environments/conda.yml
conda activate biobakery_vis
```

**Option D — Native** (already installed):
```bash
# Ensure anadama2 ≥ 0.10.0 (required for PweaveDocument.print_title)
pip install "anadama2>=0.10.0"
biobakery_workflows --help && pandoc --version && Rscript -e "library(vegan)"
```

---

## Quick Start

### WGX (metagenomics) report

```bash
nextflow run main.nf \
    -profile docker \
    --input  /data/project1/wmgx_output \
    --output results \
    --input-metadata /data/project1/metadata.tsv \
    --project-name "My Cohort Study" \
    --author-name  "Jane Smith" \
    --format pdf
```

### 16S amplicon report

```bash
nextflow run main.nf \
    -profile singularity \
    --input  /data/project1/16s_output \
    --output results \
    --input-metadata /data/project1/metadata.tsv \
    --format pdf
```

### With SLURM HPC

```bash
nextflow run main.nf \
    -profile slurm,singularity \
    --input  /data/project1/wmgx_output \
    --output results \
    --input-metadata /data/project1/metadata.tsv
```

---

## Sample Run — HMP2 IBD Cohort

The following is a **real run** performed on the Human Microbiome Project Phase 2 (HMP2)
tutorial dataset included in `examples/tutorial/stats_vis/input/`.
The dataset contains 165 samples from the Inflammatory Bowel Disease (IBD) cohort
with MetaPhlAn taxonomic profiles, HUMAnN pathway abundances, and KneadData QC counts.

### Command

```bash
# Equivalent Nextflow command:
nextflow run nextflow_vis/main.nf \
    --input  examples/tutorial/stats_vis/input/ \
    --output vis_results \
    --project-name  "HMP2 Tutorial — IBD Cohort" \
    --author-name   "bioBakery Team" \
    --introduction-text "Analysis of the HMP2 IBD cohort using bioBakery workflows (Nextflow edition)." \
    --format pdf \
    --min-abundance 0.01 \
    --min-samples 3 \
    --max-sets-heatmap 25 \
    --max-sets-barplot 15

# Underlying vis.py command executed by the Nextflow process:
biobakery_workflows vis \
    --input  examples/tutorial/stats_vis/input/ \
    --output vis_out \
    --project-name  "HMP2 Tutorial — IBD Cohort" \
    --author-name   "bioBakery Team" \
    --format pdf --local-jobs 1
```

### Run log

```
(Apr 22 06:22:30)  [0/2 -   0.00%] **Ready    ** Task 0: document
(Apr 22 06:22:30)  [0/2 -   0.00%] **Started  ** Task 0: document
(Apr 22 06:26:12)  [1/2 -  50.00%] **Completed** Task 0: document
(Apr 22 06:26:12)  [1/2 -  50.00%] **Ready    ** Task 2: archive
(Apr 22 06:26:12)  [1/2 -  50.00%] **Started  ** Task 2: archive
(Apr 22 06:26:12)  [2/2 - 100.00%] **Completed** Task 2: archive
Run Finished
```

### Run summary

| Item | Value |
|------|-------|
| Date | 2026-04-22 |
| Wall time | **3 min 46 sec** |
| Input samples | 165 (HMP2 IBD cohort) |
| Workflow type | WGX (metagenomics) |
| Report format | PDF |
| Report size | 1.1 MB |
| Archive size | 2.5 MB |
| Figures produced | **12 PNG files** |
| Data tables | 6 TSV files |
| Environment | Local, Python 3.8, anadama2 0.10.0 |

### Output PDF

The generated report is available at:
**[`sample_output_wmgx_report.pdf`](./sample_output_wmgx_report.pdf)**

Full archive (PDF + figures + data): **[`sample_output.zip`](./sample_output.zip)**

### Figures produced

| Figure | File | Description |
|--------|------|-------------|
| QC table (paired) | `qc_table_paired.png` | Per-sample paired-end read counts through KneadData filters |
| QC table (orphan) | `qc_table_orphan.png` | Per-sample orphan read counts |
| Microbial ratio table | `ratio_table.png` | Fraction of reads mapped to microbial vs. human reference |
| DNA read count barplot | `dna_read_count.png` | Paired reads per sample sorted by total count |
| DNA orphan barplot | `dna_read_count_orphan.png` | Orphan reads per sample |
| Taxa count table | `taxa_table.png` | Species/genera detected per sample (filtered vs. unfiltered) |
| PCoA — Species | `pcoa_species.png` | Bray-Curtis PCoA of all species |
| PCoA — Genera | `pcoa_genera.png` | Bray-Curtis PCoA of all genera |
| Species average abundance | `species_average_abundance.png` | Mean species abundance barplot |
| Genera average abundance | `genera_average_abundance.png` | Mean genera abundance barplot |
| Pathways average abundance | `pathways_average_abundance.png` | Mean pathway abundance barplot |
| Pathways table | `pathways_table.png` | Top pathways with average abundance and variance |

Sample figures are in [`sample_figures/`](./sample_figures/).

---

## Usage

```
nextflow run main.nf [options]

Required:
  --input <path>              Path to biobakery workflow output directory
                              (glob for multiple projects: "/data/*/wmgx_output")

Output:
  --output <path>             Output directory [default: vis_output]
  --format <pdf|html>         Report format [default: pdf]

Report metadata:
  --project-name <str>        Project name for report header
  --author-name <str>         Author name for report
  --header-image <path>       Image file for report header
  --introduction-text <str>   Custom introduction text

Input data options:
  --input-metadata <path>     Metadata TSV (samples as rows or columns)
  --input-file-type <str>     Override file type: 'filename,filetype' (repeat)
  --input-picard <path>       Folder of Picard quality score files

Metadata filtering:
  --metadata-categorical <f>  Treat feature as categorical (repeat)
  --metadata-continuous <f>   Treat feature as continuous (repeat)
  --metadata-exclude <f>      Exclude feature from analysis (repeat)
  --max-missing <float>       Max % missing values allowed [default: 20.0]

Abundance filtering:
  --min-abundance <float>     Min % abundance for filtering [default: 0.01]
  --min-samples <int>         Min % samples required [default: 3]

Visualization settings:
  --max-sets-heatmap <int>    Max features in heatmap [default: 25]
  --max-sets-barplot <int>    Max features in barplot [default: 15]
  --max-groups-barplot <int>  Max grouped barplots per variable [default: 5]
  --correlation-threshold <f> Min Spearman r for filtered heatmap [default: 0.7]

Nextflow options:
  -profile <name>             Execution profile (see Execution Environments)
  -resume                     Resume from last successful run
  -work-dir <path>            Nextflow work directory [default: ./work]

Profiles:
  standard          Local (default)       conda             Local + Conda
  docker            Local + Docker        slurm             SLURM + Singularity
  singularity       Local + Singularity   pbs               PBS + Singularity
                                          lsf               LSF + Singularity
                                          aws               AWS Batch + Docker
                                          gcp               GCP Batch + Docker
```

---

## Parameters

All `vis.py` parameters are supported 1-to-1:

| Nextflow param | Default | Description |
|----------------|---------|-------------|
| `--input` | required | Workflow output folder |
| `--output` | `vis_output` | Output directory |
| `--format` | `pdf` | `pdf` or `html` |
| `--project-name` | `""` | Report project name |
| `--author-name` | `""` | Report author |
| `--header-image` | `""` | Header image path |
| `--introduction-text` | `""` | Custom intro text |
| `--input-metadata` | `""` | Metadata TSV |
| `--input-file-type` | `[]` | Override file type |
| `--input-picard` | `""` | Picard QC folder |
| `--input-picard-extension` | `quality_by_cycle_metrics` | Picard extension |
| `--metadata-categorical` | `[]` | Categorical metadata features |
| `--metadata-continuous` | `[]` | Continuous metadata features |
| `--metadata-exclude` | `[]` | Excluded metadata features |
| `--max-missing` | `20.0` | Max % missing metadata |
| `--min-abundance` | `0.01` | Min % abundance |
| `--min-samples` | `3` | Min % samples |
| `--max-sets-heatmap` | `25` | Max features in heatmap |
| `--max-sets-barplot` | `15` | Max features in barplot |
| `--max-groups-barplot` | `5` | Max grouped barplots |
| `--correlation-threshold` | `0.7` | Spearman filter threshold |
| `--use-template` | `""` | Custom Pweave template |

---

## Execution Environments

### Local (no container)

```bash
# anadama2 ≥ 0.10.0 and biobakery_workflows must be installed locally
pip install "anadama2>=0.10.0"
nextflow run main.nf --input /path/to/output --output results
```

### Docker

```bash
nextflow run main.nf -profile docker \
    --input /path/to/output --output results
```

Override image:
```bash
nextflow run main.nf -profile docker \
    --container myregistry/my-biobakery:v4 --input /path/to/output
```

### Singularity / HPC

```bash
nextflow run main.nf -profile singularity \
    --input /path/to/output --output results

# Pre-pull and cache (avoids re-pulling on every job):
export NXF_SINGULARITY_CACHEDIR=/scratch/singularity_cache
nextflow run main.nf -profile singularity ...
```

### Conda

```bash
conda env create -f environments/conda.yml
nextflow run main.nf -profile conda \
    --input /path/to/output --output results
```

### SLURM

```bash
nextflow run main.nf \
    -profile slurm,singularity \
    --slurm_queue highmem \
    --slurm_account my_hpc_account \
    --input /path/to/output --output results
```

Add to `~/.nextflow/config` for site-wide SLURM defaults:
```groovy
params.slurm_queue   = 'short'
params.slurm_account = 'mylab'
```

### PBS / Torque

```bash
nextflow run main.nf -profile pbs,singularity \
    --pbs_queue batch --input /path/to/output
```

### IBM LSF

```bash
nextflow run main.nf -profile lsf,singularity \
    --lsf_queue normal --input /path/to/output
```

### AWS Batch

```bash
export AWS_ACCESS_KEY_ID=... AWS_SECRET_ACCESS_KEY=... AWS_DEFAULT_REGION=us-east-1
nextflow run main.nf \
    -profile aws \
    -work-dir s3://my-bucket/nextflow-work \
    --aws_batch_queue my-batch-queue \
    --input  s3://my-bucket/project1/wmgx_output \
    --output s3://my-bucket/project1/results
```

### Google Cloud Batch

```bash
gcloud auth application-default login
nextflow run main.nf \
    -profile gcp \
    -work-dir gs://my-bucket/nextflow-work \
    --gcp_project my-gcp-project \
    --input  gs://my-bucket/project1/wmgx_output \
    --output gs://my-bucket/project1/results
```

---

## Output Description

### Directory Structure

```
results/
└── <project_dir_name>/
    ├── vis_out/
    │   ├── wmgx_report.pdf         ← main report (or 16S_report.pdf)
    │   ├── figures/
    │   │   ├── qc_table_paired.png
    │   │   ├── qc_table_orphan.png
    │   │   ├── ratio_table.png
    │   │   ├── dna_read_count.png
    │   │   ├── dna_read_count_orphan.png
    │   │   ├── taxa_table.png
    │   │   ├── pcoa_species.png
    │   │   ├── pcoa_genera.png
    │   │   ├── heatmap_species_spearman.png     (WGX)
    │   │   ├── heatmap_species_bray_curtis.png  (WGX)
    │   │   ├── heatmap_species_zscore.png       (WGX)
    │   │   ├── heatmap_species_filtered_spearman.png (WGX)
    │   │   ├── heatmap_genus_spearman.png       (WGX)
    │   │   ├── heatmap_genus_bray_curtis.png    (WGX)
    │   │   ├── heatmap_genus_zscore.png         (WGX)
    │   │   ├── heatmap_genus_filtered_spearman.png  (WGX)
    │   │   ├── species_barchart_taxonomy.png    (WGX)
    │   │   ├── genera_barchart_taxonomy.png     (WGX)
    │   │   ├── species_average_abundance.png
    │   │   ├── genera_average_abundance.png
    │   │   ├── heatmap_pathways_spearman.png    (WGX with pathways)
    │   │   ├── heatmap_pathways_bray_curtis.png (WGX with pathways)
    │   │   ├── heatmap_pathways_zscore.png      (WGX with pathways)
    │   │   ├── pathways_average_abundance.png   (WGX with pathways)
    │   │   ├── pathways_table.png               (WGX with pathways)
    │   │   ├── heatmap_ecs_spearman.png         (WGX with ECs)
    │   │   └── alpha_diversity/                 (if metadata provided)
    │   │       ├── <var>_boxplot.png
    │   │       └── <var>_scatterplot.png
    │   └── data/
    │       ├── taxa_counts_table.tsv
    │       ├── qc_counts_pairs_table.tsv
    │       ├── qc_counts_orphans_table.tsv
    │       ├── microbial_counts_table.tsv
    │       ├── top_average_pathways_names.tsv
    │       └── alpha_diversity_with_metadata.txt (if metadata)
    └── vis_out.zip
pipeline_info/
    ├── report.html      ← Nextflow execution report
    ├── timeline.html    ← Nextflow timeline
    └── trace.tsv        ← per-task resource usage
```

### PDF Report Contents

| Section | Condition | Key figures |
|---------|-----------|-------------|
| Introduction | Always | Project overview, workflow diagram |
| Quality Control | KneadData files present | Read count tables, barplots, microbial ratio |
| Alpha Diversity | Metadata provided | Boxplots (categorical), scatterplots (continuous) |
| Taxonomic Count Table | WGX | Species/genera per sample before/after filter |
| PCoA Ordination | Always | Species PCoA, Genera PCoA, metadata-colored variants |
| Heatmaps — Taxa | Always | Spearman, Bray-Curtis, Z-score, Filtered Spearman |
| Stacked Barplots | Always | Species, genera, grouped by metadata |
| Functional Profiling | Pathway/EC files present | Pathway heatmaps, EC heatmaps, barplots, feature scatter plots |
| Pathway Table | Pathway files present | Top N pathways with abundance + variance |
| Workflow Info | Log file present | Tool versions, database versions |

---

## Monitoring (Web Dashboard)

The workflow includes a lightweight **real-time monitoring dashboard** for tracking
running jobs, resource usage, and pipeline progress.

### Start the dashboard

```bash
# Install monitoring dependencies (one-time)
pip install fastapi uvicorn psutil

# Launch the dashboard before or during a run
python monitoring/dashboard.py --port 8050

# Open in browser:
open http://localhost:8050
```

### Features

| Panel | Description |
|-------|-------------|
| **Pipeline status** | Real-time progress (tasks queued / running / completed / failed) |
| **Nextflow trace** | Per-task CPU %, memory, wall time from `trace.tsv` |
| **SLURM jobs** | Live view of `squeue` — job ID, state, node, user, time remaining |
| **Resource timeline** | CPU and memory over time for each task |
| **Log tail** | Last 50 lines of `.nextflow.log` with auto-refresh |

See [`monitoring/README.md`](./monitoring/README.md) for full setup instructions.

---

## Multi-Project Runs

```bash
# Process all project output directories simultaneously
nextflow run main.nf \
    -profile slurm,singularity \
    --input  "/data/cohorts/*/wmgx_output" \
    --output results

# Results are organized by project:
# results/project_A/vis_out/wmgx_report.pdf
# results/project_B/vis_out/wmgx_report.pdf
```

---

## Troubleshooting

### "No data files found in the input folder"

The `--input` directory must contain files matching bioBakery naming conventions.
Use `--input-file-type` for non-standard filenames:
```bash
--input-file-type "my_taxa.tsv,wmgx_taxonomy"
```

### "Cannot determine workflow type (16S or WGX)"

The workflow auto-detects type from file presence:
- `*_taxonomic_profiles.tsv` → WGX
- `*_otu_table.*` or `*_asv_table.tsv` → 16S

### "When running without a log file, please provide the introduction text"

```bash
--introduction-text "Analysis of 50 metagenomic samples from the XYZ cohort."
```

### `AttributeError: 'PweaveDocument' object has no attribute 'print_title'`

Upgrade anadama2 to ≥ 0.10.0:
```bash
pip install "anadama2>=0.10.0"
```

### PDF generation fails (LaTeX error)

- Use the Docker container (all LaTeX packages included).
- Try `--format html` to isolate whether the issue is data or LaTeX.
- Increase memory: `conf/base.config` → `memory = '64 GB'`.

### SLURM OOM (exit 137)

```groovy
// conf/slurm.config
withName: GENERATE_VIS_REPORT { memory = '64 GB' }
```

### Resume a partial run

```bash
nextflow run main.nf -resume [same parameters]
```

---

## Citation

If you use this workflow, please cite:

> McIver LJ, Abu-Ali G, Franzosa EA, et al. **bioBakery: a meta'omic analysis environment**. *Bioinformatics* 34(7):1235–1237, 2018. doi: 10.1093/bioinformatics/btx754

> Di Tommaso P, Chatzou M, Floden EW, et al. **Nextflow enables reproducible computational workflows**. *Nature Biotechnology* 35:316–319, 2017. doi: 10.1038/nbt.3820

---

## License

MIT License — same as bioBakery workflows. See the root [LICENSE](../LICENSE) file.
