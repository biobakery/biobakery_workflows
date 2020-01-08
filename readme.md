bioBakery Workflows
===================

bioBakery workflows is a collection of workflows and tasks for executing
common microbial community analyses using standardized, validated tools
and parameters. Quality control and statistical summary reports are
automatically generated for most data types, which include 16S
amplicons, metagenomes, and metatranscriptomes. Workflows are run
directly from the command line and tasks can be imported to create your
own custom workflows. The workflows and tasks are built with
[AnADAMA2](https://bitbucket.org/biobakery/anadama2) which allows for
parallel task execution locally and in a grid compute environment.

For additional information, see the
[bioBakery workflows tutorial](https://bitbucket.org/biobakery/biobakery/wiki/biobakery_workflows).

Table of contents
  * [Getting Started](#getting-started)
    + [Requirements](#requirements)
    + [Installation](#installation)
      - [Install software](#install-software)
      - [Install databases](#install-databases)
    + [How to Run](#how-to-run)
      - [Basic Usage](#basic-usage)
      - [Data Processing Workflows](#data-processing-workflows)
      - [Visualization Workflows](#visualization-workflows)
    + [Parallelization Options](#parallelization-options)
  * [Data Processing Workflows](#data-processing-workflows-1)
    + [Whole Metagenome Shotgun (wmgx)](#whole-metagenome-shotgun-wmgx)
    + [Whole Metagenome and Metatranscriptome Shotgun (wmgx\_wmtx)](#whole-metagenome-and-metatranscriptome-shotgun-wmgx_wmtx)
    + [16S rRNA (16s)](#16s-rrna-16s)
    + [Isolate Assembly (isolate_assembly)](#isolate-assembly-isolate_assembly)
  * [Visualization Workflows](#visualization-workflows-1)
    + [Visualization for Whole Metagenome Shotgun (wmgx_vis)](#visualization-for-whole-metagenome-shotgun-wmgx_vis)
    + [Visualization for Whole Metagenome and Metatranscriptome Shotgun (wmgx_wmtx_vis)](#visualization-for-whole-metagenome-and-metatranscriptome-shotgun-wmgx_wmtx_vis)
    + [Visualization for 16S (16s_vis)](#visualization-for-16s-16s_vis)

------------------------------------------------------------------------

## Getting Started
------------------

### Requirements

1.  [AnADAMA2](https://bitbucket.org/biobakery/anadama2) (installed
    automatically)
2.  Python v2.7+
3.  See individual workflows and tasks for additional software
    requirements.

### Installation

#### Install software

bioBakery workflows can be installed with Conda, Docker, or pip.

To install with Conda:

    $ conda install -c biobakery biobakery_workflows
    
-   Before installing the tool, configure your channels so biobakery is
    at the top of the list.
-   This will install almost all of the dependencies for all workflows
    (ie Kneaddata, MetaPhlan2, etc.) excluding those dependencies that
    have licenses.

To install and run with Docker:

    $ docker run -it biobakery/workflows bash

-   The image will include all dependencies for all workflows (ie
    Kneaddata, MetaPhlan2, etc.) excluding those dependencies that have
    licenses.

To install with pip:

    $ pip install biobakery_workflows

-   This will only install the core software and dependencies. It will
    not install the dependencies for the workflows (ie KneadData,
    MetaPhlan2, etc.).

#### Install databases

*Install automatically*

Once the software and dependencies are installed, the databases can be
installed automatically.

Run the following command to install the databases required for a
workflow:

`$ biobakery_workflows_databases --install $WORKFLOW`

-   Replace `$WORKFLOW` with the workflow name (ie wmgx, 16s,
    wmgx\_wmtx, or wmgx\_demo, isolate\_assembly).
-   The databases will be installed at
    `$HOME/biobakery_workflow_databases/` or
    `/opt/biobakery_workflow_databases/` depending on permissions.
-   To install to a custom location add the option `--location $FOLDER`.
    With this option you will also need to set the environment variable
    `$BIOBAKERY_WORKFLOWS_DATABASES` to the folder so the workflows can
    find the installed databases.
-   The database install requires some of the dependencies from the
    corresponding workflow to build and install the databases. For
    example, installing the wmgx databases requires HUMAnN2, KneadData,
    StrainPhlAn and bowtie2. Please install these dependencies prior to
    installing the databases. Depending on the method used to install
    the workflows you might need to install these dependencies in
    addition to the workflow.

*Install manually*

Alternatively the databases can be installed manually and then
referenced with environment variables. The shotgun data processing
workflows require Kneaddata (human, human transcriptome, and SILVA),
HUMAnN2 (utility mapping, nucleotide, and protein databases), and
StrainPhlAn (reference and marker) databases while the 16s data
processing workflow requires the GreenGenes fasta, taxonomy, and usearch
formatted files.

When manually installing the databases, the following environment
variables need to be set.

-   Shotgun workflows: `KNEADDATA_DB_HUMAN_GENOME`,
    `KNEADDATA_DB_RIBOSOMAL_RNA`, `KNEADDATA_DB_HUMAN_TRANSCRIPTOME`,
    `STRAINPHLAN_DB_REFERENCE`, and `STRAINPHLAN_DB_MARKERS`.
-   16s workflows: `GREEN_GENES_USEARCH_DB`, `GREEN_GENES_FASTA_DB`, and
    `GREEN_GENES_TAXONOMY_DB`.
    
### How to Run

#### Basic Usage

All workflows follow the general command format:

`$ biobakery_workflows $WORKFLOW --input $INPUT --output $OUTPUT`

For a list of all available workflows, run:

`$ biobakery_workflows --help`

For specific options for a workflow, run:

`$ biobakery_workflows $WORKFLOW --help`

#### Data Processing Workflows

The basic command to run a data processing workflow, replacing
`$WORKFLOW` with the workflow name, is:

`$ biobakery_workflows $WORKFLOW --input $INPUT_DIR --output $DATA_OUTPUT_DIR `

This command will run the workflow on the files in the input folder
(`$INPUT_DIR` to be replaced with the path to the folder containing
fastq files). It will write files to the output folder
(`$DATA_OUTPUT_DIR` to be replaced with the folder to write output
files).

#### Visualization Workflows

A visualization workflow exists corresponding to each data processing
workflow. The basic command to run a visualization workflow, replacing
`$WORKFLOW_VIS` with the visualization workflow name, is:

`$ biobakery_workflows $WORKFLOW_VIS --input $DATA_OUTPUT_DIR --project-name $PROJECT --output $OUTPUT_DIR `

The input folder (`$DATA_OUTPUT_DIR` to be replaced with the path to the
folder) in this command is the output folder from the data processing
workflow. The folder (`$OUTPUT_DIR` to be replaced with the path to the
output folder) will contain the output files from the visualization
workflow. The project name should replace `$PROJECT` in the command so
the report can include the name.

### Parallelization Options

When running any workflow you can add the following command line options
to make use of existing computing resources:

-   `--local-jobs <1>` : Run multiple tasks locally in parallel. Provide
    the max number of tasks to run at once. The default is one task
    running at a time.
-   `--grid-jobs <0>` : Run multiple tasks on a grid in parallel.
    Provide the max number of grid jobs to run at once. The default is
    zero tasks are submitted to a grid resulting in all tasks running
    locally.
-   `--grid <slurm>` : Set the grid available on your machine. This will
    default to the grid found on the machine with options of slurm and
    sge.
-   `--partition <serial_requeue>` : Jobs will be submitted to the
    partition selected. The default partition selected is based on the
    default grid.

For additional workflow options, see the
[AnADAMA2](https://bitbucket.org/biobakery/anadama2) user manual.

------------------------------------------------------------------------

## Data Processing Workflows
-------------------------

bioBakery workflows includes a collection of workflows for shotgun
sequences and 16s data processing. Most workflows can be run on the
command line with the following syntax:

`$ biobakery_workflows $WORKFLOW --input $INPUT --output $OUTPUT`

See the section on
[parallelization options](https://bitbucket.org/biobakery/biobakery_workflows/wiki/Home#!parallelization-options)
to optimize the workflow run based on your computing resources.

### Whole Metagenome Shotgun (wmgx)

![](https://github.com/biobakery/biobakery_workflows/blob/master/images/135637619-111018-wms_workflow.png)

**Super Tasks**

1.  [Quality control](https://bitbucket.org/biobakery/biobakery_workflows/wiki/Home#!shotgun-quality-control)
2.  [Taxonomic profiling](https://bitbucket.org/biobakery/biobakery_workflows/wiki/Home#!shotgun-taxonomic-profiling)
3.  [Functional profiling](https://bitbucket.org/biobakery/biobakery_workflows/wiki/Home#!shotgun-functional-profiling)
4.  Strain profiling
5.  Assembly (not run by default)

**Requirements**

1.  [KneadData](https://bitbucket.org/biobakery/kneaddata) (v0.7.0+)
    1.  Install with: `$ conda install -c biobakery kneaddata` OR
        `$ pip install kneaddata`
2.  [MetaPhlAn2](https://bitbucket.org/biobakery/metaphlan2)
    1.  Install with: `$ conda install -c bioconda metaphlan2`
3.  [HUMAnN2](https://bitbucket.org/biobakery/humann2)
    1.  Install with: `$ conda install -c biobakery humann2` OR
        `$ pip install humann2`
4.  [StrainPhlAN](http://segatalab.cibio.unitn.it/tools/strainphlan/)
    1.  Install with: `$ conda install -c bioconda strainphlan`
5.  [Prokka](https://github.com/tseemann/prokka) (Only required if
    running in assembly mode)
6.  [MegaHit](https://github.com/voutcn/megahit) (Only required if
    running in assembly mode)
7.  [seqtk](https://github.com/lh3/seqtk) (Only required if running in
    assembly mode)

**Inputs**

1.  A set of fastq (or fastq.gz) files (single-end or paired-end). The
    files are expected to be named
    `$SAMPLE.fastq.gz`,`$SAMPLE.R1.fastq.gz`, or `$SAMPLE.R2.fastq.gz`
    where `$SAMPLE` is the sample name or identifier corresponding to
    the sequences. `$SAMPLE` can contain any characters except spaces or
    periods.

The workflow will detect if paired-end files are present. By default the
workflow identifies paired end reads based on file names containing
".R1" and ".R2" strings. If your paired end reads have different
identifiers, use the option `--pair-identifier .R1` to provide the
identifier string for the first file in the set of pairs.

The workflow by default expects input files with the extension
"fastq.gz". If your files are not gzipped, run with the option
`--input-extension fastq`.

**To run the workflow**

-   `$ biobakery_workflows wmgx --input $INPUT --output $OUTPUT `
-   In the command replace `$INPUT` with the path to the folder
    containing your fastq input files and `$OUTPUT` with the path to the
    folder to write output files.
-   See the section on
    [parallelization options](https://bitbucket.org/biobakery/biobakery_workflows/wiki/Home#!parallelization-options)
    to optimize the workflow run based on your computing resources.
-   The workflow runs with the default settings for all main tool
    subtasks. These settings will work for most data sets. However, if
    you need to customize your workflow settings for the KneadData and
    StrainPhlAn subtasks please read the documentation for each tool to
    determine the optimum settings. Then apply these settings by using
    options for each tool. For example, `--qc-options="$OPTIONS"` will
    modify the default settings when running the KneadData subtask and
    `--strain-profiling-options="$OPTIONS"` will modify the options when
    running the StrainPhlAn subtask (replacing the `$OPTIONS` in each
    with your selected settings).
-   Add the option `--run-assembly` to add the tasks to run assembly.

**To run a demo**

-   Single-end
    -   `$ biobakery_workflows wmgx --input examples/wmgx/single/ --output workflow_output `
-   Paired-end
    -   `$ biobakery_workflows wmgx --input examples/wmgx/paired/ --output workflow_output `
-   Demo input files can be found in the biobakery\_workflow source
    `examples` folder.

### Whole Metagenome and Metatranscriptome Shotgun (wmgx\_wmtx)

![](https://github.com/biobakery/biobakery_workflows/blob/master/images/4071512812-wms_wtx_workflow.jpg)

**Super Tasks**

1.  [Quality control](https://bitbucket.org/biobakery/biobakery_workflows/wiki/Home#!shotgun-quality-control)
2.  [Taxonomic profiling](https://bitbucket.org/biobakery/biobakery_workflows/wiki/Home#!shotgun-taxonomic-profiling)
3.  [Functional profiling](https://bitbucket.org/biobakery/biobakery_workflows/wiki/Home#!shotgun-functional-profiling)
4.  Strain profiling

**Requirements**

1.  [KneadData](https://bitbucket.org/biobakery/kneaddata) (v0.7.0+)
    1.  Install with: `$ conda install -c biobakery kneaddata` OR
        `$ pip install kneaddata`
2.  [MetaPhlAn2](https://bitbucket.org/biobakery/metaphlan2)
    1.  Install with: `$ conda install -c bioconda metaphlan2`
3.  [HUMAnN2](https://bitbucket.org/biobakery/humann2)
    1.  Install with: `$ conda install -c biobakery humann2` OR
        `$ pip install humann2`
4.  [StrainPhlAN](http://segatalab.cibio.unitn.it/tools/strainphlan/)
    1.  Install with: `$ conda install -c bioconda strainphlan`

**Inputs**

1.  Two sets of fastq (or fastq.gz) files (single-end or paired-end).
    One set is of whole metagenome shotgun data and the other is whole
    metatranscriptome shotgun data. The files are expected to be named
    `$SAMPLE.fastq.gz`,`$SAMPLE.R1.fastq.gz`, or `$SAMPLE.R2.fastq.gz`
    where `$SAMPLE` is the sample name or identifier corresponding to
    the sequences. `$SAMPLE` can contain any characters except spaces or
    periods.
2.  Optionally, provide a mapping file. This file will have two columns
    and be tab delimited. The first column is the sample names for the
    metatranscriptomes and the second is the corresponding metagenome
    sample. See the demo mapping file for an example.

The workflow will detect if paired-end files are present. By default the
workflow identifies paired end reads based on file names containing
".R1" and ".R2" strings. If your paired end reads have different
identifiers, use the option `--pair-identifier .R1` to provide the
identifier string for the first file in the set of pairs.

The workflow by default expects input files with the extension
"fastq.gz". If your files are not gzipped, run with the option
`--input-extension fastq`.

**To run the workflow**

-   `$ biobakery_workflows wmgx_wmtx --input-metagenome $INPUT_WMS --input-metatranscriptome $INPUT_WTS --input-mapping $INPUT_MAPPING --output $OUTPUT `
-   In the command replace `$INPUT_WMS` with the path to the folder
    containing your whole metagenome shotgun fastq.gz input files,
    `$INPUT_WTS` with the path to the folder containing your whole
    metatranscriptome shotgun fastq.gz input files, and `$OUTPUT` with
    the path to the folder to write output files. Replace
    `$INPUT_MAPPING` with your file of mapping between the metagenome
    and metatranscriptome samples.
-   See the section on
    [parallelization options](https://bitbucket.org/biobakery/biobakery_workflows/wiki/Home#!parallelization-options)
    to optimize the workflow run based on your computing resources.
-   The workflow runs with the default settings for all main tool
    subtasks. These settings will work for most data sets. However, if
    you need to customize your workflow settings for the KneadData
    subtasks please read the documentation for KneadData to determine
    the optimum settings. Then apply these settings by using the option
    `--qc-options="$OPTIONS"`.

**To run a demo**

-   `$ biobakery_workflows wmgx_wmtx --input-metagenome examples/wmgx_wmtx/wms/ --input-metatranscriptome examples/wmgx_wmtx/wts/ --input-mapping examples/wmgx_wmtx/mapping.tsv --output workflow_output `
-   Demo input files can be found in the biobakery\_workflow source
    `examples` folder.

### 16S rRNA (16s)

The 16s workflow has two methods that can be used: UPARSE (with either
USEARCH or VSEARCH (default)) and DADA2. All methods perform quality
control and generate taxonomic tables.

**Workflow diagrams**

![](https://github.com/biobakery/biobakery_workflows/blob/master/images/1618113365-16s_workflow_with_fasttree_merge_truncate_filter.jpg)
 
![](https://github.com/biobakery/biobakery_workflows/blob/master/images/4005903235-DADA2Workflow_75percent.jpg)
 
**Super Tasks**

1.  Demultiplex (only for raw reads)
2.  Merge Samples and Rename
3.  Quality Control (Learn Error Rates for DADA2 method)
4.  Taxonomic Profiling
5.  Functional Profiling (with PICRUSt v1 or v2; v2 requires python3)

**Requirements**

1.  [Vsearch](https://github.com/torognes/vsearch)
2.  [Usearch](http://www.drive5.com/usearch/) (v9.0.2132) (Not required
    for DADA2 or VSEARCH method)
3.  [PICRUSt v1.1 (or v2)](http://picrust.github.io/picrust/)
    1.  Install with: `$ conda install -c bioconda picrust`
    2.  OR Install with: `$ conda install -c bioconda picrust2`
4.  [BIOM v2](http://biom-format.org/) (Not required for DADA2 method)
    1.  Install with: `$ pip install biom-format`
5.  [Clustal Omega](http://www.ebi.ac.uk/Tools/msa/clustalo/)
    1.  Install with: `$ conda install -c bioconda clustal-omega`
6.  [EA-Utils](https://expressionanalysis.github.io/ea-utils/) (Only
    required for raw reads)
    1.  Install with: `$ conda install -c bioconda ea-utils`
7.  [FastTree](http://www.microbesonline.org/fasttree/)
    1.  Install with: `$ conda install -c bioconda fasttree`
8.  R and the following packages: dada2, gridExtra, tools, ggplot2,
    seqinr (Only required for DADA2 method)
    1.  Follow the
        [DADA2 install directions](https://benjjneb.github.io/dada2/dada-installation.html)

**Inputs**

1.  A set of fastq (or fastq.gz) files (single-end or paired-end; only
    pair-end for DADA2) raw or demultiplexed. If the files are
    demultiplexed, the files are expected to be named
    `$SAMPLE.fastq.gz`,`$SAMPLE_R1_001.fastq.gz`, or
    `$SAMPLE_R2_001.fastq.gz` where `$SAMPLE` is the sample name or
    identifier corresponding to the sequences. `$SAMPLE` can contain any
    characters except spaces or periods.
2.  A barcode file (only required for raw reads).

The workflow will detect if paired-end files are present. By default the
workflow identifies paired end reads based on file names containing
"\_R1\_001" and "\_R2\_001" strings. If your paired end reads have
different identifiers, use the option `--pair-identifier .R1.` to
provide the identifier string for the first file in the set of pairs.

The workflow by default expects input files with the extension
"fastq.gz". If your files are not gzipped, run with the option
`--input-extension fastq`.

**To run the workflow**

-   `$ biobakery_workflows 16s --input $INPUT --output $OUTPUT `
-   In the command replace `$INPUT` with the path to the folder
    containing your fastq input files and `$OUTPUT` with the path to the
    folder to write output files.
-   By default the workflow implements the UPARSE workflow using VSEARCH
    plus additional tasks. To run with USEARCH as the executable add the
    option `--method dada2`.
-   See the section on
    [parallelization options](https://bitbucket.org/biobakery/biobakery_workflows/wiki/Home#!parallelization-options)
    to optimize the workflow run based on your computing resources.
-   The workflow runs with default settings for the subtasks. Depending
    on the lengths of your reads you might want to change the setting
    `--trunc-len-max 200`, if running the VSEARCH/USEARCH method, to a
    smaller value. Reading through the maxee table will help to
    determine the length to use for trimming based on the joined reads
    and their quality scores. For other default settings, please run the
    workflow with the `--help` option. All of the other settings will
    work for most data sets. If there are any you would like to change,
    please review the usearch documentation to determine the optimal
    settings.
-   Add the option `--method dada2` to run the DADA2 method instead of
    VSEARCH.

### Isolate Assembly (isolate_assembly)

This workflow will assemble and annotate sequenced microbial isolate
genomes. It runs the raw sequences through quality control, assembly
(with SPAdes), annotation (with Prokka), functional annotation, quality
assessment, and then creates a final annotated contig file.

**Workflow diagram**

![](https://github.com/biobakery/biobakery_workflows/blob/master/images/2782656604-isolate_workflow.jpg)

**Super Tasks**

1.  Quality control
2.  Assembly
3.  Annotation
4.  Functional annotation
5.  Quality assessment

**Requirements**

1.  [KneadData](https://bitbucket.org/biobakery/kneaddata) (v0.7.0+)
    1.  Install with: `$ conda install -c biobakery kneaddata` OR
        `$ pip install kneaddata`
2.  [SPAdes](http://cab.spbu.ru/software/spades/)
    1.  Install with: `$ conda install -c bioconda spades`
3.  [Prokka](https://github.com/tseemann/prokka)
    1.  Install with: `$ conda install -c bioconda prokka`
4.  [Quast](http://quast.sourceforge.net/quast)
    1.  Install with: `$ conda install -c bioconda quast`
5.  [CheckM](https://ecogenomics.github.io/CheckM/) (plus databases)
6.  [run\_DBcan](https://github.com/linnabrown/run_dbcan) (plus
    databases) 
7.  [Emapper](https://github.com/jhcepas/eggnog-mapper) (version 2+)
    (databases installed with biobakery utility,
    biobakery\_workflows\_databases)

**Inputs**

1.  A set of fastq (or fastq.gz) files (single-end or paired-end) raw.
    This workflow does not allow for multiplexed files as input. The
    files are expected to be named
    `$SAMPLE.fastq.gz`,`$SAMPLE_R1_001.fastq.gz`, or
    `$SAMPLE_R2_001.fastq.gz` where `$SAMPLE` is the sample name or
    identifier corresponding to the sequences. `$SAMPLE` can contain any
    characters except spaces or periods.

The workflow will detect if paired-end files are present. By default the
workflow identifies paired end reads based on file names containing
"\_R1\_001" and "\_R2\_001" strings. If your paired end reads have
different identifiers, use the option `--pair-identifier .R1.` to
provide the identifier string for the first file in the set of pairs.

**To run the workflow**

-   `$ biobakery_workflows isolate_assembly --input $INPUT --species-name $SPECIES --output $OUTPUT `
-   In the command replace `$INPUT` with the path to the folder
    containing your fastq input files, `$SPECIES` with the name of the
    isolate sequenced and `$OUTPUT` with the path to the folder to write
    output files. The `$SPECIES` input string is used as the basename of
    the contig files.
-   See the section on
    [parallelization options](https://bitbucket.org/biobakery/biobakery_workflows/wiki/Home#!parallelization-options)
    to optimize the workflow run based on your computing resources.

Visualization Workflows
-----------------------

bioBakery workflows includes a collection of visualization workflows for
shotgun sequences and 16s data. Most workflows can be run on the command
line with the following syntax:

`$ biobakery_workflows $WORKFLOW --input $INPUT --project-name $PROJECT --output $OUTPUT`

The `$OUTPUT` folder of a data processing workflow can be used as the
`$INPUT` folder to the corresponding visualization workflow. For
detailed information on the input files required for each visualization
workflow, see the help message for the workflow by running the command:

`$ biobakery_workflows $WORKFLOW --help`

### Visualization for Whole Metagenome Shotgun (wmgx_vis)

This workflow generates a document of tables, bar plots, a PCoA plot,
scatter plots, and heatmaps using the output of the wmgx workflow as
input.

**Requirements**

1.  [Pweave](http://mpastell.com/pweave) (installed automatically)
2.  [NumPy and SciPy](https://docs.scipy.org)
    1.  Install with: `$ pip install numpy ` AND `$ pip install scipy `
3.  [Matplotlib](http://matplotlib.org/)
    1.  Install with: `$ pip install matplotlib `
4.  [LaTeX](https://www.latex-project.org/)
5.  [Pandoc](http://pandoc.org/) (&lt;version2 required)
    1.  Install with: `$ conda install pandoc `
6.  [Hclust2](https://bitbucket.org/nsegata/hclust2)
    1.  Install with: `$ conda install -c biobakery hclust2 `
7.  [R](https://www.r-project.org/) with the
    [vegan](https://cran.r-project.org/web/packages/vegan/index.html)
    package

**Inputs**

1.  An input folder containing the final products from the wmgx data
    workflow.
    1.  A file of the KneadData read counts for the wmgx samples (single
        or paired end).
    2.  A file of the merged taxonomic profile.
    3.  A file of the merged pathway abundances (normalized).
    4.  A file of the HUMAnN2 alignment counts.
    5.  A file of the HUMAnN2 feature counts.
    6.  The log file from the corresponding data processing workflow.
2.  The project name.
3.  Introduction text (Optional).
4.  The report format (Options: pdf/html, pdf is default).

**Outputs**

1.  A pdf (or html) report.
2.  A directory of figures included in the report.
    1.  Quality control section
        1.  Paired end read count table and barchart
        2.  Orphan read count table and barchart
        3.  Microbial read proportion table
    2.  Taxonomy section
        1.  Species count table
        2.  Ordination (PCoA)
        3.  Heatmap of top species abundances
        4.  Stacked barplot of top species abundances
    3.  Functional profiling section
        1.  Pathway abundance
            1.  Top pathways by abundance (heatmap and table)
            2.  Top pathways by variance (heatmap and table) \#\#
                Feature counts section
        2.  Scatter plot of the aligned reads (nucleotide/translated)
        3.  Scatter plot of the gene families counts
            (nucleotide/translated)
        4.  Scatter plot of the ECs counts (nucleotide/translated)
        5.  Scatter plot of the pathway counts (nucleotide/translated)
3.  A directory of read count tables
4.  A zip archive containing the pdf (or html) report, figures, and data
    files

**To run the workflow**

-   `$ biobakery_workflows wmgx_vis --input $INPUT --project-name $PROJECT --output $OUTPUT `
-   In the command replace `$INPUT` with the output folder created by
    running the wmgx data processing workflow, `$PROJECT` with the name
    of the project, and `$OUTPUT` with the path to the folder to write
    output files.

### Visualization for Whole Metagenome and Metatranscriptome Shotgun (wmgx_wmtx_vis)

This workflow generates a document of tables, bar plots, a PCoA plot,
and heatmaps using the output of the wmgx\_wmtx workflow as input.

**Requirements**

1.  [Pweave](http://mpastell.com/pweave) (installed automatically)
2.  [NumPy and SciPy](https://docs.scipy.org)
    1.  Install with: `$ pip install numpy ` AND `$ pip install scipy `
3.  [Matplotlib](http://matplotlib.org/)
    1.  Install with: `$ pip install matplotlib `
4.  [LaTeX](https://www.latex-project.org/)
5.  [Pandoc](http://pandoc.org/) (&lt;version2 required)
    1.  Install with: `$ conda install pandoc `
6.  [Hclust2](https://bitbucket.org/nsegata/hclust2)
    1.  Install with: `$ conda install -c biobakery hclust2 `
7.  [R](https://www.r-project.org/) with the
    [vegan](https://cran.r-project.org/web/packages/vegan/index.html)
    package

**Inputs**

Please note this workflow is currently only for paired end reads.

1.  An input folder containing the final products from the wmgx\_wmtx
    data workflow.
    1.  A file of the KneadData read counts for the wmgx samples (paired
        end).
    2.  A file of the KneadData read counts for the wmtx samples (paired
        end).
    3.  A file of the merged taxonomic profile for the wmgx samples.
    4.  A file of the merged pathway abundances for the wmgx samples
        (normalized).
    5.  The log file from the corresponding data processing workflow.
    6.  A file of the HUMAnN2 alignment counts for the wmgx samples
        (Optional).
    7.  A file of the HUMAnN2 feature counts for the wmgx samples
        (Optional).
    8.  A file of the HUMAnN2 alignment counts for the wmtx samples
        (Optional).
    9.  A file of the HUMAnN2 feature counts for the wmtx samples
        (Optional).
    10. A file of the gene families RNA/DNA normalized (Optional).
    11. A file of the ECs RNA/DNA normalized (Optional).
    12. A file of the pathway abundances RNA/DNA normalized (Optional).
2.  The project name.
3.  Introduction text (Optional).
4.  The report format (Options: pdf/html, pdf is default).

**Outputs**

1.  A pdf (or html) report.
2.  A directory of figures included in the report.
    1.  Quality control section
        1.  Paired end read count table and barchart (for DNA and RNA)
        2.  Orphan read count table and barchart (for DNA and RNA)
        3.  Microbial read proportion table (for DNA and RNA)
    2.  Taxonomy section
        1.  Species count table
        2.  Ordination (PCoA)
        3.  Heatmap of top species abundances
        4.  Stacked barplot of top species abundances
    3.  Functional profiling section
        1.  Pathway abundance
            1.  Top pathways by abundance (heatmap and table)
            2.  Top pathways by variance (heatmap and table)
        2.  Feature counts section (for DNA and RNA, optional)
            1.  Scatter plot of the aligned reads
                (nucleotide/translated)
            2.  Scatter plot of the gene families counts
                (nucleotide/translated)
            3.  Scatter plot of the ECs counts (nucleotide/translated)
            4.  Scatter plot of the pathway counts
                (nucleotide/translated)
3.  A zip archive containing the pdf (or html) report, figures, and data
    files

### Visualization for 16S (16s_vis)

This workflow generates a document of bar plots and a PCoA plot using
the output of the 16S workflow as input.

**Requirements**

1.  [Pweave](http://mpastell.com/pweave) (installed automatically)
2.  [NumPy and SciPy](https://docs.scipy.org)
    1.  Install with: `$ pip install numpy ` AND `$ pip install scipy `
3.  [Matplotlib](http://matplotlib.org/)
    1.  Install with: `$ pip install matplotlib `
4.  [LaTeX](https://www.latex-project.org/)
5.  [Pandoc](http://pandoc.org/) (&lt;version2 required)
    1.  Install with: `$ conda install pandoc `
6.  [Hclust2](https://bitbucket.org/nsegata/hclust2)
    1.  Install with: `$ conda install -c biobakery hclust2 `
7.  [R](https://www.r-project.org/) with the
    [vegan](https://cran.r-project.org/web/packages/vegan/index.html)
    package

**Inputs**

1.  An input folder containing the final products from the 16s data
    workflow.
    1.  A file of the closed reference OTU table.
    2.  A file of the read counts per sample (including total reads,
        classified, and unclassified).
    3.  A file of the eestats for all samples.
    4.  The log file from the corresponding data processing workflow.
2.  The project name.
3.  Introduction text (Optional).
4.  The report format (Options: pdf/html, pdf is default)

**Outputs**

1.  A pdf (or html) report.
2.  A directory of figures included in the report.
    1.  Read counts per sample barchart
    2.  OTU counts per sample barchart
    3.  Stacked barchart of the top 15 genera by average abundance
    4.  Stacked barchart of the top 15 terminal taxa by average
        abundance
    5.  Ordination (PCoA) of terminal taxa
3.  A zip archive containing the pdf (or html) report, figures, and data
    files

**To run the workflow**

-   `$ biobakery_workflows 16s_vis --input $INPUT --project-name $PROJECT --output $OUTPUT `
-   In the command replace `$INPUT` with the output folder created by
    running the 16s data processing workflow, `$PROJECT` with the name
    of the project, and `$OUTPUT` with the path to the folder to write
    output files.


