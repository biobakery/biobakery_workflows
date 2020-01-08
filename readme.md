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

