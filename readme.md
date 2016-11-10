# bioBakery Workflows #

bioBakery Workflows is a collection of microbiome workflows. Each workflow is run with [AnADAMA2](https://bitbucket.org/biobakery/anadama2).

## Getting Started with bioBakery Workflows ##

### Requirements ###
1. [AnADAMA2](https://bitbucket.org/biobakery/anadama2) (installed automatically when installing bioBakery workflows)

### Installation ###

1. Download [biobakery_workflows.tar.gz](https://bitbucket.org/biobakery/biobakery_workflows/get/tip.tar.gz)

2. Install the software
    * ``$ tar zxvf biobakery_workflows.tar.gz``
    * ``$ cd biobakery_workflows``
    * ``$ python setup.py install``

## How to Run ##

### Basic Usage ###

``$ biobakery_workflows whole_metagenome_shotgun --input $INPUT_DIR --output $OUTPUT_DIR --kneadddata-db $DB``

This command will run the whole metagenome shotgun workflow on the files in the input folder ($INPUT_DIR to be replaced with the path to the folder containing fastq files). It will write files to the output folder ($OUTPUT_DIR to be replaced with the folder to write output files). It requires the location of the KneadData database ($DB to be replaced with the folder containing the KneadData database). 

For a list of all available workflows, run:

``$ biobakery_workflows --help``

For specific options for a workflow, run:

``$ biobakery_workflows whole_metagenome_shotgun --help``


