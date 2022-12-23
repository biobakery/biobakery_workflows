# assembly_workflow

> This repository contains code for Will's assembly workflow based on Aaron Walsh's pets assembly workflow and Segata et al. 2019. For queries, contact Will Nickols (email: <willnickols@college.harvard.edu>).

# Installation

The workflow can be installed with the following commands.  It is okay if the Checkm2 database install throws an error starting with `File "checkm2/bin/checkm2", line 244, in <module>...`.  The PhyloPhlAn run will intentionally fail in order to download the database but will take a while to download the database.
```
git clone https://github.com/WillNickols/assembly_workflow
cd assembly_workflow
conda env create -f assembly_environment.yml
conda activate biobakery_assembly
git clone --recursive https://github.com/chklovski/checkm2.git
checkm2/bin/checkm2 database --download --path databases/checkm/
export CHECKM_DATA_PATH=$(pwd)/databases/checkm/
mkdir tmp_to_delete && mkdir -p databases/phylophlan && touch tmp_to_delete/tmp.fa && phylophlan_metagenomic -d SGB.Jul20 --database_folder databases/phylophlan/ -i tmp_to_delete/; rm -r tmp_to_delete/
export PHYLOPHLAN_PATH=$(pwd)/databases/phylophlan/
```

Run the following commands to install the necessary R packages.
```
R
install.packages(c("docopt", "dplyr", "data.table", "stringr", "doParallel", "tidyr"))
q()
```

Once the conda environment is created and you are in the `assembly_workflow` directory, you can activate the environment with these commands:
```
conda activate biobakery_assembly
export CHECKM_DATA_PATH=$(pwd)/databases/checkm/
export PHYLOPHLAN_PATH=$(pwd)/databases/phylophlan/
```

# Example runs

This runs the workflow on a set of paired-end files that have already been cleaned with Kneaddata.  The `by_sample` flag specifies that abundance estimates of the sample's taxa should come from aligning that sample's reads only to its contigs.  Alternatively, the `by_dataset` flag can be used to align a sample's reads against all contigs generated from the input folder when the samples likely share the same taxa.
```
python assembly_workflow.py \
  -i example/input/sample_1/ \
  -o example/output/sample_1/ \
  --abundance-type by_sample \
  --input-extension fastq.gz \
  --paired paired \
  --pair-identifier _R1 \
  --cores 8 \
  --local-jobs 12 \
  --remove-intermediate-output
```

The output, `example/output/sample_1/final_profile_by_sample.tsv` is a MetaPhlAn-like output with the abundance of each taxonomic group identified in the sample.  The output shows that Pseudoalteromonas marina is present in the sample along with a species, `sgb_01` that represents a new species genome bin not within a Mash distance of 0.05 of any known SGB, 0.15 of any known GGB, or 0.3 of any known FGB.

This runs the workflow on a concatenated file with pre-created contigs in `example/output/sample_2/assembly/main/sample_2/sample_2.contigs.fa`.  Such a file could be created from the `biobakery_workflows wmgx` workflow with `--run-assembly`.  Note that the original reads are still required for the abundance calculation.  Additionally, grid options and scratch space information can be provided as seen below.
```
python assembly_workflow.py \
  -i example/input/sample_2/ \
  -o example/output/sample_2/ \
  --abundance-type by_sample \
  --input-extension fastq.gz \
  --paired concatenated \
  --skip-contigs \
  --cores 8 \
  --local-jobs 12 \
  --grid-scratch /n/holyscratch01/nguyen_lab/wnickols/assembly/example/sample_2/ \
  --grid-partition 'shared' --grid-jobs 8 --cores 8 --mem 20000
```

# Extra checks

This command runs the workflow without taxonomically placing the MAGs (it runs only assembly, binning, and quality checking).
```
python assembly_workflow.py \
  -i /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/mags_and_sgbs_pipeline_testing/test_inputs/single_end/ \
  -o /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/mags_and_sgbs_pipeline_testing/test_outputs/single_end/ \
  --abundance-type by_sample --input-extension fastq.gz --paired unpaired \
  --local-jobs 12 \
  --skip-placement \
  --remove-intermediate-output
```

This command runs a single-end `fastq.gz` file.
```
python assembly_workflow.py \
  -i /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/mags_and_sgbs_pipeline_testing/test_inputs/single_end/ \
  -o /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/mags_and_sgbs_pipeline_testing/test_outputs/single_end/ \
  --abundance-type by_sample --input-extension fastq.gz --paired unpaired \
  --local-jobs 12 \
  --remove-intermediate-output
```

This command runs a paired-end `fastq` file.  Read headers should end with "/1" or "/2" if the files are paired (e.g. `@read_57/1` and `@read_57/2`).
```
python assembly_workflow.py \
  -i /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/mags_and_sgbs_pipeline_testing/test_inputs/paired_end/ \
  -o /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/mags_and_sgbs_pipeline_testing/test_outputs/paired_end/ \
  --abundance-type by_sample --input-extension fastq --paired paired \
  --local-jobs 12 \
  --remove-intermediate-output \
  --cores 8
```

This command runs two concatenated `fastq.gz` files, one of which is single-end and one of which is paired-end.  These read headers should also end with "/1" and "/2" to indicate pairing.  Files from Kneaddata automatically satisfy this requirement.
```
python assembly_workflow.py \
  -i /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/mags_and_sgbs_pipeline_testing/test_inputs/concat/ \
  -o /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/mags_and_sgbs_pipeline_testing/test_outputs/concat/ \
  --abundance-type by_sample \
  --input-extension fastq.gz \
  --paired concatenated \
  --cores 8 \
  --local-jobs 12 \
  --remove-intermediate-output
```

These commands run the `biobakery wmgx` assembly and then this pipeline from the assembled contigs.  The `biobakery_workflows wmgx` command with `--run-assembly` fails in the Prokka step (unrelated to this workflow), but enough of the assembly happens beforehand that the assembly workflow can proceed afterwards.
```
hutlab load centos7/python3/biobakery_workflows/3.0.0-beta-devel-dependsUpdate
biobakery_workflows wmgx \
  --input /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/mags_and_sgbs_pipeline_testing/test_inputs/contigs_int_kneaddata/ \
  --output /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/mags_and_sgbs_pipeline_testing/test_outputs/contigs_int/ \
  --bypass-quality-control \
  --threads 8 \
  --bypass-functional-profiling \
  --bypass-strain-profiling \
  --bypass-taxonomic-profiling \
  --run-assembly \
  --grid-jobs 8 \
  --grid-scratch /n/holyscratch01/nguyen_lab/wnickols/mags_and_sgbs_pipeline_testing/contigs_int/ \
  --grid-partition shared \
  --input-extension fastq \
  --grid-options="--account=nguyen_lab"
  
hutlab unload
conda activate biobakery_assembly
export CHECKM_DATA_PATH=$(pwd)/databases/checkm/
export PHYLOPHLAN_PATH=$(pwd)/databases/phylophlan/

python assembly_workflow.py \
  -i /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/mags_and_sgbs_pipeline_testing/test_inputs/contigs_int_kneaddata/ \
  -o /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/mags_and_sgbs_pipeline_testing/test_outputs/contigs_int/ \
  --abundance-type by_sample --input-extension fastq --paired concatenated \
  --grid-scratch /n/holyscratch01/nguyen_lab/wnickols/mags_and_sgbs_pipeline_testing/contigs_int/ \
  --grid-partition 'shared' --grid-jobs 96 --cores 8 --time 10000 --mem 20000 \
  --local-jobs 12 \
  --skip-contigs \
  --remove-intermediate-output
```
