# assembly_workflow

> This repository contains code for Will's assembly workflow based on Aaron Walsh's pets assembly workflow and Segata et al. 2019. For queries, contact Will Nickols (email: <willnickols@college.harvard.edu>).

# Installation

The workflow can be installed with the following commands.  It is okay if the Checkm2 database install throws an error starting with `File "checkm2/bin/checkm2", line 244, in <module>...`.  The PhyloPhlAn run will intentionally fail in order to download the database but will take a while to download the database.
```
git clone https://github.com/WillNickols/assembly_workflow
cd assembly_workflow
conda env create -f assembly_environment.yml
conda activate biobakery_assembly
checkm2 database --download --path databases/checkm/
export CHECKM_DATA_PATH=$(pwd)/databases/checkm/
mkdir tmp_to_delete && mkdir -p databases/phylophlan && touch tmp_to_delete/tmp.fa && phylophlan_metagenomic -d SGB.Jul20 --database_folder databases/phylophlan/ -i tmp_to_delete/; rm -r tmp_to_delete/
export PHYLOPHLAN_PATH=$(pwd)/databases/phylophlan/
```

Once the conda environment is created and you are in the `assembly_workflow` directory, you can activate the environment with these commands:
```
conda activate biobakery_assembly
export CHECKM_DATA_PATH=$(pwd)/databases/checkm/
export PHYLOPHLAN_PATH=$(pwd)/databases/phylophlan/
```

# Example runs

The example files can be downloaded and unzipped using the following commands:
```
wget http://huttenhower.sph.harvard.edu/biobakery_demo/biobakery_workflows/assembly/example.tar.gz
tar -xf example.tar.gz
```

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
  --remove-intermediate-output
```
