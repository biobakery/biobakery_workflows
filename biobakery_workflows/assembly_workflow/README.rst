**Assembly workflow**
========================
This assembly workflow produces quality-controlled metagenomic assembled genomes (MAGs) and profiles the composition of microbial communities from metagenomic shotgun sequencing. The pipeline is based on Aaron Walsh's pets assembly workflow and `Segata et al. 2019 <https://doi.org/10.1016/j.cell.2019.01.001>`_. For queries, contact Will Nickols (willnickols@college.harvard.edu).

**Overview**
................
Often, a metagenomic sample contains new microbes not present in reference taxonomy databases. To characterize these new microbes, it is often helpful to determine their genomes in order to assess their functional abilities and taxonomy relative to known microbes. This workflow provides a method of building such genomes from cleaned metagenomic sequencing files by doing the following:

#. Assemble cleaned reads into contigs with MEGAHIT.
#. Align the reads to these contigs with Bowtie2 to determine contig coverage.
#. Bin the contigs into MAGs with MetaBAT 2.
#. Calculate the per-sample MAG abundance from the contig coverage.
#. Quality control the MAGs with CheckM2.
#. Taxonomically assign MAGs if possible with PhyloPhlAn 3.
#. Recluster unassigned high- and medium-quality MAGs into SGBs with Mash.
#. Create an abundance table from the taxonomically assigned MAGs and unassigned SGBs.

-------

.. contents:: **Table of Contents**

-------

**Installation**
................

The workflow can be installed with the following commands.  It is okay if the Checkm2 database install throws an error starting with ``File "checkm2", line 244, in <module>...``.  The PhyloPhlAn run will intentionally fail in order to download the database, but it takes a while silently preparing the database once downloaded.

::
  
    git clone https://github.com/WillNickols/assembly_workflow
    cd assembly_workflow
    conda env create -f assembly_environment.yml
    conda activate biobakery_assembly
    checkm2 database --download --path databases/checkm/
    export CHECKM_DATA_PATH=$(pwd)/databases/checkm/
    mkdir tmp_to_delete && mkdir -p databases/phylophlan && touch tmp_to_delete/tmp.fa && phylophlan_metagenomic -d SGB.Jul20 --database_folder databases/phylophlan/ -i tmp_to_delete/; rm -r tmp_to_delete/
    export PHYLOPHLAN_PATH=$(pwd)/databases/phylophlan/

Once the conda environment is created and you are in the ``assembly_workflow`` directory, you can activate the environment with these commands:

::

    conda activate biobakery_assembly
    export CHECKM_DATA_PATH=$(pwd)/databases/checkm/
    export PHYLOPHLAN_PATH=$(pwd)/databases/phylophlan/

**Tutorial**
................
  
**Example download**
--------------------
The example files can be downloaded and unzipped using the following commands:

::

    wget http://huttenhower.sph.harvard.edu/biobakery_demo/biobakery_workflows/assembly/example.tar.gz
    tar -xf example.tar.gz

**Example 1**
-------------
The following command runs the workflow on a set of paired-end fastq files that have already been cleaned with Kneaddata.  The ``by_sample`` flag specifies that abundance estimates of the sample's taxa should come from aligning that sample's reads only to its contigs.  Alternatively, the ``by_dataset`` flag can be used to align a sample's reads against all contigs generated from the input folder. Using the ``by_dataset`` flag can be useful when the samples likely share the same taxa but not every sample will assemble the shared genomes.

::

    python assembly_workflow.py \
      -i example/input/sample_1/ \
      -o example/output/sample_1/ \
      --abundance-type by_sample \
      --input-extension fastq.gz \
      --paired paired \
      --pair-identifier _R1 \
      --cores 4 \
      --local-jobs 8 \
      --remove-intermediate-output

**Abundance output**
^^^^^^^^^^^^^^^^^^^^
The output, ``example/output/sample_1/final_profile_by_sample.tsv`` is a MetaPhlAn-like output with the abundance of each taxonomic group identified in the sample.  We can inspect it with:

::

    head example/output/sample_1/final_profile_by_sample.tsv

**Phylogenetic output**
^^^^^^^^^^^^^^^^^^^^
The output shows that *Pseudoalteromonas marina* is present in the sample along with the assembled genome ``sgb_01`` representing a new species genome bin (SGB). We can check the PhyloPhlAn placement of this new SGB by examining ``example/output/sample_1/main/phylophlan/phylophlan_relab.tsv``:

::

    head example/output/sample_1/main/phylophlan/phylophlan_relab.tsv

This confirms that the closest SGB, GGB, and FGB have Mash distances of more than 0.05, 0.15, and 0.3 respectively. 

**Quality output**
^^^^^^^^^^^^^^^^^^^^
We can also check the Checkm2 report to determine the quality of these bins:

::

    head example/output/sample_1/main/checkm/qa/checkm_qa_and_n50.tsv

This file will contain the quality metrics for all assembled bins regardless of their quality. However, the bins used for the final abundance table are only the medium- and high-quality bins.

**Genome bins**
^^^^^^^^^^^^^^^^^^^^
As seen in the output file tree below, the bins are in ``example/output/sample_1/main/bins/sample_1/bins/``. 

**Example 2**
-------------
We might want to create genome bins after running a standard biobakery workflow. In this case, we can run the SGB workflow on pre-created contigs such as from the ``biobakery_workflows wmgx`` workflow with the ``--run-assembly`` flag. Here, we'll start from the contigs in ``example/output/sample_2/assembly/main/sample_2/sample_2.contigs.fa``. Note that the original read files are still required since we need to perform alignment for the abundance calculation.

::
  
    python assembly_workflow.py \
      -i example/input/sample_2/ \
      -o example/output/sample_2/ \
      --abundance-type by_sample \
      --input-extension fastq.gz \
      --paired concatenated \
      --skip-contigs \
      --cores 4 \
      --local-jobs 8 \
      --remove-intermediate-output

**Example 3**
-------------
In the ``tutorial`` folder of this GitHub, the ``tutorial/animal_guts_profile.tsv`` file is an example output from a set of 62 diverse animal stool samples.

**New SGBs**
^^^^^^^^^^^^
  
To find the number of new SGBs, We can check the number of times 'sgb' appears in the first column:

::

    awk -F'\t' '$1 ~ /sgb/ {count++} END {print count}' tutorial/animal_guts_profile.tsv

We see there were 93 new SGBs.

We can also see that some new SGBs show up in multiple samples:

::
                            
    awk -F'\t' '$1 == "\"sgb_89\"" {print}' tutorial/animal_guts_profile.tsv

Samples 2 and 4 had MAGs that were close enough that they were merged into the same novel SGB. In fact, both of these samples came from the same fin whale.

**Proportion unknown**
^^^^^^^^^^^^^^^^^^^^^^
Finally, we can visualize how much of each sample's abundance is made of known microbes, new SGBs, and unknown microbes. The following script will produce a ``figures`` folder in the ``tutorial`` folder, from which you can examine the unknown abundance.

::

    cd tutorial/
    Rscript abundance_script.R

We can see that the vast majority of most samples consists of unknown genetic material. Partially, this is due to the fact that wild animal guts are not very well characterized, but it is also due to the fact that assembly methods tend to have low recall. 

**Output file tree**
................

The folder specified by ``-o`` will have the following important files:

::
                            
    - anadama.log (log of commands and outputs)
    - final_profile_by_[sample/dataset].tsv (MetaPhlAn-like abundance table)
    - main/
      - abundance_by_[sample/dataset]/
        - [sample_name].abundance.tsv (abundance estimates of MAGs in this sample)
        - [sample_name].coverage.tsv (per-congig coverage in this sample)
        - [sample_name].mapped_read_num.txt (number of reads mapping to contigs in this sample)
        - [sample_name].total_read_num.txt (total reads in this sample)
      - assembly/
        - main/
          - [sample_name]/
            - [sample_name].final.contigs.fa (fasta file of contigs for this sample)
      - bins/
        - [sample_name]/
          - bins/
            - [sample_name].bin.[bin number].fa (one MAG from this sample)
      - checkm/
        - qa/
          - checkm_qa_and_n50.tsv (Checkm2 quality information for each MAG)
      - phylophlan/
        - phylophlan_relab.tsv (PhyloPhlAn taxonomic information for each MAG)
      - sgbs/ (for MAGs not assigned by PhyloPhlAn)
        - sgbs/
          - SGB_info.tsv (Information on which bins are in which SGBs and which bin represents the SGB)
          - sgb_[SGB number].fa (SGB representative genome)

