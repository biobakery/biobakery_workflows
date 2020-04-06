# MTX workflow

A mtx workflow based on the bioBakery wmgx workflows.

For more information about the bioBakery wmgx workflows, including a detailed diagram of tasks,  please see the user manual sections (wmgx and wmgx_vis): https://bitbucket.org/biobakery/biobakery_workflows/wiki/Home

### Inputs

The workflow has five required inputs and seven optional inputs. 

#### Required inputs
The workflow requires ten inputs for each run. Five inputs can be modified for each project where as five would only be modified with software version changes.
* ProjectName : The name of the sequencing project. The final output report and zip archive will use this name (only alphanumeric characters allowed).
* InputExtension : The extension for all of the input files (example ".fastq.gz")
* InputRead1Identifier : The identifier in the file name for those files that are read1 (example ".R1")
* InputRead2Identifier : The identifier in the file name for those files that are read2 (example ".R2")
* InputRead1Files : A file path (in google bucket) with a list of all of the read1 files. This file must have the full paths to all of the files and is only expected to include the read1 files (not those for read2). The names for each of the samples will be computed based on the read pair identifier and the input file extension provided. For example a file named SAMPLE1.R1.fastq.gz would have a sample name of "SAMPLE1", a read1 identifier of ".R1". and an extension of ".fastq.gz". It is expected that each sample with have two files (one file for each read of the pair). 

To generate a file to use as input for InputRead1Files, follow the Terra instructions https://support.terra.bio/hc/en-us/articles/360033353952-Creating-a-list-file-of-reads-for-input-to-a-workflow , adding to command #2 the InputRead1Identifier and the InputExtension. For example with InputRead1Identifier = ".R1" and InputExtension = ".fastq.gz" command #2 would now be 
`gsutil ls gs:/your_data_Google_bucket_id/ | grep ".fastq.gz" | grep ".R1" > ubams.list` . Also since for this workflow we are looking for fastq or fastq.gz input files you might change the name of the file list in this command from "ubams.list" to "fastq_list.txt" . 

These five required inputs would only be modified if the versions of Kneaddata and HUMAnN v2 change. These are databases that are specifically tied to the software version.
* versionSpecificChocophlan : The Chocophlan database used by HUMAnN. This is located at `databases/humann/full_chocophlan_plus_viral.v0.1.1.tar.gz` in this workspace google bucket.
* versionSpecifichumanDB : The human reference database used by Kneaddata. This is located at `databases/kneaddata/Homo_sapiens_hg37_human_contamination_Bowtie2_v0.1.tar.gz` in this workspace google bucket.
* versionSpecificrrnaDB : The human rrna reference database used by Kneaddata. This is located at `databases/kneaddata/Homo_sapiens_hg38_transcriptome_Bowtie2_v0.1.tar.gz` in this workspace google bucket. 
* versionSpecificUniRef90 : The uniref90 reference database used by HUMAnN. This is located at `databases/humann/uniref90_annotated_1_1.tar.gz` in this workspace google bucket.
* versionSpecificUtilityMapping : The utility mapping database used by HUMAnN. This is located at `databases/humann/full_utility_mapping_1_1.tar.gz` in this workspace google bucket.

#### Optional inputs
There are an additional seven optional inputs for each workflow run. These are not required. If not set, the default values will be used.
* bypassFunctionalProfiling (Default = false): This set to true will bypass running functional profiling and all of the downstream tasks including normalization and merging of the functional data products including gene families, ecs, and pathways.
* dataType (Default = "mtx"): This is the sequencing data type (mtx or mgx).
* inputMetadataFile (Default = None) : This file is used with the visualization task to annotate the figures with metadata.
* MaxMemGB_FunctionalProfileTasks (Default = 24 Gb): This is the max memory to request for each HUMAnN v2.0 task. This might need to be increased for larger input files.
* MaxMemGB_QualityControlTasks (Default = 8 Gb): This is the max memory to request for each Kneaddata task. This might need to be increased for larger input files.
* MaxMemGB_TaxonomicProfileTasks (Default = 24 Gb): This is the max memory to request for each MetaPhlAn v2.0 task. This might need to be increased for larger input files.
* preemptibleAttemptsOverride (Default = 2): This setting determines how many times to rerun one of the main compute tasks (HUMAnN2 v2.0, Kneaddata, and MetaPhlAn v2.0) on pre-emptible instances. If set to zero a non-pre-emptible instance will be used.

### Outputs

The workflow has several intermediate outputs and a final zip archive that includes a report of exploratory figures plus compiled data tables. Each task has its own folder in the google bucket output folder with a sub-folder for each time it is run. The outputs of interest, including their respective folders, are described below. `$SAMPLE_NAME` is the name of the sample included in the original raw files. For example, SAMPLE1.R1.fastq.gz would have a sample name of "SAMPLE1".

* call-QualityControl (with sub-folders for each sample, in each subfolder there are the following files)
  * `$SAMPLE_NAME.fastq.gz` : This is the file of reads after running through QC.
  * `$SAMPLE_NAME.log` : This is the log from Kneaddata that includes read counts.
  * `glob*/$SAMPLE_NAME_DB_contam*.fastq.gz` : These are the reads that mapped to the reference database (with name `$DB`) for this sample.
  * `glob*/$SAMPLE_NAME_[R1|R2].[html|zip]` : These are the output files from running fastqc on read1 and read2 prior to running quality control.
* call-FunctionalProfile (with sub-folders for each sample, in each subfolder there are the following files)
  * `$SAMPLE_NAME.log` : This is the log from the HUMAnN v2.0 that includes read alignment counts.
  * `glob*/$SAMPLE_NAME_bowtie2_unaligned.fa` : These are the unaligned reads from running the nucleotide search.
  * `glob*/$SAMPLE_NAME_diamond_unaligned.fa` : These are the unaligned reads from running the translated search.
* call-VisualizationReport
  * `$PROJECT_NAME_visualization.zip` : This folder contains a visualization report plus final compiled data tables.
    * `wmgx_report.pdf` : This is the exploratory report of tables and figures.
    * `data/humann2_feature_counts.tsv` : This contains the feature counts (pathways, gene families, ecs) for each sample.
    * `data/humann2_read_and_species_counts.tsv` : This contains the counts of reads aligning at each step plus the total number of species identified for each sample.
    * `data/kneaddata_read_count_table.tsv` : This contains the read counts (split into pairs and orphans) for each step in the quality control process for each sample.
    * `data/metaphlan2_taxonomic_profiles.tsv` : This contains the merged taxonomic profiles for all samples.
    * `data/microbial_counts_table.tsv` : This table includes counts ratios for each step of the quality control process for all samples.
    * `data/pathabundance_relab.tsv` : This is a merged table of the pathway abundances for all samples normalized to relative abundance.
    * `data/qc_counts_orphans_table.tsv` : This is table with the total number of orphan reads not aligning to each of the reference tables.
    * `data/qc_counts_pairs_table.tsv` : This is table with the total number of paired reads not aligning to each of the reference tables.
    * `data/taxa_counts_table.tsv`: This table includes the total number of species and general before and after filtering.
    * `data/top_average_pathways_names.tsv` : This table includes the top pathways by average abundance, with their full names, including average abundance and variance.

### Run a demo

A demo data set is included in this workspace. The demo set includes six paired samples (three MTX and three MGX) from IBDMDB plus a small metadata file. Using preemptive instances, this demo set will cost about $5 to run.

IBDMDB (6 sample) demo run configuration:
* ProjectName : `"ibdmdb_demo"` (this can be any string you would like)
* InputExtension : `".fastq.gz"`
* InputRead1Identifier : `"_R1"`
* InputRead2Identifier : `"_R2"`
* InputRead1Files : `"gs://fc-7130738a-5cde-4238-b00a-e07eba6047f2/IBDMDB/ibdmdb_file_list.txt"`
* inputMetadataFile : `"gs://fc-7130738a-5cde-4238-b00a-e07eba6047f2/IBDMDB/ibdmdb_demo_metadata.txt"`

Refer to the section above for descriptions of the output files generated by running the workflow.

Example output files from running the IBDMDB data set with metadata can be found in this workspace in the folder `IBDMDB/final_outputs/ibdmdb_demo_visualizations.zip`.


