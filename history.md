
# bioBakery workflows History #

## v0.13.0 05-24-2019 ##

* Add panphlan optional tasks to wmgx workflow.
* Add ITS option to 16s workflow.

## v0.12.1 01-11-2019 ##

* Add colorbar for continuous data to PCoA plots.

## v0.12.0 10-22-2018 ##

* In 16s usearch/vsearch methods combine truncate and filter (needed for data sets that require more filtering).
* Add option to bypass msa generation in 16s workflows. 

## v0.11.0 10-18-2018 ##

* Add vsearch as 16s workflow option (now the default).
* Add assembly as option to wmgx workflow.
* Add fasttree task to 16s workflow.
* Add dual indexing option to 16s and wmgx workflows.

## v0.10.0 08-21-2018 ##

* Add DADA2 16s workflow option.
* Add option to bypass taxonomic profiling for the wmgx workflow.
* Strainphlan option now selects top species by average abundance for profiling.
* Add genera visualizations to both 16s and shotgun workflows.
* Add average grouped metadata plots for relative abundance.

## v0.9.0 10-19-2017 ##

* Add metadata input option to wmgx and 16s visualization workflows.
* Add multiple sequence alignment task for closed reference sequences to 16s workflow.
* Include 16s data products in 16s report archive.
* Update 16s visualization workflow to pull variables from data processing workflow log to write report introduction.
* Add optional picard input files to 16s visualization workflow.
* Add options to bypass quality control and functional profiling to wmgx workflow.
* Improve error message in database install script printed when required dependencies are not found.
* Add counter to track PyPI download stats.
* Update dependencies to require anadama2 v0.4.0.

## v0.3.1 07-27-2017 ##

* Change the default size of heatmaps in reports based on format. Increase the size in pdf and decrease in html.

## v0.3.0 07-25-2017 ##

* Added utility to automatically install databases for each data processing workflow.
* Reorganized shotgun output file locations so the products are stored with respect to the software used to generate them.

## v0.2.0 07-18-2017 ##

* Shotgun task names updated to include sample names.
* Strain profiling was added to the shotgun data processing workflows.
* Executables are tracked for tasks in the data processing workflows.
* Tutorial data files were added to the examples folder.
* Eestats table added to 16s report.
* Discordant alignments are now allowed in qc shotgun tasks (kneaddata v0.6.1+ now required).
* Removed intermediate output option added to shotgun data processing workflows to reduce output size.
* Table of contents added to reports.
* Added workflow information (ie commands, software versions) to reports.
* Reports updated to allow for large numbers of samples.
* Updated qc tables to use serial filtering to show serial filtering in report.
* Added rna/dna norm to data processing and visualization workflows.
* Refactored visualization input to use output of data processing to allow for more input files.
* Rna/dna renorm script refactored to increase speed (minutes vs days).
* Time/memory equations added for compute intense tasks.
* Reorganizing and polishing of reports.
* Allow for different template formats and report formats.
* Allow for custom comtaminate databases in the shotgun reports.
* Added archive generation to visualization workflows.
* Allow for compressed 16s paired end files as input.
* Allow for fasta shotgun files as input.
* Initial 16s workflow added.
* Initial visualization workflows added.

## v0.0.1 12-15-2016 ##

* Initial shotgun data processing workflows added.
