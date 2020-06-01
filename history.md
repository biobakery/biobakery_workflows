
# bioBakery workflows History #

## v3.0.0-alpha.4 06-01-2020 ##

* Add utility visualization scripts, two from breadcrumbs, and new exe to access Rscripts directly.

## v3.0.0-alpha.3 05-28-2020 ##

* One change in utility function (map to list) for python 3 compatibility with visualization workflows.

## v3.0.0-alpha.2 05-19-2020 ##

* Fix join taxonomic tables script to allow for scientific notation in MetaPhlAn v3 outputs.

## v3.0.0-alpha.1 05-06-2020 ##

* Modifications to work with StrainPhlAn v3.

## v3.0.0-alpha 05-05-2020 ##

* Updated to run with HUMAnN v3 and MetaPhlAn v3 (including StrainPhlAn).
* Modified to run with python 3.

## v0.15.0 04-24-2020 ##

* Added WDL workflow.
* Increased default strainphlan slurm memory request to allow for larger runs.
* Change formatting of panphlan task command to allow for grid runs.

## v0.14.3 03-31-2020 ##

* Allow for single QC database (fix ratio table in wmgx vis).
* Add target for strainphlan best tree if created with generic clade name.

## v0.14.2 03-23-2020 ##

* Fix error when running wmgx_vis with metadata for ec heatmaps.

## v0.14.1 03-20-2020 ##

* Added heatmaps for ecs to wmgx vis.
* Changed functional profiling outputs to optional in wmgx vis.
* In demultiplex script make the index file optional.
* Add script to pull out reads mapping to metaphlan2 species by marker.
* Add two utility scripts for renaming tables and fastq files to sample ids.
* Add script to update anadama2 database with new files.
* Add workflow burst script. 

## v0.14.0 01-08-2020 ##

* Added isolate workflow.
* Increase kneaddata tasks time/memory to allow for new kneaddata feature which reorders sequences.
* Add option to use bz2 input files for wmgx workflows.
* Update dada2 ASV taxonomy table to sync ASV ids.

## v0.13.5 11-08-2019 ##

* Remove the taxonomy from the picrust2 input file to resolve the int/str error.

## v0.13.4 11-07-2019 ##

* Add PICRUSt v2 option (now also included as a task for the dada2 method).
* Update 16s workflow for python3 compatibility.
* Modify extract_markers.py strainphlan task to use new folder naming convention (to work with latest metaphlan2/strainphlan packages).

## v0.13.3 10-16-2019 ##

* Change kneaddata tasks to work with latest version as to not overwrite final pairs output file.
* Data2 workflow new options added: minoverlap and maxmismatch.
* Require min length for cutadapt to prevent reads of zero length passed to dada2 tasks which will cause an error.
* New option to allow wmgx workflow to just run panphlan.

## v0.13.2 08-12-2019 ##

* Default kneaddata (QC tasks) no longer use the rRNA database for filtering
* Default kneaddata runs with trf to filter repeats
* For database install, the metaphlan2 folder name has been updated to work with the latest metaphlan2 version
* Users can now run without any filtering databases for kneaddata tasks
* Users can now pass custom arguments to humann2 tasks

## v0.13.1 06-17-2019 ##

* Add option to provide a list of strains to run (instead of default of top 25 by average abundance).

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
