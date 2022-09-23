
# bioBakery workflows History #

## v3.1 TBD ###
* Make the search for the strainphlan database folder version agnostic
* Allow for bam extension in the wmgx workflow plus add qc-scratch option (requires kneaddata v0.12.0)
* Allow for primers to be removed in the DADA2 method
* Updates for the new kneaddata with options to specify read pairs (v0.11.0)
* Update MetaPhlAn grid time to sync with newer version
* Increase kneaddata run time to allow for new gzip addition of final target
* Add a utility script to compute species of interest from raw reads and taxonomic profile
* Allow for kneaddata final target to exist and then overwrite
* Add DADA2 option to increase min folder to make chimera calling more strict
* Update to support the latest panphlan version
* Resolve issue with QC ratios GT 1 and add exceptions
* Add option to exclude a subset of reads in the subset data file script
* Allow for bypass of taxonomic profiling for the running on strain based profiling methods
* Add check to 16s workflow for max trunc len based on average read length
* Add check to make sure the fixed and random effect variables provided are include in the metadata
* For 16s, add check for many unclassified ASVs, check for empty input files, and option to allow reverse complement strands
* Allow for database installs to fail with a warning
* Visualization changes (includes changes to merge vis to a single "universal vis")
  - Update to allow for search of inputs files in input folder instead of requiring exact paths
  - Allow for quotes in the data passed to the Mantel test function
  - Allow for larger feature names
  - Decrease the min samples used for filtering
  - Allow for spaces and forward slack in chars in table headers
  - Only filter the top subset of the top features to reduce runtime
  - Fix floating point error in correlation function
  - Add filtered heatmaps to functional templates
  - Add correlation threshold option
  - Add filtered heatmaps to taxonomy template
  - Add functions to filter based on spearman correlation
  - Add heatmaps with zscores
  - Reduce margins and legends to display larger plots
  - Add metadata file to the final zip archive
  - Allow for multiple 1s in the pair identifiers
  - Update the qc counts figures to sync the names and text
  - Shorten reference database names for plots
  - Add check for spaces in metadata variables to prevent issues downstream
  - Add more detailed methods to the wmgx vis intro
  - Update the default usearch intro to include sample number
  - Allow for taxa without the full taxonomy to be used in the filter taxa function
  - Update the 16s taxa plots to sort the legend in reverse
  - Allow for taxa without the fill taxonomy to be included
  - Remove unexpected chars from tables and pathways bullet list
  - Update maaslin2 tile generation to allow for numbers
  - Make halla output heatmaps optional if not enough associations are found
  - Add option to always generate stratified outputs from picrust2 task
  - Remove pairwise comparision if there is just a single metadata variable
  - Convert metadata samples as columns in stats workflow if needed
  - If an ec file is provided without the names, add the names
  - Add max missing option
  - Add the covariate equation to the vis to match that used for the script
  - Track samples missing metadata
  - Add barcharts for pathways and ecs
  - Make author and project name optional  

## v3.0.0-alpha.7 02-19-2021 ##

* Add initial (alpha) stats workflow.
* Update stats workflow to include mantel tests and HAllA plus reorder sections. Also tile MaAsLin2 plots to fit more on one page. Add author option. Add option to add image to header. Update beta diversity table to scale with number of metadata variables. Allow for different orientations for data/metadata files for beta diversity script plus allow for special characters in sample names. In beta diversity script only consider NA values as "missing" and not zeros in metadata. Update maaslin plots from jpeg to png.
* Update panphlan tasks to work with panphlan v3.
* Update WDL to include strainphlan tasks (minus references and selection of clades by average abundance).
* Add third custom database option to WDL.
* In dada2 workflow, add option to set reverse trunc len to different value then forward and also allow for zero maxee.
* In ITS workflow, update database selection to find ITS database installed in picrust2 library folder.
* In usearch workflow, increase qmax for eestats for old illumina format and add ascii option for older format.
* Added a new utility script to remove primers. 
* In wmgx workflow, add option to provide options for taxonomic profiling. 
* In 16s workflows, change format of picrust2 inputs for asv tables to also remove taxonomic profiling information to resolve error with latest picrust2.

## v3.0.0-alpha.6 07-10-2020 ##

* Add custom usearch db option plus log to vis report for 16s.
* Small change to allow for python 3 compatibility with vis with metadata (average barplots).

## v3.0.0-alpha.5 07-06-2020 ##

* Change default version run for PICRUSt from 1 to 2.
* Update reverse comp barcodes script for python 3.
* Update gzip file read for python 3.
* Add check for feature with single metadata type to vis.
* Update strainphlan vis script for latest demo inputs (Thanks, Kelsey!).
* For PICRUSt 2 inputs remove taxonomy column to fit new requirements.

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

## v0.15.1 06-09-2020 ##

* Added check in visualization to throw error if feature only has one type.
* Increase metaphlan2 and strainphlan task times for larger grid runs.
* Added option for custom databases for 16s usearch/vsearch methods.

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
