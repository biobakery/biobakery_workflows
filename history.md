
# bioBakery workflows History #

## v0.2.0 07-18-2017 ##

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
