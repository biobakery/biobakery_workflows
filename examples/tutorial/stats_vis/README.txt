
# Run the following commands to generate the output visualizations.

# This data set was generated from the bioBakery v3 output files provided by the IBDMDB site: https://ibdmdb.org/ .

# The metadata was scaled down to only a couple variables. Samples were removed from the data set if they were missing fields for any of the selected metadata variables. The taxonomic abundance was reduced to just those species and genera included in the top 15 most abundant based on average abundance plus any taxa included in the top 50 associations from running MaAsLin2. The pathways abundance was reduced to just include the top 50 pathways based on significant associations from running MaAsLin2. The Kneaddata read counts file was generated from the kneaddata logs provided on the site. The log file was generated from a demo of running five samples through the standard shotgun workflow.

$ biobakery_workflows stats --input input/ --input-metadata input/HMP2_metadata_subset.tsv --fixed-effects="diagnosis,antibiotics,age" --random-effects="site,subject" --project-name HMP2 --output HMP2_stats_output --static-covariates="age" --permutations 10 --maaslin-options="reference='diagnosis,nonIBD'" --author-name LM

$ biobakery_workflows vis --input input/ --input-metadata input/HMP2_metadata_subset.tsv --project-name HMP2 --output HMP2_vis_output --author-name LM

# The data is further reduced to just include samples for subjects at the first timepoint. 

$ biobakery_workflows stats --input input_subset/ --input-metadata input_subset/HMP2_metadata_subset_multivariable.tsv --fixed-effects="diagnosis,antibiotics,age" --project-name HMP2 --output HMP2_stats_subset_output --static-covariates="age" --maaslin-options="reference='diagnosis,nonIBD'" --author-name LM
