
# Run the following commands to generate the output visualizations

$ biobakery_workflows stats --input input/ --input-metadata input/HMP2_metadata.tsv --fixed-effects="diagnosis,dysbiosisnonIBD,dysbiosisUC,dysbiosisCD,antibiotics,age" --random-effects="site,subject" --project-name HMP2 --output HMP2_stats_example_output --static-covariates="age" --permutations 10 --maaslin-options="reference='diagnosis,nonIBD'" --author-name LM

$ biobakery_workflows stats --input input/ --input-metadata input/HMP2_metadata_subset.tsv --fixed-effects="diagnosis,antibiotics,age" --random-effects="site,subject" --project-name HMP2 --output HMP2_stats_example_output --static-covariates="age" --permutations 10 --maaslin-options="reference='diagnosis,nonIBD'" --author-name LM

$ biobakery_workflows vis --input input/ --input-metadata input/HMP2_metadata_subset.tsv --project-name HMP2 --output HMP2_vis_output --author-name LM
