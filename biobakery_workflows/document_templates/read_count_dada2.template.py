#' # Read Count

#+ echo=False
total_reads=[row[0] for row in data]
sorted_samples, sorted_total_reads = utilities.sort_data(total_reads, samples)
sorted_all_read_data = [data[samples.index(sample)] for sample in sorted_samples]

filtered_reads = [row[1] for row in sorted_all_read_data]
merged_reads = [row[3] for row in sorted_all_read_data]
tabled_reads = [row[4] for row in sorted_all_read_data]
nochim_reads = [row[5] for row in sorted_all_read_data]

# plot the read counts
document.plot_grouped_barchart([sorted_total_reads,filtered_reads,merged_reads,tabled_reads,nochim_reads],
       ["Original","Filtered","Merged","Tabled","Nochimera"], sorted_samples,
       title="Read counts by Sample", xlabel="Samples", ylabel="Total Reads")

# plot grouped taxonomy for all categorical data provided
if visualizations.metadata_provided(vars):
    categorical_metadata, ordered_sorted_data, ordered_metadata, samples_found = visualizations.merge_categorical_metadata(vars, sorted_samples,
        [total_reads,filtered_reads,merged_reads,tabled_reads,nochim_reads])
    for cat_metadata in categorical_metadata:
        visualizations.plot_grouped_taxonomy_subsets(document, ordered_sorted_data, cat_metadata,
            ["total","filtered","merged","tabled","nochimera"], samples_found,
            title="Read counts by sample", ylabel="Total Reads", legend_title="")

#' <% print(visualizations.Sixteen_S.captions["dada2countsinfo"]) %>

#' <% if pdf_format: print("\clearpage") %>
