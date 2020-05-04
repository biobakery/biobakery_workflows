#!/usr/bin/env python

# Using the gene families or pathways files for a set of rna/dna samples,
# compute the normalization. This is the starr method.

import sys
import os
import argparse

def parse_arguments(args):
    """ 
    Parse the arguments from the user
    """
    parser = argparse.ArgumentParser(
        description= "RNA/DNA normalization\n",
        formatter_class=argparse.RawTextHelpFormatter,
        prog="rna_dna_norm")
    parser.add_argument(
        "--input-rna", 
        help="a tab-delimited abundance table of RNA data", 
        metavar="<input_rna.tsv>", 
        required=True)
    parser.add_argument(
        "--input-dna", 
        help="a tab-delimited abundance table of DNA data", 
        metavar="<input_dna.tsv>", 
        required=True)
    parser.add_argument(
        "--output", 
        help="a folder to write the output files", 
        metavar="<output_directory>", 
        required=True)
    parser.add_argument(
        "--mapping", 
        help="a tab-delimited mapping file of RNA to DNA samples (required if RNA/DNA sample names differ)", 
        metavar="<mapping.tsv>")
    parser.add_argument(
        "--min-abundance", 
        help="the minimum abundance value, data below this value are considered zero (default: 1)", 
        default=1)
    parser.add_argument(
        "--reduce-sample-name",
        help="remove the extra strings from the sample names",
        action='store_true')
    
    return parser.parse_args()
    
def read_table(file, min_abundance, reduce_sample_name, delimiter="\t"):
    """ Read the table from a text file with the first line the column
    names and the first column the row names. Filter humann un-identifiers and split
    data into stratified categories. """

    def format_data(x):
        """ Format the data as float, filtering based on min abundance """
        try:
            x = float(x)
        except ValueError:
            x = 0.0
            
        return x if x >= min_abundance else 0.0
    
    if reduce_sample_name:
        format_sample_name = lambda x: x.split("_Abundance")[0]
    else:
        format_sample_name = lambda x: x
    
    stratified_features={}
    unstratified_features={}
    strat_no_unclass_features={}
    
    stratified_data=[]
    unstratified_data=[]
    strat_no_unclass_data=[]
    
    stratified_index=0
    unstratified_index=0
    strat_no_unclass_index=0    
    
    with open(file) as file_handle:
        sample_names = {format_sample_name(name):i for i, name in enumerate(file_handle.readline().rstrip().split(delimiter)[1:])}
        for line in file_handle:
            line=line.rstrip().split(delimiter)
            if not "UNINTEGRATED" in line[0] and not "UNMAPPED" in line[0] and not "UNGROUPED" in line[0]:
                feature_name=line.pop(0)
                formatted_data=[format_data(i) for i in line]
                if "|" in feature_name:
                    # this is a stratified feature with species information
                    if not "|unclassified" in feature_name.lower():
                        strat_no_unclass_features[feature_name]=strat_no_unclass_index
                        strat_no_unclass_index+=1
                        strat_no_unclass_data.append(formatted_data)
                    stratified_features[feature_name]=stratified_index
                    stratified_index+=1
                    stratified_data.append(formatted_data)
                else:
                    unstratified_features[feature_name]=unstratified_index
                    unstratified_index+=1
                    unstratified_data.append(formatted_data)
                                    
    return ( unstratified_features, unstratified_data, stratified_features, 
             stratified_data, strat_no_unclass_features, strat_no_unclass_data,
             sample_names )

def read_mapping(file, delimiter="\t"):
    """ Read the mapping file """
        
    data=[]
    with open(file) as file_handle:
        for line in file_handle:
            data.append([item.strip() for item in line.rstrip().split(delimiter)])
                
    return data

def get_sample_feature(sample,feature,sample_labels,feature_labels,data):
    """ Return the data set for the sample, return zeros if name not found """
    
    # get sample index
    sample_index = sample_labels.get(sample,"NA")
    feature_index = feature_labels.get(feature,"NA")
    
    try:
        value=data[feature_index][sample_index]
    except (TypeError, IndexError):
        value=0
        
    return value
                
def divide_by_sample_total_abundance(data):
    """ For each sample in the data, divide all values by the total abundance 
    for the sample. """

    for j in range(len(data[0])):
        column_sum=sum(row[j] for row in data)
        if not column_sum == 0:
            for i in range(len(data)):
                data[i][j] = data[i][j] / column_sum
            
def compute_rna_dna_norm(samples, rna_features, rna_samples, rna_data,
    dna_features, dna_samples, dna_data):
    """ Divide the rna value by the corresponding dna value. First
    normalize by sample abundance. """

    # divide each value by the total abundance for the sample
    print("Normalize DNA")
    divide_by_sample_total_abundance(dna_data)
    print("Normalize RNA")
    divide_by_sample_total_abundance(rna_data)

    # compute norm
    norm_features=list(set(rna_features).union(dna_features))
    # create a list of empty lists, can't use [[]]*N because the empty list
    # are all pointers to the same list
    norm_data=[]
    for i in range(len(norm_features)):
        norm_data.append([])
        
    # compute the norm for all features for each sample
    for sample_name in samples:
        for index, feature in enumerate(norm_features):
            rna_value=get_sample_feature(sample_name,feature,rna_samples,rna_features,rna_data)
            dna_value=get_sample_feature(sample_name,feature,dna_samples,dna_features,dna_data)
            try:
                norm=rna_value/dna_value
            except ZeroDivisionError:
                if rna_value == 0:
                    norm="NaN"
                else:
                    norm="Inf"
            norm_data[index].append(norm)
        
    return norm_data, norm_features

def write_file(column_labels, row_labels, data, file):
    """ Write the data to a tab delimited file """
    
    # first order the rows alphabetically
    row_order=sorted(range(len(row_labels)),key=lambda i: row_labels[i])

    with open(file, "wb") as file_handle:
        file_handle.write("\t".join(column_labels)+"\n")
        for i in row_order:
            file_handle.write("\t".join([row_labels[i]]+[str(j) for j in data[i]])+"\n")
    
def main():
    # Parse arguments from command line
    args=parse_arguments(sys.argv)
    
    # read in the rna table, organizing by stratification
    print("Reading RNA table")
    ( rna_unstrat_features, rna_unstrat_data, rna_strat_features, rna_strat_data, 
      rna_strat_no_unclass_features, rna_strat_no_unclass_data, rna_samples ) = read_table(args.input_rna, args.min_abundance, args.reduce_sample_name)
    
    # read in the dna table, organizing by stratification
    print("Reading DNA table")
    ( dna_unstrat_features, dna_unstrat_data, dna_strat_features, dna_strat_data, 
      dna_strat_no_unclass_features, dna_strat_no_unclass_data, dna_samples ) = read_table(args.input_dna, args.min_abundance, args.reduce_sample_name)

    # apply sample map if provided
    if args.mapping:
        map_data = read_mapping(args.mapping)
        # convert the dna/rna sample names to the map name
        # map should be rna sample then dna sample names
        for data in map_data:
            # ignore comments or headers
            if not data[0].strip().startswith("#"):
                # find the name in the dna columns and replace
                try:
                    rna_name, dna_name = data
                except ValueError:
                    continue
                try:
                    dna_samples[rna_name]=dna_samples[dna_name]
                    del dna_samples[dna_name]
                except (ValueError, KeyError):
                    print("Mapping name "+dna_name+" is not included in the dna input file.")
    
    # get the intersection of the samples
    # apply list to set ordering
    samples=list(set(rna_samples).intersection(dna_samples))
    
    # check for matching sample names
    if len(samples) == 0:
        message="The rna/dna sample names do not match. "
        if args.mapping is None:
            message+=" Please provide a mapping file."  
        else:
            message+="Please check the formatting of the mapping file."
        sys.exit(message)
            
    # compute the norm for three sets of values: unstratified, all stratified, and
    # the stratified minus the unclassified set
    print("Compute unstratified features")
    norm_unstrat_data, norm_unstrat_features=compute_rna_dna_norm(samples, rna_unstrat_features, rna_samples,
        rna_unstrat_data, dna_unstrat_features, dna_samples, dna_unstrat_data)
    print("Compute stratified features")
    norm_strat_data, norm_strat_features=compute_rna_dna_norm(samples, rna_strat_features, rna_samples,
        rna_strat_data, dna_strat_features, dna_samples, dna_strat_data)
    print("Compute only classified features")
    norm_strat_no_unclass_data, norm_strat_no_unclass_features=compute_rna_dna_norm(samples,  
        rna_strat_no_unclass_features, rna_samples, rna_strat_no_unclass_data, 
        dna_strat_no_unclass_features, dna_samples, dna_strat_no_unclass_data)
    
    # create the output folder if needed
    if not os.path.isdir(args.output):
        print("Creating output directory: " + args.output)
        try:
            os.mkdir(args.output)
        except EnvironmentError:
            sys.exit("Unable to create output directory.")
    
    # write the output files
    output_unstrat_file=os.path.join(args.output,"rna_dna_relative_expression_unstratified.tsv")
    output_strat_file=os.path.join(args.output,"rna_dna_relative_expression.tsv")
    output_strat_no_unclass_file=os.path.join(args.output,"rna_dna_relative_expression_no_unclassifed.tsv")
    
    print("Writing unstratified table")
    write_file(["# features"]+samples, norm_unstrat_features, norm_unstrat_data,
        output_unstrat_file)
    print("Output file written: "+output_unstrat_file)
    
    print("Writing stratified table")
    write_file(["# features"]+samples, norm_unstrat_features+norm_strat_features,
        norm_unstrat_data+norm_strat_data, output_strat_file)
    print("Output file written: "+output_strat_file)
    
    print("Writing only classified table")
    write_file(["# features"]+samples, norm_unstrat_features+norm_strat_no_unclass_features,
        norm_unstrat_data+norm_strat_no_unclass_data, output_strat_no_unclass_file)
    print("Output file written: "+output_strat_no_unclass_file)
        
if __name__ == "__main__":
    main()
         
        
