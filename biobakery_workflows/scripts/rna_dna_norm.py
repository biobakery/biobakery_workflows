#!/usr/bin/env python

# Using the gene families or pathways files for a set of rna/dna samples,
# compute the normalization. This is the starr method.

import sys
import os
import argparse
import copy

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
        "--remove-zeros", 
        help="remove features where all samples have zero abundance from the output files", 
        action="store_true")
    
    return parser.parse_args()
    
def read_table(file, labels=True, filter=False, delimiter="\t", apply_float=True):
    """ Read the table from a text file with the first line the column
    names and the first column the row names. If no labels are provided,
    then read in a single data matrix. Filter humann2 un-identifiers if set. """
        
    data=[]
    row_names=[]
    with open(file) as file_handle:
        if labels:
            column_names = file_handle.readline().rstrip().split(delimiter)[1:]
        for line in file_handle:
            line=line.rstrip().split(delimiter)
            if labels and filter:
                if not "UNINTEGRATED" in line[0] and not "UNMAPPED" in line[0] and not "UNGROUPED" in line[0]:
                    row_names.append(line.pop(0))
                    if apply_float:
                        data.append([float(i) for i in line])
                    else:
                        data.append(line)
            else:
                if labels:
                    row_names.append(line.pop(0))
                if apply_float:
                    data.append([float(i) for i in line])
                else:
                    data.append(line)
                
    if labels:
        return column_names, row_names, data
    else:
        return data
    
def get_feature(name,feature_labels,data):
    """ Return the data set for the feature, return zeros if name not found """
    
    try:
        index=feature_labels.index(name)
        values=data[index]
    except ValueError:
        values=[0]*len(data[0])
        
    return values

def get_sample(name,sample_labels,data):
    """ Return the data set for the sample, return zeros if name not found """
    
    try:
        index=sample_labels.index(name)
        values=[row[index] for row in data]
    except ValueError:
        values=[0]*len(data)
        
    return values

def filter_min_abundance(data, min_abundance):
    """ Replace values below min abundance with zero.
    Edit data matrix directly so no real need to return data (done for clarity). """
    
    for i in range(len(data)):
        for j in range(len(data[0])):
            if data[i][j] < min_abundance:
                data[i][j] = 0
                
    return data

def filter_features_from_list(data, current_features, filtered_features_list):
    """ Remove the features from the data which are not included in the
    filtered features list """
    
    data_filtered=[]
    features_filtered=[]
    for feature, row in zip(current_features, data):
        if feature in filtered_features_list:
            features_filtered.append(feature)
            data_filtered.append(row)
            
    return data_filtered, features_filtered
                
def divide_by_sample_total_abundance(data):
    """ For each sample in the data, divide all values by the total abundance 
    for the sample. """
    
    norm_data=copy.deepcopy(data)
    
    for j in range(len(data[0])):
        column_sum=sum(row[j] for row in data)
        for i in range(len(data)):
            norm_data[i][j] = data[i][j] / column_sum
            
    return norm_data
            
def compute_rna_dna_norm(samples, rna_features, rna_samples, rna_data,
    dna_features, dna_samples, dna_data):
    """ Divide the rna value by the corresponding dna value. First
    normalize by sample abundance. """

    # divide each value by the total abundance for the sample
    norm_dna_data=divide_by_sample_total_abundance(dna_data)
    norm_rna_data=divide_by_sample_total_abundance(rna_data)

    # compute norm
    norm_features=list(set(rna_features).union(dna_features))
    # create a list of empty lists, can't use [[]]*N because the empty list
    # are all pointers to the same list
    norm_data=[]
    for i in range(len(norm_features)):
        norm_data.append([])
        
    # compute the norm for all features for each sample
    for s in samples:
        # get the values for all features for this sample
        rna_sample_data=dict([(feature,value) for feature,value in zip(rna_features, get_sample(s, rna_samples, norm_rna_data))])
        dna_sample_data=dict([(feature,value) for feature,value in zip(dna_features, get_sample(s, dna_samples, norm_dna_data))])
        for index, feature in enumerate(norm_features):
            rna_value=rna_sample_data.get(feature,0)
            dna_value=dna_sample_data.get(feature,0)
            try:
                norm=rna_value/dna_value
            except ZeroDivisionError:
                norm=0
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
            
def split_stratified(feature_labels, data):
    """ Split the features by stratification, Return features/data for
    the unstratified features, the stratified feature, and the stratified features
    without the unclassified groups. """
    
    stratified_features=[]
    unstratified_features=[]
    strat_no_unclass_features=[]
    stratified_data=[]
    unstratified_data=[]
    strat_no_unclass_data=[]
    
    for feature, row in zip(feature_labels, data):
        if "|" in feature:
            # this is a stratified feature with species information
            if not "|unclassified" in feature.lower():
                strat_no_unclass_features.append(feature)
                strat_no_unclass_data.append(row)
            stratified_features.append(feature)
            stratified_data.append(row)
        else:
            unstratified_features.append(feature)
            unstratified_data.append(row)
            
    return ( unstratified_features, unstratified_data, stratified_features, 
             stratified_data, strat_no_unclass_features, strat_no_unclass_data )
    
def main():
    # Parse arguments from command line
    args=parse_arguments(sys.argv)
    
    # read in the rna table, removing humann2 "un-features"
    rna_samples, rna_features, rna_data = read_table(args.input_rna, filter=True)
    
    # read in the dna table, removing humann2 "un-features"
    dna_samples, dna_features, dna_data = read_table(args.input_dna, filter=True)

    # apply sample map if provided
    if args.mapping:
        map_data = read_table(args.mapping, labels=False, apply_float=False)
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
                    dna_samples[dna_samples.index(dna_name)]=rna_name
                except ValueError:
                    print("Mapping name "+dna_name+" is not included in the dna input file.")
    
    # get the intersection of the samples
    # apply list to set ordering
    samples=list(set(rna_samples).intersection(dna_samples))
    
    # get the union of the features
    features=list(set(rna_features).union(dna_features))
    
    # filter data to reduce those items below min abundance to zero
    dna_data=filter_min_abundance(dna_data, args.min_abundance)
    rna_data=filter_min_abundance(rna_data, args.min_abundance)
    
    # remove features from the set that are all zero in both dna and rna
    filtered_features=[]
    for f in features:
        dna_sum=sum(get_feature(f, dna_features, dna_data))
        rna_sum=sum(get_feature(f, rna_features, rna_data))
            
        if dna_sum + rna_sum > 0:
            if args.remove_zeros:
                # check if they both are not zero
                if dna_sum > 0 and rna_sum > 0:
                    filtered_features.append(f)
            else:
                filtered_features.append(f)
            
    # remove features from the rna/dna sets that are not included in the filtered list
    dna_data_filtered, dna_features_filtered = filter_features_from_list(dna_data,
        dna_features, filtered_features)
    rna_data_filtered, rna_features_filtered = filter_features_from_list(rna_data,
        rna_features, filtered_features)
    
    # split the data into three groups based on stratification
    ( dna_unstrat_features, dna_unstrat_data, dna_strat_features, dna_strat_data, 
      dna_strat_no_unclass_features, dna_strat_no_unclass_data ) = split_stratified(dna_features_filtered, dna_data_filtered)

    ( rna_unstrat_features, rna_unstrat_data, rna_strat_features, rna_strat_data, 
      rna_strat_no_unclass_features, rna_strat_no_unclass_data ) = split_stratified(rna_features_filtered, rna_data_filtered)      
            
    # compute the norm for three sets of values: unstratified, all stratified, and
    # the stratified minus the unclassified set
    norm_unstrat_data, norm_unstrat_features=compute_rna_dna_norm(samples, rna_unstrat_features, rna_samples,
        rna_unstrat_data, dna_unstrat_features, dna_samples, dna_unstrat_data)
    norm_strat_data, norm_strat_features=compute_rna_dna_norm(samples, rna_strat_features, rna_samples,
        rna_strat_data, dna_strat_features, dna_samples, dna_strat_data)
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
    
    
    write_file(["# features"]+samples, norm_unstrat_features, norm_unstrat_data,
        output_unstrat_file)
    print("Output file written: "+output_unstrat_file)
    
    write_file(["# features"]+samples, norm_unstrat_features+norm_strat_features,
        norm_unstrat_data+norm_strat_data, output_strat_file)
    print("Output file written: "+output_strat_file)
    
    write_file(["# features"]+samples, norm_unstrat_features+norm_strat_no_unclass_features,
        norm_unstrat_data+norm_strat_no_unclass_data, output_strat_no_unclass_file)
    print("Output file written: "+output_strat_no_unclass_file)
        
if __name__ == "__main__":
    main()
         
        
