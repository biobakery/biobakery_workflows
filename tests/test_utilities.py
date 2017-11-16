
import unittest
import tempfile
import os
import sys

from biobakery_workflows import utilities

# write to a temp file
def write_temp(data):
    """ Write the data to a temp file """
    
    handle, file = tempfile.mkstemp(prefix="biobakery_workflows_test")
    os.write(handle,data)
    os.close(handle)
    
    return file

class TestUtiltiesFunctions(unittest.TestCase):
    """ Test the functions found in the biobakery workflows utilities module """
    
    def test_taxa_remove_unclassified(self):
        """ Test the taxa remove unclassified function """
        
        taxa=["k__k1;p__p1;c__","k__k1;p__p1;c__c1;o__o1;f__f1;g__g1;s__s1",
              "k__k1;p__p1;c__c1;o__o1;f__f1;g__g1;s__"]
        expected_taxa=["k__k1;p__p1","k__k1;p__p1;c__c1;o__o1;f__f1;g__g1;s__s1",
              "k__k1;p__p1;c__c1;o__o1;f__f1;g__g1"]
        
        self.assertEqual(list(utilities.taxa_remove_unclassified(taxa)), expected_taxa)
        
    def test_taxa_remove_unclassified_spaces(self):
        """ Test the taxa remove unclassified function with spaces in names """
        
        taxa=["k__k1;p__p1 ; c__ ",
            "k__k1;p__p1;c__c1 ; o__o1 ;f__f1;g__g1 ; s__s1 ",
            "k__k1;p__p1;c__c1;o__o1;f__f1;g__g1;s__"]
        expected_taxa=["k__k1;p__p1",
            "k__k1;p__p1;c__c1;o__o1;f__f1;g__g1;s__s1",
            "k__k1;p__p1;c__c1;o__o1;f__f1;g__g1"]
        
        self.assertEqual(list(utilities.taxa_remove_unclassified(taxa)), expected_taxa)

    def test_taxa_by_level(self):
        """ Test the taxa by level function """
        
        taxa=["k__k1;p__p1;c__c1","k__k1;p__p1;c__c2","k__k1;p__p2;c__c3"]
        data=[[1,2],[1,2],[1,2]]
        expected_taxa=["k__k1;p__p1","k__k1;p__p2"]
        expected_data=[[2,4],[1,2]]
        
        actual_taxa, actual_data = utilities.taxa_by_level(taxa, data, 1)
        
        self.assertEqual(actual_taxa, expected_taxa)
        self.assertEqual(actual_data, expected_data)
        
    def test_taxa_by_level_kingdom(self):
        """ Test the taxa by level function set to merge to kingdom level"""
        
        taxa=["k__k1;p__p1;c__c1","k__k1;p__p1;c__c2","k__k1;p__p2;c__c3"]
        data=[[1,2],[1,2],[1,2]]
        expected_taxa=["k__k1"]
        expected_data=[[3,6]]
        
        actual_taxa, actual_data = utilities.taxa_by_level(taxa, data, 0)
        
        self.assertEqual(actual_taxa, expected_taxa)
        self.assertEqual(actual_data, expected_data)
        
    def test_taxa_by_level_strain(self):
        """ Test the taxa by level function set to merge to strain level"""
        
        taxa=["k__k1;p__p1;c__c1;o__o1;f__f1;g__g1;s__s1;t__t1","k__k1;p__p1;c__c2"]
        data=[[1,2],[1,2]]
        expected_taxa=["k__k1;p__p1;c__c1;o__o1;f__f1;g__g1;s__s1;t__t1"]
        expected_data=[[1,2]]
        
        actual_taxa, actual_data = utilities.taxa_by_level(taxa, data, 7)
        
        self.assertEqual(actual_taxa, expected_taxa)
        self.assertEqual(actual_data, expected_data)
        
    def test_taxa_by_level_strain_missing(self):
        """ Test the taxa by level function set to merge to strain level
            though there is no strain level input data """
        
        taxa=["k__k1;p__p1;c__c1;o__o1;f__f1;g__g1;s__s1","k__k1;p__p1;c__c2"]
        data=[[1,2],[1,2]]
        expected_taxa=[]
        expected_data=[]
        
        actual_taxa, actual_data = utilities.taxa_by_level(taxa, data, 7)
        
        self.assertEqual(actual_taxa, expected_taxa)
        self.assertEqual(actual_data, expected_data)

    def test_taxa_by_level_unknown(self):
        """ Test the taxa by level function with one taxon that does not include that level"""
        
        taxa=["k__k1;p__p1","k__k1;p__p1;c__c2","k__k1;p__p2;c__c3"]
        data=[[1,2],[1,2],[1,2]]
        expected_taxa=["k__k1;p__p1;c__c2","k__k1;p__p2;c__c3"]
        expected_data=[[1,2],[1,2]]
        
        actual_taxa, actual_data = utilities.taxa_by_level(taxa, data, 2)
        
        self.assertEqual(sorted(actual_taxa), sorted(expected_taxa))
        self.assertEqual(actual_data, expected_data)
        
    def test_relative_abundance(self):
        """ Test the relative abundance function """
        
        data=[[1,2,3],[1,2,3],[1,1,1]]
        expected_relab=[[1/3.0,2/5.0,3/7.0],[1/3.0,2/5.0,3/7.0],[1/3.0,1/5.0,1/7.0]]
        
        self.assertEqual(utilities.relative_abundance(data),expected_relab)

    def test_relative_abundance_zeros(self):
	    """ Test the relative abundance function with a sample with all zeros"""
	    
	    data=[[1,2,3,0],[1,2,3,0],[1,1,1,0]]
	    expected_relab=[[1/3.0,2/5.0,3/7.0,0],[1/3.0,2/5.0,3/7.0,0],[1/3.0,1/5.0,1/7.0,0]]
	    
	    self.assertEqual(utilities.relative_abundance(data),expected_relab)
        
    def test_taxa_shorten_name(self):
        """ Test the taxa shorten name function """
        
        taxa=["k__k1;p__p1;c__c1","k__k1;p__p1;c__c1;o__o1","k__k1;p__p1"]
        expected_taxa=["p__p1","p__p1","p__p1"]
        
        self.assertEqual(utilities.taxa_shorten_name(taxa, 1), expected_taxa)
        
    def test_taxa_shorten_name_remove_identifier(self):
        """ Test the taxa shorten name function with removing the identifier"""
        
        taxa=["k__k1;p__p1;c__c1","k__k1;p__p1;c__c1;o__o1","k__k1;p__p1"]
        expected_taxa=["p1","p1","p1"]
        
        self.assertEqual(utilities.taxa_shorten_name(taxa, 1, remove_identifier=True), expected_taxa)
        
    def test_terminal_taxa(self):
        """ Test the terminal taxa function """
        
        taxa=["k__k1;p__p1","k__k1;p__p1;c__c1","k__k2;p__p2","k__k3;p__p3;c__c2;o__o3"]
        data=[[1],[2],[3],[4]]
        
        # it is expected that order of the taxa will stay the same as the original input
        expected_taxa=["k__k1;p__p1;c__c1","k__k2;p__p2","k__k3;p__p3;c__c2;o__o3"]
        expected_data=[[2],[3],[4]]
        
        actual_taxa, actual_data = utilities.terminal_taxa(taxa, data)
        
        self.assertEqual(actual_taxa,expected_taxa)
        self.assertEqual(actual_data,expected_data)
        
    def test_terminal_taxa_duplicate(self):
        """ Test the terminal taxa function with two duplicate terminal taxa"""
        
        taxa=["k__k1;p__p1","k__k1;p__p1;c__c1","k__k2;p__p2","k__k3;p__p3;c__c2;o__o3","k__k3;p__p3;c__c2;o__o3"]
        data=[[1],[2],[3],[4],[5]]
        
        # it is expected that order of the taxa will stay the same as the original input
        expected_taxa=["k__k1;p__p1;c__c1","k__k2;p__p2","k__k3;p__p3;c__c2;o__o3"]
        expected_data=[[2],[3],[9]]
        
        actual_taxa, actual_data = utilities.terminal_taxa(taxa, data)
        
        self.assertEqual(actual_taxa,expected_taxa)
        self.assertEqual(actual_data,expected_data)

    def test_terminal_taxa_unclassified_family(self):
        """ Test the terminal taxa function with a taxon with unclassified at family level"""
        
        taxa=["k__k1;p__p1","k__k1;p__p1;c__c1","k__k2;p__p2","k__k3;p__p3;c__c2;o__o3;f__;g__;s__","k__k3;p__p3;c__c2;o__o3;f__f1"]
        data=[[1],[2],[3],[4],[5]]
        
        # it is expected that order of the taxa will stay the same as the original input
        expected_taxa=["k__k1;p__p1;c__c1","k__k2;p__p2","k__k3;p__p3;c__c2;o__o3;f__;g__;s__","k__k3;p__p3;c__c2;o__o3;f__f1"]
        expected_data=[[2],[3],[4],[5]]
        
        actual_taxa, actual_data = utilities.terminal_taxa(taxa, data)
        
        self.assertEqual(actual_taxa,expected_taxa)
        self.assertEqual(actual_data,expected_data)
        
    def test_terminal_taxa_unclassified_species(self):
        """ Test the terminal taxa function with a taxon with unclassified at species level"""
        
        taxa=["k__k3;p__p3;c__c2;o__o3;f__f1;g__g1;s__","k__k3;p__p3;c__c2;o__o3;f__f1;g__g1;s__s1"]
        data=[[1],[2]]
        
        # it is expected that order of the taxa will stay the same as the original input
        expected_taxa=["k__k3;p__p3;c__c2;o__o3;f__f1;g__g1;s__","k__k3;p__p3;c__c2;o__o3;f__f1;g__g1;s__s1"]
        expected_data=[[1],[2]]
        
        actual_taxa, actual_data = utilities.terminal_taxa(taxa, data)
        
        self.assertEqual(actual_taxa,expected_taxa)
        self.assertEqual(actual_data,expected_data)
        
    def test_filter_zero_rows(self):
        """ Test the filter zero rows function """
        
        taxa=["k__k1;p__p1","k__k1;p__p1;c__c1","k__k2;p__p2","k__k3;p__p3;c__c2;o__o3","k__k3;p__p3;c__c2;o__o3"]
        data=[[1],[0],[3],[0],[5]]
        
        # it is expected that order of the taxa will stay the same as the original input
        expected_taxa=["k__k1;p__p1","k__k2;p__p2","k__k3;p__p3;c__c2;o__o3"]
        expected_data=[[1],[3],[5]]
        
        actual_taxa, actual_data = utilities.filter_zero_rows(taxa, data)
        
        self.assertEqual(actual_taxa,expected_taxa)
        self.assertEqual(actual_data,expected_data)     
        
    def test_sort_data(self):
        """ Test the sort data function """
        
        samples=["s1","s3","s2"] 
        data=[1,3,2] 
        
        expected_samples=["s3","s2","s1"]
        expected_data=[3,2,1]
        
        actual_samples, actual_data = utilities.sort_data(data, samples)
        
        self.assertEqual(actual_data, expected_data)
        self.assertEqual(actual_samples, expected_samples)
        
    def test_sort_data_list_of_lists(self):
        """ Test the sort data function with data as a list of lists """
        
        samples=["s1","s3","s2"] 
        data=[[1],[3],[2]] 
        
        expected_samples=["s3","s2","s1"]
        expected_data=[3,2,1]
        
        actual_samples, actual_data = utilities.sort_data(data, samples)

        self.assertEqual(actual_data, expected_data)
        self.assertEqual(actual_samples, expected_samples)
        
    def test_name_files(self):
        """ Test the name files function """
        
        samples=["sample1","sample2"]
        
        self.assertEqual(sorted(utilities.name_files(samples,"/tmp")), ["/tmp/sample1","/tmp/sample2"])
        
    def test_name_files_single_file(self):
        """ Test the name files function with a single file to make sure a string instead of list is returned """
        
        self.assertEqual(utilities.name_files("sample1","/tmp"), "/tmp/sample1")
        
    def test_taxonomy_trim(self):
        """ Test the taxonomy trim function """
        
        taxa = ["k__k3;p__p3;c__c2;o__o3;f__;g__;s__","k__k3;p__p3;c__c2;o__o3;f__f1;g__g1;s__",
                "k__k3;p__p3;c__c2;o__o3;f__f1;g__g1;s__s1","k__k3;p__p3;c__c2;o__o3;f__f1;g__;s__"]
        
        expected_taxa = ["o__o3.f__.g__.s__","g__g1.s__","g__g1.s__s1","f__f1.g__.s__"]
        
        self.assertEqual(utilities.taxonomy_trim(taxa), expected_taxa)
        
    def test_sample_names(self):
        """ Test the sample names function """
        
        files=["s1.fastq","s2.fastq"]
        
        self.assertEqual(utilities.sample_names(files,".fastq"), ["s1","s2"])
        
    def test_sample_names_extension(self):
        """ Test the sample names function without leading period"""
        
        files=["s1.fastq","s2.fastq"]
        
        self.assertEqual(utilities.sample_names(files,"fastq"), ["s1","s2"])
        
    def test_sample_names_period(self):
        """ Test the sample names function with period included in name"""
        
        files=["s1.1.fastq","s2.1.fastq"]
        
        self.assertEqual(utilities.sample_names(files,"fastq"), ["s1.1","s2.1"])
        
    def test_sample_names_pair_identifier(self):
        """ Test the sample names function with a pair identifier """
        
        files=["s1.R1.fastq","s2.R1.fastq"]
        
        self.assertEqual(utilities.sample_names(files,".fastq",".R1"), ["s1","s2"])
        
    def test_sample_names_pair_identifier_not_found(self):
        """ Test the sample names function with a pair identifier that is not included in the names """
        
        files=["s1.R1.fastq","s2.R1.fastq"]
        
        self.assertEqual(utilities.sample_names(files,".fastq","_R1"), ["s1.R1","s2.R1"])
        
    def test_paired_files(self):
        """ Test the paired files function """
        
        files=["s1.R1.fastq","s1.R2.fastq","s2.R1.fastq","s2.R2.fastq"]
        
        expected_pairs = [["s1.R1.fastq","s2.R1.fastq"],["s1.R2.fastq","s2.R2.fastq"]]
        
        actual_pairs = utilities.paired_files(files, ".fastq", pair_identifier=".R1")
        
        self.assertEqual(expected_pairs[0],actual_pairs[0])
        self.assertEqual(expected_pairs[1],actual_pairs[1])

    def test_paired_files_period(self):
        """ Test the paired files function with periods in sample name"""
        
        files=["s1.1.R1.fastq","s1.1.R2.fastq","s2.1.R1.fastq","s2.1.R2.fastq"]
        
        expected_pairs = [["s1.1.R1.fastq","s2.1.R1.fastq"],["s1.1.R2.fastq","s2.1.R2.fastq"]]
        
        actual_pairs = utilities.paired_files(files, ".fastq", pair_identifier=".R1")
        
        self.assertEqual(expected_pairs[0],actual_pairs[0])
        self.assertEqual(expected_pairs[1],actual_pairs[1])
        
    def test_sample_names_pair_identifier_duplicate(self):
        """ Test the sample names function with a pair identifier included in the sample name"""
        
        files=["s_1_1.fastq.gz","s1_1.fastq.gz"]
        
        self.assertEqual(utilities.sample_names(files,".fastq.gz","_1"), ["s_1","s1"])
        
    def test_paired_files_duplicate_identifier_1(self):
        """ Test the paired files function with a first identifier duplicated"""
        
        files=["s_1_1.fastq.gz","s_1_2.fastq.gz","s_1_3_1.fastq.gz","s_1_3_2.fastq.gz"]
        
        expected_pairs = [["s_1_1.fastq.gz","s_1_3_1.fastq.gz"],["s_1_2.fastq.gz","s_1_3_2.fastq.gz"]]
        
        actual_pairs = utilities.paired_files(files, ".fastq.gz", pair_identifier="_1")
        
        self.assertEqual(expected_pairs[0],actual_pairs[0])
        self.assertEqual(expected_pairs[1],actual_pairs[1])
        
    def test_paired_files_duplicate_identifier_2(self):
        """ Test the paired files function with a second identifier duplicated
            Also test extension without leading period. """
        
        files=["s_2_1.fastq","s_2_2.fastq","s_2_3_1.fastq","s_2_3_2.fastq"]
        
        expected_pairs = [["s_2_1.fastq","s_2_3_1.fastq"],["s_2_2.fastq","s_2_3_2.fastq"]]
        
        actual_pairs = utilities.paired_files(files, "fastq", pair_identifier="_1")
        
        self.assertEqual(expected_pairs[0],actual_pairs[0])
        self.assertEqual(expected_pairs[1],actual_pairs[1])
        
    def test_paired_files_duplicate_identifier_3(self):
        """ Test the paired files function with pair identified duplicated """
        
        files=["MR100.R1.fastq","MR100.R2.fastq","MR200.R1.fastq","MR200.R2.fastq"]
        
        expected_pairs = [["MR100.R1.fastq","MR200.R1.fastq"],["MR100.R2.fastq","MR200.R2.fastq"]]
        
        actual_pairs = utilities.paired_files(files, "fastq", pair_identifier="R1")
        
        self.assertEqual(expected_pairs[0],actual_pairs[0])
        self.assertEqual(expected_pairs[1],actual_pairs[1])
        
    def test_paired_files_identifier_not_found(self):
        """ Test the paired files function with an identifier that is not found"""
        
        files=["sample-1.R1.fastq","sample-1.R2.fastq","sample-2.R1.fastq","sample-2.R2.fastq"]
        
        expected_pairs = [[],[]]
        
        actual_pairs = utilities.paired_files(files, ".fastq", pair_identifier="_R1.")
        
        self.assertEqual(expected_pairs[0],actual_pairs[0])
        self.assertEqual(expected_pairs[1],actual_pairs[1])
        
    def test_paired_files_identifier_includes_extension(self):
        """ Test the paired files function with an identifier that is not found because
            it includes the period from the file extension"""
        
        files=["sample-1.R1.fastq","sample-1.R2.fastq","sample-2.R1.fastq","sample-2.R2.fastq"]
        
        expected_pairs = [[],[]]
        
        actual_pairs = utilities.paired_files(files, ".fastq", pair_identifier="R1.")
        
        self.assertEqual(expected_pairs[0],actual_pairs[0])
        self.assertEqual(expected_pairs[1],actual_pairs[1])
        
    def test_read_metadata_columns(self):
        """ Test the read metadata function. Metadata file has samples as columns. """
        
        # create a temp taxonomy file (for sample names)
        taxon_file = write_temp("# taxon\tsample1\tsample2\nbug1\t1\t2\n")
        # create a temp metadata file
        metadata_file = write_temp("# feature\tsample1\tsample2\nfeature1\tA\tB\n")
        
        # get the metadata values
        values = utilities.read_metadata(metadata_file,taxon_file)
        
        # remove the temp files
        os.remove(taxon_file)
        os.remove(metadata_file)
        
        expected_values=[["# feature","sample1","sample2"],["feature1","A","B"]]
        self.assertEqual(values,expected_values)
        
    def test_read_metadata_rows(self):
        """ Test the read metadata function. Metadata file has samples as rows. """
        
        # create a temp taxonomy file (for sample names)
        taxon_file = write_temp("# taxon\tsample1\tsample2\nbug1\t1\t2\n")
        # create a temp metadata file
        metadata_file = write_temp("# sample\tfeature1\nsample1\tA\nsample2\tB\n")
        
        # get the metadata values
        values = utilities.read_metadata(metadata_file,taxon_file)
        
        # remove the temp files
        os.remove(taxon_file)
        os.remove(metadata_file)
        
        expected_values=[["# sample","sample1","sample2"],["feature1","A","B"]]
        self.assertEqual(values,expected_values)   
        
    def test_read_metadata_remove_feature(self):
        """ Test the read metadata function. Metadata file has samples as columns.
            Remove a single feature."""
        
        # create a temp taxonomy file (for sample names)
        taxon_file = write_temp("# taxon\tsample1\tsample2\nbug1\t1\t2\n")
        # create a temp metadata file
        metadata_file = write_temp("# feature\tsample1\tsample2\nfeature1\tA\tB\nfeature2\tC\tD\n")
        
        # get the metadata values
        values = utilities.read_metadata(metadata_file,taxon_file,ignore_features=["feature1"])
        
        # remove the temp files
        os.remove(taxon_file)
        os.remove(metadata_file)
        
        expected_values=[["# feature","sample1","sample2"],["feature2","C","D"]]
        self.assertEqual(values,expected_values)  
        
    def test_label_metadata(self):
        """ Test the label_metadata function. Test default labels. """
        
        data=[["# samples","sample1","sample2"],["feature1","1","2.2"],["feature2","A","B"]]
        labels, new_data=utilities.label_metadata(data)   
        
        expected_data=[["# samples","sample1","sample2"],["feature1",1.0,2.2],["feature2","A","B"]]
        expected_labels={"feature1":"con","feature2":"cat"}
        
        self.assertEqual(new_data,expected_data)
        self.assertEqual(labels, expected_labels)
        
    def test_label_metadata_categorical(self):
        """ Test the label_metadata function. Set categorical label. """
        
        data=[["# samples","sample1","sample2"],["feature1","1","2.2"],["feature2","A","B"]]
        labels, new_data=utilities.label_metadata(data, categorical=["feature1"])   
        
        expected_data=[["# samples","sample1","sample2"],["feature1","1","2.2"],["feature2","A","B"]]
        expected_labels={"feature1":"cat","feature2":"cat"}
        
        self.assertEqual(new_data,expected_data)
        self.assertEqual(labels, expected_labels)
        
    def test_label_metadata_continuous(self):
        """ Test the label_metadata function. Set continuous label. """
        
        data=[["# samples","sample1","sample2"],["feature1","1","2.2"],["feature2","A","B"]]
        labels, new_data=utilities.label_metadata(data, continuous=["feature2"])   
        
        expected_data=[["# samples","sample1","sample2"],["feature1",1.0,2.2],["feature2","A","B"]]
        expected_labels={"feature1":"con","feature2":"con"}
        
        self.assertEqual(new_data,expected_data)
        self.assertEqual(labels, expected_labels)
        
    def test_merge_metadata(self):
        """ Test the merge metadata function."""
        
        metadata=[["# samples","s1","s2"],["feature1","A","B"],["feature2",1,2]]
        samples=["s1","s2"]
        values=[["bug1",1,2],["bug2",2,4]]
        
        merged, new_samples=utilities.merge_metadata(metadata, samples, values)
        
        expected=[["feature1","A","B"],["feature2",1,2],["bug1",1,2],["bug2",2,4]]
        expected_samples=["s1","s2"]
        
        self.assertEqual(merged, expected) 
        self.assertEqual(new_samples, expected_samples)    
        
    def test_merge_metadata_values_without_names(self):
        """ Test the merge metadata function. Test values without names option."""
        
        metadata=[["# samples","s1","s2"],["feature1","A","B"],["feature2",1,2]]
        samples=["s1","s2"]
        values=[[1,2],[2,4]]
        
        merged, new_samples=utilities.merge_metadata(metadata, samples, values, values_without_names=True)
        
        expected=[["feature1","A","B"],["feature2",1,2],[1,2],[2,4]]
        expected_samples=["s1","s2"]
        
        self.assertEqual(merged, expected) 
        self.assertEqual(new_samples, expected_samples)  
        
    def test_merge_metadata_reorder(self):
        """ Test the merge metadata function. Test with reordering."""
        
        metadata=[["# samples","s1","s2"],["feature1","A","B"],["feature2",1,2]]
        samples=["s2","s1"]
        values=[["bug1",1,2],["bug2",2,4]]
        
        merged, new_samples=utilities.merge_metadata(metadata, samples, values)
        
        expected=[["feature1","A","B"],["feature2",1,2],["bug1",2,1],["bug2",4,2]]
        expected_samples=["s1","s2"]
        
        self.assertEqual(merged, expected) 
        self.assertEqual(new_samples, expected_samples) 
        
    def test_merge_metadata_less_metadata(self):
        """ Test the merge metadata function. Test with only a single sample metadata."""
        
        metadata=[["# samples","s1"],["feature1","A"],["feature2",1]]
        samples=["s1","s2"]
        values=[["bug1",1,2],["bug2",2,4]]
        
        # suppress warning message
        # Redirect stdout
        original_stdout=sys.stdout
        sys.stdout=open(os.devnull,"w")
        
        merged, new_samples=utilities.merge_metadata(metadata, samples, values)
        
        # Redirect stdout
        sys.stdout=original_stdout
        
        expected=[["feature1","A"],["feature2",1],["bug1",1],["bug2",2]]
        expected_samples=["s1"]
        
        self.assertEqual(merged, expected)  
        self.assertEqual(new_samples, expected_samples)   
        
    def test_merge_metadata_nomatch(self):
        """ Test the merge metadata function. Test with metadata not matching."""
        
        metadata=[["# samples","s3"],["feature1","A"],["feature2",1]]
        samples=["s1","s2"]
        values=[["bug1",1,2],["bug2",2,4]]
        
        # suppress warning message
        # Redirect stdout
        original_stdout=sys.stdout
        sys.stdout=open(os.devnull,"w")
        
        merged, new_samples=utilities.merge_metadata(metadata, samples, values)
        
        # Redirect stdout
        sys.stdout=original_stdout
        
        expected=[["bug1",1,2],["bug2",2,4]]
        
        self.assertEqual(merged, expected) 
        self.assertEqual(new_samples, samples)  
        
    def test_group_samples_by_metadata(self):
        """ Test the group samples by metadata function """
        
        metadata=["feature1","A","B","A","A","B"]
        data=[[1,2,3,4,5],[11,22,33,44,55],[111,222,333,444,555]]
        samples=["s1","s2","s3","s4","s5"]
        
        sorted_data_grouped, sorted_samples_grouped = utilities.group_samples_by_metadata(metadata, data, samples)
        
        expected_samples={"A":["s1","s3","s4"], "B":["s2","s5"]}
        expected_data={"A":[[1,3,4],[11,33,44],[111,333,444]], "B":[[2,5],[22,55],[222,555]]}
        
        self.assertEqual(sorted_data_grouped, expected_data)
        self.assertEqual(sorted_samples_grouped, expected_samples)
        
    def test_read_picard(self):
        """ Test the read picard function. Test with all values passing threshold """
        
        file_text="\n".join(["## htsjdk.samtools.metrics.StringHeader",
            "# picard.analysis.CollectMultipleMetrics INPUT=AYBNT.1.unmapped.bam ASSUME_SORTED=true",
            "# Started on: Thu Jan 01 01:12:33 EST 2017",
            "",
            "",
            "## HISTOGRAM    java.lang.Integer",
            "CYCLE\tMEAN_QUALITY",
            "1\t39.111874",
            "2\t35.614163",
            "3\t35.887343",
            "4\t42.882504",
            "5\t52.643778"])
        
        expected_data=[(1,39.111874),(2,35.614163),(3,35.887343),(4,42.882504),(5,52.643778)]
        
        tmp_file=write_temp(file_text)
        
        data, below_threshold = utilities.read_picard(tmp_file)
        os.remove(tmp_file)
        
        self.assertFalse(below_threshold)
        self.assertEqual(data, expected_data)
        
    def test_read_picard_threshold(self):
        """ Test the read picard function. Test with some values passing threshold """
        
        file_text="\n".join(["## htsjdk.samtools.metrics.StringHeader",
            "# picard.analysis.CollectMultipleMetrics INPUT=AYBNT.1.unmapped.bam ASSUME_SORTED=true",
            "# Started on: Thu Jan 01 01:12:33 EST 2017",
            "",
            "",
            "## HISTOGRAM    java.lang.Integer",
            "CYCLE\tMEAN_QUALITY",
            "1\t39.111874",
            "2\t35.614163",
            "3\t35.887343",
            "4\t42.882504",
            "5\t52.643778"])
        
        expected_data=[(1,39.111874),(2,35.614163),(3,35.887343),(4,42.882504),(5,52.643778)]
        
        tmp_file=write_temp(file_text)
        
        data, below_threshold = utilities.read_picard(tmp_file, threshold=37)
        os.remove(tmp_file)
        
        self.assertTrue(below_threshold)
        self.assertEqual(data, expected_data)
        
    def test_rank_species_average_abundance(self):
        """ Test the rank species average abundance function with a small merged taxonomy file """
        
        file_contents = "\n".join(["#SampleID\tS1\tS2\tS3",
        "k__Archaea\t0.5\t0.5\t0.5", "k__Archaea|p__Euryarchaeota\t0.5\t0.5\t0.5",
        "k__Viruses|p__Viruses_noname|c__Viruses_noname|o__Viruses_noname|f__Retroviridae|g__Gammaretrovirus|s__Murine_osteosarcoma_virus\t0.1\t0.2\t0.3",
        "k__Viruses|p__Viruses_noname|c__Viruses_noname|o__Viruses_noname|f__Retroviridae|g__Gammaretrovirus|s__Murine_osteosarcoma_virus|t__PRJNA14655\t0.2\t0.2\t0.2",
        "g__Alistipes|s__Alistipes_indistinctus\t0.3\t0.4\t0.5",
        "g__Alistipes|s__Alistipes_onderdonkii\t0.4\t0.5\t0.6",
        "g__Alistipes|s__Alistipes_putredinis\t0.1\t0.2\t0.1"])

        expected_species = ["s__Alistipes_onderdonkii","s__Alistipes_indistinctus","s__Murine_osteosarcoma_virus","s__Alistipes_putredinis"]

        temp_file = write_temp(file_contents)
        
        ranked_species = utilities.rank_species_average_abundance(temp_file)
        os.remove(temp_file)
        
        self.assertEqual(ranked_species,expected_species)
        
    def test_order_clade_list(self):
        """ Test the order clade list function """
        
        taxonomy_file_contents = "\n".join(["#SampleID\tS1\tS2\tS3",
        "k__Archaea\t0.5\t0.5\t0.5", "k__Archaea|p__Euryarchaeota\t0.5\t0.5\t0.5",
        "k__Viruses|p__Viruses_noname|c__Viruses_noname|o__Viruses_noname|f__Retroviridae|g__Gammaretrovirus|s__Murine_osteosarcoma_virus\t0.1\t0.2\t0.3",
        "k__Viruses|p__Viruses_noname|c__Viruses_noname|o__Viruses_noname|f__Retroviridae|g__Gammaretrovirus|s__Murine_osteosarcoma_virus|t__PRJNA14655\t0.2\t0.2\t0.2",
        "g__Alistipes|s__Alistipes_indistinctus\t0.3\t0.7\t0.5",
        "g__Alistipes|s__Alistipes_onderdonkii\t0.4\t0.5\t0.6",
        "g__Alistipes|s__Alistipes_putredinis\t0.1\t0.1\t0.1"])

        clades_list_file_contents = "\n".join(["s__Alistipes_putredinis", 
            "s__Alistipes_indistinctus", "s__Murine_osteosarcoma_virus", "s__Acidaminococcus_sp_D21"])

        expected_species = ["s__Alistipes_indistinctus","s__Murine_osteosarcoma_virus","s__Alistipes_putredinis"]

        taxonomy_temp_file = write_temp(taxonomy_file_contents)
        clade_temp_file = write_temp(clades_list_file_contents)
        output_file = write_temp("")
        
        ranked_species = utilities.order_clade_list(None,clade_temp_file,taxonomy_temp_file,output_file)
        os.remove(taxonomy_temp_file)
        os.remove(clade_temp_file)
        
        with open(output_file) as file_handle:
            ranked_species=[line.rstrip() for line in file_handle.readlines()]
        os.remove(output_file)
        
        self.assertEqual(ranked_species,expected_species)        
        
        

