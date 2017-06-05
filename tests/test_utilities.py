
import unittest

from biobakery_workflows import utilities

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
        
        files=["s1.fastq","s1.fastq.gz"]
        
        self.assertEqual(utilities.sample_names(files), ["s1","s1"])
        
    def test_sample_names_pair_identifier(self):
        """ Test the sample names function with a pair identifier """
        
        files=["s1.R1.fastq","s1.R1.fastq.gz"]
        
        self.assertEqual(utilities.sample_names(files,".R1"), ["s1","s1"])
        
    def test_sample_names_pair_identifier_not_found(self):
        """ Test the sample names function with a pair identifier that is not included in the names """
        
        files=["s1.R1.fastq","s1.R1.fastq.gz"]
        
        self.assertEqual(utilities.sample_names(files,"_R1"), ["s1.R1","s1.R1"])
        
    def test_paired_files(self):
        """ Test the paired files function """
        
        files=["s1.R1.fastq","s1.R2.fastq","s2.R1.fastq.gz","s2.R2.fastq.gz"]
        
        expected_pairs = [["s1.R1.fastq","s2.R1.fastq.gz"],["s1.R2.fastq","s2.R2.fastq.gz"]]
        
        actual_pairs = utilities.paired_files(files, pair_identifier=".R1")
        
        self.assertEqual(expected_pairs[0],actual_pairs[0])
        self.assertEqual(expected_pairs[1],actual_pairs[1])
        
    def test_sample_names_pair_identifier_duplicate(self):
        """ Test the sample names function with a pair identifier included in the sample name"""
        
        files=["s_1_1.fastq","s1_1.fastq.gz"]
        
        self.assertEqual(utilities.sample_names(files,"_1"), ["s_1","s1"])
        
    def test_paired_files_duplicate_identifier_1(self):
        """ Test the paired files function with a first identifier duplicated"""
        
        files=["s_1_1.fastq","s_1_2.fastq","s_1_3_1.fastq.gz","s_1_3_2.fastq.gz"]
        
        expected_pairs = [["s_1_1.fastq","s_1_3_1.fastq.gz"],["s_1_2.fastq","s_1_3_2.fastq.gz"]]
        
        actual_pairs = utilities.paired_files(files, pair_identifier="_1")
        
        self.assertEqual(expected_pairs[0],actual_pairs[0])
        self.assertEqual(expected_pairs[1],actual_pairs[1])
        
    def test_paired_files_duplicate_identifier_2(self):
        """ Test the paired files function with a second identifier duplicated"""
        
        files=["s_2_1.fastq","s_2_2.fastq","s_2_3_1.fastq.gz","s_2_3_2.fastq.gz"]
        
        expected_pairs = [["s_2_1.fastq","s_2_3_1.fastq.gz"],["s_2_2.fastq","s_2_3_2.fastq.gz"]]
        
        actual_pairs = utilities.paired_files(files, pair_identifier="_1")
        
        self.assertEqual(expected_pairs[0],actual_pairs[0])
        self.assertEqual(expected_pairs[1],actual_pairs[1])
        
    def test_paired_files_identifier_not_found(self):
        """ Test the paired files function with an identifier that is not found"""
        
        files=["sample-1.R1.fastq","sample-1.R2.fastq","sample-2.R1.fastq.gz","sample-2.R2.fastq.gz"]
        
        expected_pairs = [[],[]]
        
        actual_pairs = utilities.paired_files(files, pair_identifier="_R1.")
        
        self.assertEqual(expected_pairs[0],actual_pairs[0])
        self.assertEqual(expected_pairs[1],actual_pairs[1])
        
    def test_paired_files_identifier_includes_extension(self):
        """ Test the paired files function with an identifier that is not found because
            it includes the period from the file extension"""
        
        files=["sample-1.R1.fastq","sample-1.R2.fastq","sample-2.R1.fastq.gz","sample-2.R2.fastq.gz"]
        
        expected_pairs = [[],[]]
        
        actual_pairs = utilities.paired_files(files, pair_identifier="R1.")
        
        self.assertEqual(expected_pairs[0],actual_pairs[0])
        self.assertEqual(expected_pairs[1],actual_pairs[1])
        

