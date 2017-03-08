
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
        
        
        
