#!/usr/bin/env python

from rpy2 import robjects
from rpy2.robjects.packages import importr
import sys
import argparse

base = importr('base')
dada2 = importr('dada2')

def parse_arguments(args):
    """ Parse the arguments"""

    parser = argparse.ArgumentParser(
        description="Assign taxonomy arguments\n",
        formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument(
        "--seqtab_rds",
        help="sequence table rds file\n[REQUIRED]",
    )
    parser.add_argument(
        "--refdb_path",
        help="path to reference db\n[REQUIRED]",
    )
    parser.add_argument(
        "--refdb_species_path",
        help="path to reference species db\n[REQUIRED]",
    )
    parser.add_argument(
        "--otu_closed_ref_path",
        help = "path to otu closed reference file\n[REQUIRED]",
    )
    parser.add_argument(
        "--threads",
        help="number of threads\n[REQUIRED]",
    )

    return parser.parse_args()


def main():

    # parse arguments from the user
    args = parse_arguments(sys.argv)

    robjects.globalenv["seqtab_rds"] = args.seqtab_rds
    robjects.globalenv["refdb_species_path"] = args.refdb_species_path
    robjects.globalenv["refdb_path"] = args.refdb_path
    robjects.globalenv["otu_closed_ref_path"] = args.otu_closed_ref_path
    robjects.globalenv["threads"] = args.threads

    robjects.r('''
    
    
    # Helper function to replace NAs in taxonomy assignment with space
    removeNA.in.assignedTaxonomy <- 
      function( tax.table ) {
        prefix <- c( ' ', ' ', ' ', ' ', ' ', ' ', ' ' )
    
        for( i in 1 : length( colnames( tax.table ) ) ) {
          tax.table[ ,i ] <- 
            ifelse( is.na(tax.table[ ,i ] ), 
                    prefix[i],
                    tax.table[ ,i ]
            )
        }
        rm(i)
        return( tax.table )
      }
    
    # Helper function to replace NAs in taxonomy assignment with prefix
    replaceNA.in.assignedTaxonomy <- 
      function( tax.table ) {
        prefix <- c( 'k__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__' )
    
        for( i in 1 : length( colnames( tax.table ) ) ) {
          tax.table[ ,i ] <- 
            ifelse( is.na(tax.table[ ,i ] ), 
                    prefix[i],
                    tax.table[ ,i ]
            )
        }
        rm(i)
        return( tax.table )
      }
    
    
    refdb.path <- normalizePath(refdb_path )
    seqtab.nochim <- readRDS(seqtab_rds)
    
    ## Asign SILVA or  RDP taxonomies and merge with OTU table
    taxa.refdb <- dada2::assignTaxonomy(seqtab.nochim, refdb.path, multithread=as.numeric(threads))
    
    if (!identical(refdb_species_path, "None")) {
    refdb.species.path <- normalizePath(refdb_species_path )
    # Append species. Note that appending the argument 'allowMultiple=3' will return up to 3 different matched
    # species, but if 4 or more are matched it returns NA.
    taxa.refdb.species <- addSpecies(taxa.refdb, refdb.species.path)
    
    # Replace NAs in taxonomy assignment table with prefix corresponding to tax rank
    taxa.refdb.species.2 <- removeNA.in.assignedTaxonomy(taxa.refdb.species )
    } else {
    # No need to add species, just replacing NAs
    taxa.refdb.species.2 <- replaceNA.in.assignedTaxonomy(taxa.refdb)
    }
    
    # Print first 6 rows of taxonomic assignment
    unname(head(taxa.refdb.species.2))
    
    # Merge with OTU table and save to file
    otu.refdb.tax.table <- merge(t(seqtab.nochim), taxa.refdb.species.2, by = 'row.names' )
    rownames(otu.refdb.tax.table) <- otu.refdb.tax.table[, 1]
    otu.refdb.tax.table <- otu.refdb.tax.table[, -1]
    
    otu.refdb.tax.table_taxcombined <- cbind(otu.refdb.tax.table)
    colnum <- length(otu.refdb.tax.table_taxcombined[1,])
    
    otu.refdb.tax.table_taxcombined <- otu.refdb.tax.table_taxcombined[, -c((colnum - 6): colnum)]
    
    taxonomy <- vector()
    
    if (!identical(refdb_species_path, "None")) {
    taxonomy <- paste0("k__", as.character(otu.refdb.tax.table$Kingdom), "; ",
    "p__", as.character(otu.refdb.tax.table$Phylum), "; ",
    "c__", as.character(otu.refdb.tax.table$Class), "; ",
    "o__", as.character(otu.refdb.tax.table$Order), "; ",
    "f__", as.character(otu.refdb.tax.table$Family), "; ",
    "g__", as.character(otu.refdb.tax.table$Genus), "; ",
    "s__", as.character(otu.refdb.tax.table$Species), "; ")
    } else {
    taxonomy <- paste0(as.character(otu.refdb.tax.table$Kingdom), "; ",
    as.character(otu.refdb.tax.table$Phylum), "; ",
    as.character(otu.refdb.tax.table$Class), "; ",
    as.character(otu.refdb.tax.table$Order), "; ",
    as.character(otu.refdb.tax.table$Family), "; ",
    as.character(otu.refdb.tax.table$Genus), "; ",
    as.character(otu.refdb.tax.table$Species))
    }
    
    otu.refdb.tax.table_taxcombined <- cbind(otu.refdb.tax.table_taxcombined, taxonomy)
    
    write.table(otu.refdb.tax.table_taxcombined, otu_closed_ref_path,
        sep = "\t", eol = "\n", quote = F, col.names = NA)
    write.table(otu.refdb.tax.table, paste0(gsub(".tsv", "", otu_closed_ref_path), "_taxcolumns.tsv"),
        sep = "\t", eol = "\n", quote = F, col.names = NA)
    
    ''')


if __name__ == "__main__":
        main()