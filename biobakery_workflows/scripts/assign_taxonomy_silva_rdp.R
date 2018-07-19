#!/usr/bin/env Rscript

# load packages
library(dada2); packageVersion("dada2")


# Helper function to replace NAs in taxonomy assignment table with prefix corresponding to tax rank
replaceNA.in.assignedTaxonomy <- 
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


## Collect arguments
args <- commandArgs(TRUE)

## Parse arguments (we expect the form --arg=value)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
args.df <- as.data.frame(do.call("rbind", parseArgs(args)))

args.list <- as.list(as.character(args.df$V2))
names(args.list) <- args.df$V1

## Arg1 default
if(is.null(args.list$output_dir)) {
 stop("At least one argument must be supplied (output folder).\n", call.=FALSE)
}

# Print args list to STDOUT
for( i in names(args.list) ) {
  cat( i, "\t", args.list[[i]], "\n")
}


output.dir <- ifelse( is.null(args.list$output_dir), "output", args.list$output_dir )
output.path <- normalizePath( args.list$output_dir )
refdb.path <- normalizePath( args.list$refdb_path )
refdb.species.path <- normalizePath( args.list$refdb_species_path )

seqtab.nochim <- readRDS(args.list$seqtab_file_path)

## Asign SILVA or  RDP taxonomies and merge with OTU table
taxa.refdb <- dada2::assignTaxonomy(seqtab.nochim, refdb.path, multithread = TRUE)

# Print first 6 rows of taxonomic assignment
unname(head(taxa.refdb))

# Append species. Note that appending the argument 'allowMultiple=3' will return up to 3 different matched
# species, but if 4 or more are matched it returns NA.
taxa.refdb.species <- addSpecies(taxa.refdb, refdb.species.path)

# Replace NAs in taxonomy assignment table with prefix corresponding to tax rank
taxa.refdb.species.2 <- replaceNA.in.assignedTaxonomy(taxa.refdb.species )

# Merge with OTU table and save to file
otu.refdb.tax.table <- merge( t(seqtab.nochim), taxa.refdb.species.2, by = 'row.names' )
rownames( otu.refdb.tax.table ) <- otu.refdb.tax.table[,1]
otu.refdb.tax.table <- otu.refdb.tax.table[,-1]

otu.refdb.tax.table_taxcombined <- cbind(otu.refdb.tax.table)
colnum <- length(otu.refdb.tax.table_taxcombined[1,])

otu.refdb.tax.table_taxcombined <- otu.refdb.tax.table_taxcombined[, -c((colnum-6): colnum)]

taxonomy <- vector()
taxonomy<- paste0("k__",as.character(otu.refdb.tax.table$Kingdom),"; ",
                  "p__",as.character(otu.refdb.tax.table$Phylum),"; ",
                  "c__",as.character(otu.refdb.tax.table$Class),"; ",
                  "o__",as.character(otu.refdb.tax.table$Order),"; ",
                  "f__",as.character(otu.refdb.tax.table$Family),"; ",
                  "g__",as.character(otu.refdb.tax.table$Genus),"; ",
                  "s__",as.character(otu.refdb.tax.table$Species),"; ")


otu.refdb.tax.table_taxcombined <- cbind(otu.refdb.tax.table_taxcombined,taxonomy)


write.table(otu.refdb.tax.table_taxcombined, args.list$otu_closed_ref_path , sep = "\t", eol = "\n", quote = F, col.names = NA)
write.table(otu.refdb.tax.table, paste0(gsub(".tsv", "", args.list$otu_closed_ref_path),"_taxcolumns.tsv") , sep = "\t", eol = "\n", quote = F, col.names = NA)

