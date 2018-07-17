#!/usr/bin/env Rscript

# load packages
library(dada2); packageVersion("dada2")

# Helper function to replace NAs in taxonomy assignment table with prefix corresponding to tax rank
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
gg.path <- normalizePath( args.list$refdb_path )

seqtab.nochim <- readRDS(args.list$seqtab_file_path)

# Assign taxonomy:
taxa.gg13_8 <- dada2::assignTaxonomy(seqtab.nochim, gg.path, multithread=TRUE, tryRC=TRUE)

# Print first 6 rows of taxonomic assignment
unname(head(taxa.gg13_8))

# Replace NAs in taxonomy assignment table with prefix corresponding to tax rank
taxa.gg13_8.2 <- replaceNA.in.assignedTaxonomy(taxa.gg13_8 )

# Write taxa table to file
write.table( taxa.gg13_8.2, args.list$allsamples_gg_tax_path, sep = "\t", eol = "\n", quote = F, col.names = NA )

## Merge OTU and GG13-8 taxonomy tables
otu.gg.tax.table <- merge( t(seqtab.nochim), taxa.gg13_8.2, by = 'row.names' )
rownames( otu.gg.tax.table ) <- otu.gg.tax.table[,1]
otu.gg.tax.table <- otu.gg.tax.table[,-1]

otu.gg.tax.table_taxcombined <- cbind(otu.gg.tax.table)
colnum <- length(otu.gg.tax.table_taxcombined[1,])

otu.gg.tax.table_taxcombined <- otu.gg.tax.table_taxcombined[, -c((colnum-6): colnum)]

taxonomy <- vector()
taxonomy<- paste0(as.character(otu.gg.tax.table$Kingdom),"; ",
                 as.character(otu.gg.tax.table$Phylum),"; ",
                 as.character(otu.gg.tax.table$Class),"; ",
                 as.character(otu.gg.tax.table$Order),"; ",
                 as.character(otu.gg.tax.table$Family),"; ",
                 as.character(otu.gg.tax.table$Genus),"; ",
                 as.character(otu.gg.tax.table$Species))


otu.gg.tax.table_taxcombined <- cbind(otu.gg.tax.table_taxcombined,taxonomy)

write.table(otu.gg.tax.table, paste0(gsub(".tsv", "", args.list$otu_closed_ref_path),"_taxcolumns.tsv"), sep = "\t", eol = "\n", quote = F, col.names = NA)
write.table(otu.gg.tax.table_taxcombined, args.list$otu_closed_ref_path, sep = "\t", eol = "\n", quote = F, col.names = NA)

