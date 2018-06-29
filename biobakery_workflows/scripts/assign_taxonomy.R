#!/usr/bin/env Rscript

# load packages
library(dada2); packageVersion("dada2")
#library(ggplot2)
#library(msa)
#library(gridExtra)
#library(phangorn)



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
pool.samples <- args.list$pool
output.path <- normalizePath( args.list$output_dir )

seqtab.nochim <- readRDS(paste0(output.path, "/seqtab_final.rds"))

# Assign taxonomy:
taxa.gg13_8 <- assignTaxonomy(seqtab.nochim, paste0(output.path,"/../dada2_reference_databases/gg_13_8_train_set_97.fa.gz"), multithread=TRUE, tryRC=TRUE)

# Print first 6 rows of taxonomic assignment
unname(head(taxa.gg13_8))

# Replace NAs in taxonomy assignment table with prefix corresponding to tax rank
taxa.gg13_8.2 <- replaceNA.in.assignedTaxonomy( taxa.gg13_8 )

# Write taxa table to file
write.table( taxa.gg13_8.2, paste0( output.path, "/all_samples_GG13-8-taxonomy.tsv" ), sep = "\t", eol = "\n", quote = F, col.names = NA )

## Merge OTU and GG13-8 taxonomy tables
otu.gg.tax.table <- merge( t(seqtab.nochim), taxa.gg13_8.2, by = 'row.names' )
rownames( otu.gg.tax.table ) <- otu.gg.tax.table[,1]
otu.gg.tax.table <- otu.gg.tax.table[,-1]

write.table(otu.gg.tax.table, paste0( output.path, "/all_samples_SV-counts_and_GG13-8-taxonomy.tsv" ), sep = "\t", eol = "\n", quote = F, col.names = NA)
