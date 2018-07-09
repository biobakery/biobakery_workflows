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
#pool.samples <- args.list$pool
output.path <- normalizePath( args.list$output_dir )

seqtab.nochim <- readRDS(paste0(output.path, "/seqtab_final.rds"))

## Asign SILVA and RDP taxonomies and merge with OTU table

# Assign SILVA taxonomy
taxa.silva <- assignTaxonomy(seqtab.nochim, paste0(output.path,"/../dada2_reference_databases/silva_nr_v128_train_set.fa.gz"), multithread = TRUE)

# Replace NAs in taxonomy assignment table with prefix corresponding to tax rank
taxa.silva.2 <- replaceNA.in.assignedTaxonomy( taxa.silva )

# OMIT APPENDING SPECIES FOR SILVA DUE TO MEMORY CONSTRAINTS
# Append species. Note that appending the argument 'allowMultiple=3' will return up to 3 different matched
# species, but if 4 or more are matched it returns NA.
#taxa.silva.species <- addSpecies(taxa.silva, paste0(output.path,"/../dada2_reference_databases/silva_species_assignment_v128.fa.gz"))

# Merge with OTU table and save to file
otu.silva.tax.table <- merge( t(seqtab.nochim), taxa.silva.2, by = 'row.names' )
rownames( otu.silva.tax.table ) <- otu.silva.tax.table[,1]
otu.silva.tax.table <- otu.silva.tax.table[,-1]

write.table(otu.silva.tax.table, paste0(output.path, "/all_samples_taxonomy_closed_reference_silva.tsv" ), sep = "\t", eol = "\n", quote = F, col.names = NA)

# Assign RDP taxonomy
taxa.rdp <- assignTaxonomy(seqtab.nochim, paste0(output.path,"/../dada2_reference_databases/rdp_train_set_16.fa.gz"), multithread = TRUE)

# Replace NAs in taxonomy assignment table with prefix corresponding to tax rank
taxa.rdp.2 <- replaceNA.in.assignedTaxonomy( taxa.rdp )

# OMIT APPENDING SPECIES FOR RDP DUE TO MEMORY CONSTRAINTS
# Append species. Note that appending the argument 'allowMultiple=3' will return up to 3 different matched
# species, but if 4 or more are matched it returns NA.
#taxa.rdp.species <- addSpecies(taxa.rdp, paste0(output.path,"/../dada2_reference_databases/rdp_species_assignment_16.fa.gz"))


# Merge with OTU table and save to file
otu.rdp.tax.table <- merge( t(seqtab.nochim), taxa.rdp.2, by = 'row.names' )
rownames( otu.rdp.tax.table ) <- otu.rdp.tax.table[,1]
otu.rdp.tax.table <- otu.rdp.tax.table[,-1]

write.table(otu.rdp.tax.table, paste0( output.path, "/all_samples_taxonomy_closed_reference_rdp.tsv" ), sep = "\t", eol = "\n", quote = F, col.names = NA)