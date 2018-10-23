#!/usr/bin/env Rscript

# Load packages
library(data.table)
library(dada2); packageVersion("dada2")


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
    prefix <- c( 'k__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__', 'otuid__' )
    
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

# These variables are passed to the workflow
output.path <- normalizePath( args.list$output_dir )
refdb.path <- normalizePath( args.list$refdb_path )

# Read ASV table from saved file
seqtab.nochim <- readRDS(args.list$seqtab_file_path)

if (!identical(args.list$refdb_species_path,"None")) {
  
  ## Asign  SILVA or  RDP taxonomies 
  taxa.refdb <- dada2::assignTaxonomy(seqtab.nochim, refdb.path, multithread = as.numeric(args.list$threads))
 
  ## Append species if reference db is not green genes
  refdb.species.path <- normalizePath( args.list$refdb_species_path )
  taxa.refdb.species <- addSpecies(taxa.refdb, refdb.species.path)
 
  # Remove NAs in taxonomy assignment
  taxa.refdb.species.2 <- removeNA.in.assignedTaxonomy(taxa.refdb.species )
} else {
  
  ## Asign GreenGenes taxonomies with OTU IDs
   taxa.refdb <- dada2::assignTaxonomy(seqtab.nochim, refdb.path,
                                       minBoot = 0,
                                       outputBootstrap = TRUE,
                                       taxLevels = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "OTUID"),
                                       verbose = TRUE,
                                       multithread = as.numeric(args.list$threads))
    
  
  # No need to add species if green genes, just replacing NAs with prefix
  taxa.refdb.species <- taxa.refdb$tax
  taxa.refdb.species.2 <- as.data.frame(replaceNA.in.assignedTaxonomy(taxa.refdb.species))
  
  # Merge taxonomic assignments, bootstrap info, and counts in samples
  otu.refdb.tax.boot <- merge( t(seqtab.nochim), taxa.refdb, by = 'row.names' )
  write.table(otu.refdb.tax.boot, paste0(gsub(".tsv", "", args.list$otu_closed_ref_path),"_withbootstrap.tsv") , sep = "\t", eol = "\n", quote = F, row.names=FALSE)
  
}

# Merge taxa with counts in samples for ASV/OTU table
otu.refdb.tax.table <- merge( t(seqtab.nochim), taxa.refdb.species, by = 'row.names' )
rownames( otu.refdb.tax.table ) <- otu.refdb.tax.table[,1]
otu.refdb.tax.table <- otu.refdb.tax.table[,-1]

otu.refdb.tax.table_taxcombined <- cbind(otu.refdb.tax.table)
colnum <- length(otu.refdb.tax.table_taxcombined[1,])

# Combining taxonomy levels into one column 
if (!identical(args.list$refdb_species_path,"None")) {
  otu.refdb.tax.table_taxcombined <- otu.refdb.tax.table_taxcombined[, -c((colnum-6): colnum)]
}else{
  otu.refdb.tax.table_taxcombined <- otu.refdb.tax.table_taxcombined[, -c((colnum-7): (colnum-1))]  
}

# Name combined column taxonomy
taxonomy <- vector()
 
# Concatenate taxonomy levels separting by '; ' and write into taxonomy column. Append prefix if db is green genes.
if (!identical(args.list$refdb_species_path,"None")) {
  taxonomy<- paste0("k__",as.character(otu.refdb.tax.table$Kingdom),"; ",
                  "p__",as.character(otu.refdb.tax.table$Phylum),"; ",
                  "c__",as.character(otu.refdb.tax.table$Class),"; ",
                  "o__",as.character(otu.refdb.tax.table$Order),"; ",
                  "f__",as.character(otu.refdb.tax.table$Family),"; ",
                  "g__",as.character(otu.refdb.tax.table$Genus),"; ",
                  "s__",as.character(otu.refdb.tax.table$Species))
} else{
  taxonomy<- paste0(as.character(otu.refdb.tax.table$Kingdom),"; ",
                    as.character(otu.refdb.tax.table$Phylum),"; ",
                    as.character(otu.refdb.tax.table$Class),"; ",
                    as.character(otu.refdb.tax.table$Order),"; ",
                    as.character(otu.refdb.tax.table$Family),"; ",
                    as.character(otu.refdb.tax.table$Genus),"; ",
                    as.character(otu.refdb.tax.table$Species))
}

# Add taxonomy column
otu.refdb.tax.table_taxcombined <- cbind(otu.refdb.tax.table_taxcombined,taxonomy)

if (identical(args.list$refdb_species_path,"None")) {
  
  otuids <- gsub("otuid__","",otu.refdb.tax.table_taxcombined$OTUID)
  otu.refdb.tax.table_taxcombined <- otu.refdb.tax.table_taxcombined[ ,names(otu.refdb.tax.table_taxcombined) != "OTUID"]
  otu.refdb.tax.table_taxcombined <- cbind(otuids,otu.refdb.tax.table_taxcombined)
  otu.refdb.tax.table_taxcombined_uniqueids <- as.data.frame(data.table(otu.refdb.tax.table_taxcombined)[, lapply(.SD, sum), keyby = .(otuids, taxonomy)])
  taxonomy <- otu.refdb.tax.table_taxcombined_uniqueids$taxonomy
  otu.refdb.tax.table_taxcombined_uniqueids <- otu.refdb.tax.table_taxcombined_uniqueids[ ,names(otu.refdb.tax.table_taxcombined_uniqueids) != "taxonomy"]
  otu.refdb.tax.table_taxcombined_uniqueids <-  cbind(otu.refdb.tax.table_taxcombined_uniqueids, taxonomy)
 
  # Save closed reference with OTU ids to tsv file
  write.table(otu.refdb.tax.table_taxcombined, paste0(gsub(".tsv", "", args.list$otu_closed_ref_path),"_nosum.tsv") , sep = "\t", eol = "\n", quote = F, row.names=FALSE)
  write.table(otu.refdb.tax.table_taxcombined_uniqueids, args.list$otu_closed_ref_path , sep = "\t", eol = "\n", quote = F, row.names=FALSE)
 
}else{
  
  # Save closed reference to tsv file both versions: single column and mutiple taxonomy levels columns
  write.table(otu.refdb.tax.table_taxcombined, paste0(gsub(".tsv", "", args.list$otu_closed_ref_path),"_withseqs.tsv") , sep = "\t", eol = "\n", quote = F, col.names = NA)
  write.table(otu.refdb.tax.table, paste0(gsub(".tsv", "", args.list$otu_closed_ref_path),"_taxcolumns.tsv") , sep = "\t", eol = "\n", quote = F, col.names = NA)
  

  # Create version with ASV1, ASV2 ... ids instead of sequences as ids
  seqids <- c(1:length(otu.refdb.tax.table_taxcombined[,1]))
  seqids <- paste0("ASV",seqids)
  row.names(otu.refdb.tax.table_taxcombined) <- seqids 

  # Save closed reference with ASV ids to tsv file
  write.table(otu.refdb.tax.table_taxcombined, args.list$otu_closed_ref_path , sep = "\t", eol = "\n", quote = F, col.names = NA)
}