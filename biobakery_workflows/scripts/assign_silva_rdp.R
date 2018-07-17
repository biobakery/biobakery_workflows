#!/usr/bin/env Rscript

# load packages
library(dada2); packageVersion("dada2")


## Collect arguments
args <- commandArgs(TRUE)

## Parse arguments (we expect the form --arg=value)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
args.df <- as.data.frame(do.call("rbind", parseArgs(args)))

args.list <- as.list(as.character(args.df$V2))
names(args.list) <- args.df$V1

## Arg1 default
if(is.null(args.list$output_dir)) {
#  stop("At least one argument must be supplied (output folder).\n", call.=FALSE)
  args.list$output_dir = "/Users/anamailyan/outtest"  
}

# Print args list to STDOUT
for( i in names(args.list) ) {
  cat( i, "\t", args.list[[i]], "\n")
}

output.dir <- ifelse( is.null(args.list$output_dir), "output", args.list$output_dir )
output.path <- normalizePath( args.list$output_dir )
rdp.path <- normalizePath( args.list$rdp_path )
silva.path <- normalizePath( args.list$silva_path )

seqtab.nochim <- readRDS(args.list$seqtab_file_path)

## Asign SILVA and RDP taxonomies and merge with OTU table

# Assign SILVA taxonomy
taxa.silva <- dada2::assignTaxonomy(seqtab.nochim, silva.path, multithread = TRUE)

# Print first 6 rows of taxonomic assignment
unname(head(taxa.silva))


# OMIT APPENDING SPECIES FOR SILVA DUE TO MEMORY CONSTRAINTS
# Append species. Note that appending the argument 'allowMultiple=3' will return up to 3 different matched
# species, but if 4 or more are matched it returns NA.
#taxa.silva.species <- addSpecies(taxa.silva, silva.species.path)

# Merge with OTU table and save to file
otu.silva.tax.table <- merge( t(seqtab.nochim), taxa.silva, by = 'row.names' )
rownames( otu.silva.tax.table ) <- otu.silva.tax.table[,1]
otu.silva.tax.table <- otu.silva.tax.table[,-1]

otu.silva.tax.table_taxcombined <- cbind(otu.silva.tax.table)
colnum <- length(otu.silva.tax.table_taxcombined[1,])

otu.silva.tax.table_taxcombined <- otu.silva.tax.table_taxcombined[, -c((colnum-5): colnum)]

taxonomy <- vector()
taxonomy<- paste0("k__",as.character(otu.silva.tax.table$Kingdom),"; ",
                  "p__",as.character(otu.silva.tax.table$Phylum),"; ",
                  "c__",as.character(otu.silva.tax.table$Class),"; ",
                  "o__",as.character(otu.silva.tax.table$Order),"; ",
                  "f__",as.character(otu.silva.tax.table$Family),"; ",
                  "g__",as.character(otu.silva.tax.table$Genus),"; ")


otu.silva.tax.table_taxcombined <- cbind(otu.silva.tax.table_taxcombined,taxonomy)


write.table(otu.silva.tax.table_taxcombined, args.list$tax_closed_ref_silva , sep = "\t", eol = "\n", quote = F, col.names = NA)
write.table(otu.silva.tax.table, paste0(gsub(".tsv", "", args.list$tax_closed_ref_silva),"_taxcolumns.tsv") , sep = "\t", eol = "\n", quote = F, col.names = NA)

# Assign RDP taxonomy
taxa.rdp <- dada2::assignTaxonomy(seqtab.nochim,rdp.path, multithread = TRUE)


# OMIT APPENDING SPECIES FOR RDP DUE TO MEMORY CONSTRAINTS
# Append species. Note that appending the argument 'allowMultiple=3' will return up to 3 different matched
# species, but if 4 or more are matched it returns NA.
#taxa.rdp.species <- addSpecies(taxa.rdp, rdp.path)


# Merge with OTU table and save to file
otu.rdp.tax.table <- merge( t(seqtab.nochim), taxa.rdp, by = 'row.names' )
rownames( otu.rdp.tax.table ) <- otu.rdp.tax.table[,1]
otu.rdp.tax.table <- otu.rdp.tax.table[,-1]


otu.rdp.tax.table_taxcombined <- cbind(otu.rdp.tax.table)
colnum <- length(otu.rdp.tax.table_taxcombined[1,])

otu.rdp.tax.table_taxcombined <- otu.rdp.tax.table_taxcombined[, -c((colnum-5): colnum)]

taxonomy <- vector()
taxonomy<- paste0("k__",as.character(otu.rdp.tax.table$Kingdom),"; ",
                  "p__",as.character(otu.rdp.tax.table$Phylum),"; ",
                  "c__",as.character(otu.rdp.tax.table$Class),"; ",
                  "o__",as.character(otu.rdp.tax.table$Order),"; ",
                  "f__",as.character(otu.rdp.tax.table$Family),"; ",
                  "g__",as.character(otu.rdp.tax.table$Genus),"; ")


otu.rdp.tax.table_taxcombined <- cbind(otu.rdp.tax.table_taxcombined,taxonomy)


write.table(otu.rdp.tax.table_taxcombined, args.list$tax_closed_ref_rdp, sep = "\t", eol = "\n", quote = F, col.names = NA)
write.table(otu.rdp.tax.table, paste0(gsub(".tsv", "", args.list$tax_closed_ref_rdp, "_taxcolumns.tsv"), "_taxcolumns.tsv"), sep = "\t", eol = "\n", quote = F, col.names = NA)