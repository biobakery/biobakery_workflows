#!/usr/bin/env Rscript

# load packages
library(dada2); packageVersion("dada2")
#library(ggplot2)
library(msa)
#library(gridExtra)
library(phangorn)

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

# Get sequences
seqs <- getSequences(seqtab.nochim)
names(seqs) <- seqs # This propagates to the tip labels of the tree
# Multiple seqeuence alignment
mult <- msa(seqs, method="ClustalOmega", type="dna", order="input")
# Save msa to file; convert first to phangorn object
phang.align <- as.phyDat(mult, type="DNA", names=getSequences(seqtab.nochim))
write.phyDat(phang.align, format = 'fasta', file = paste0( output.path,"/msa.fasta") )

# Call FastTree (via 'system') to reconstruct phylogeny
system( paste( "FastTree -gtr -nt ", output.path, "/msa.fasta > ", output.path, "/FastTree.tre", sep = '' ) )


detach("package:phangorn", unload=TRUE)
detach("package:msa", unload=TRUE)