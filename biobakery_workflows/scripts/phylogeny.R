#!/usr/bin/env Rscript

# load packages
library(dada2); packageVersion("dada2")
library(msa)
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
output.path <- normalizePath( args.list$output_dir )


seqtab.nochim <- readRDS(args.list$seqtab_file_path)

# Get sequences
seqs <- dada2::getSequences(seqtab.nochim)
names(seqs) <- seqs # This propagates to the tip labels of the tree
# Multiple seqeuence alignment
mult <- msa::msa(seqs, method="ClustalOmega", type="dna", order="input")
# Save msa to file; convert first to phangorn object
phang.align <- phangorn::as.phyDat(mult, type="DNA", names=dada2::getSequences(seqtab.nochim))
write.phyDat(phang.align, format = 'fasta', file = args.list$msa_fasta_path)

detach("package:phangorn", unload=TRUE)
detach("package:msa", unload=TRUE)