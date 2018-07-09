#!/usr/bin/env Rscript

# load packages
library(dada2); packageVersion("dada2")
#library(ggplot2)
#library(msa)
#library(gridExtra)
#library(phangorn)
library(tools)

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


# these variables are passed to the workflow
input.path <- normalizePath( args.list$input_dir )
output.path <- normalizePath( args.list$output_dir )
#pool.samples <- args.list$pool

#Filtered files folder path
filt_path <- file.path(output.path, "filtered_input") 

# Sort ensures forward/reverse reads are in same order
fnFs <- sort(grep( "_F_filt.*\\.fastq", list.files(filt_path), value = T ) )
fnRs <- sort(grep( "_R_filt.*\\.fastq", list.files(filt_path), value = T ) )


# Extract sample names, allowing variable filenames; e.g. *_R1[_001].fastq[.gz]
sample.names <- gsub( "_F_filt.*\\.fastq", "", fnFs, perl = T)
sample.namesR <- gsub( "_R_filt.*\\.fastq", "", fnRs, perl = T)

sample.ext <- file_ext(fnFs)
if(identical("gz",sample.ext[1])){
  sample.ext <- "fastq.gz"
  # Extract sample names, allowing variable filenames
  sample.names <- gsub( "_F_filt.*\\.fastq.gz", "", fnFs, perl = T)
  sample.namesR <- gsub( "_R_filt.*\\.fastq.gz", "", fnRs, perl = T)
}


mergers <- readRDS(file.path(output.path,"mergers.rds"))

seqtab <- makeSequenceTable(mergers)
dim(seqtab)
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
# The sequence table is a matrix with rows corresponding to (and named by) the samples, and columns corresponding to (and named by) the sequence variants. # Remove chimeric sequences:
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

dim(seqtab.nochim)
# ratio of chimeric sequence reads
1 - sum(seqtab.nochim)/sum(seqtab)

# write sequence variants count table to file
write.table( t(seqtab.nochim), paste0( output.path, "/all_samples_SV-counts.tsv"), sep = "\t", eol = "\n", quote = F, col.names = NA )
# write OTU table to file
saveRDS(seqtab.nochim, paste0(output.path, "/seqtab_final.rds"))

rd.counts <- readRDS(paste0(output.path, "/Read_counts_filt.rds" ))
# remove rows with 0 reads after filtering and trimming
#rdf.counts <- rd.counts[row(rd.counts)[rd.counts$reads.out != 0],]
rdf.counts <- rd.counts

getN <- function(x) sum(getUniques(x))
track <- cbind(rdf.counts, sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "ratio", "merged", "tabled", "nonchim")
rownames(track) <- sample.names
# print table
track
# save to file
write.table( track, paste0(output.path, "/Read_counts_at_each_step.tsv" ), sep = "\t", quote = F, eol = "\n", col.names = NA )