#!/usr/bin/env Rscript

# Load packages
library(dada2); packageVersion("dada2")
library(tools)
library(seqinr)

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
output.path <- normalizePath(args.list$output_dir )

# Filtered files folder path
filt_path <- file.path(output.path, args.list$filtered_dir) 

# Sort ensures forward/reverse reads are in same order
fnFs <- sort(grep( "_F_filt.*\\.fastq", list.files(filt_path), value = T ) )
fnRs <- sort(grep( "_R_filt.*\\.fastq", list.files(filt_path), value = T ) )


# Extract sample extensions
sample.ext <- tools::file_ext(fnFs)
if(identical("gz",sample.ext[1])){
  sample.ext <- "fastq.gz"
}
# Extract sample names, allowing variable filenames
sample.names <- gsub( paste0("_F_filt.*\\.", sample.ext), "", fnFs, perl = T)
sample.namesR <- gsub( paste0("_R_filt.*\\.", sample.ext), "", fnRs, perl = T)

# Read merged reads from file
mergers <- readRDS(args.list$merged_file_path)

# Construct sequence table (ASV)
seqtab <- dada2::makeSequenceTable(mergers)
dim(seqtab)
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
# The sequence table is a matrix with rows corresponding to (and named by) the samples, and columns corresponding to (and named by) the sequence variants.

# Remove chimeric sequences:
seqtab.nochim <- dada2::removeBimeraDenovo(seqtab, method="consensus", multithread=as.numeric(args.list$threads), verbose=TRUE)
dim(seqtab.nochim)

# Ratio of chimeric sequence reads
1 - sum(seqtab.nochim)/sum(seqtab)

# Write ASV  table to tsv file
write.table( t(seqtab.nochim), paste0( output.path, "/", args.list$asv_tsv), sep = "\t", eol = "\n", quote = F, col.names = NA )
# Save ASV table to rds file for further processing
saveRDS(seqtab.nochim, args.list$seqtab_file_path)


# Get sequences, assign ids
seqs <- dada2::getSequences(seqtab.nochim)
seqids <- c(1:length(seqs))
seqids <- paste0("ASV",seqids)
names(seqs) <- seqids 

# Write sequences to fasta file for further processing
seqinr::write.fasta(sequences = as.list(seqs), names = names(seqs), file.out=args.list$seqs_fasta_path, open = "w",  nbchar = 60, as.string = FALSE)

# Save read counts on each step to rds file
rd.counts <- readRDS(paste0(output.path, "/", args.list$readcounts_rds ))

# Remove rows with 0 reads after filtering and trimming
rdf.counts <- rd.counts[which(rd.counts$reads.out != 0),]

getN <- function(x) sum(getUniques(x))
track <- cbind(rdf.counts, sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "ratio", "merged", "tabled", "nonchim")
rownames(track) <- sample.names

# Print table
track
# Save read count to tsv file
write.table( track,  args.list$read_counts_steps_path, sep = "\t", quote = F, eol = "\n", col.names = NA )