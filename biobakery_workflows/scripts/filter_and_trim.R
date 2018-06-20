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
if(is.null(args.list$input_dir)) {
  stop("At least one argument must be supplied (input folder).\n", call.=FALSE)
}

# Print args list to STDOUT
for( i in names(args.list) ) {
  cat( i, "\t", args.list[[i]], "\n")
}

# print contents of folder
cat( grep( "*\\.fastq", list.files(args.list$input_dir), value=T ), sep = "\n" )

#print args.list$input_dir
#print args.list$output_dir 

# these variables are passed to the workflow
input.path <- normalizePath( args.list$input_dir )
#input.path <- args.list$input_dir

output.dir <- ifelse( is.null(args.list$output_dir), "output", args.list$output_dir )
pool.samples <- args.list$pool


# Variable "input.path" containing path to input fastq files directory 
# is inherited from wrapper script dada2_cli.r.

input.file.list <- grep( "*fastq", list.files( input.path ), value = T )
#input.path <- normalizePath("input/")

# List of input files

# Sort ensures forward/reverse reads are in same order
fnFs <- sort(grep( "_R1.*\\.fastq", list.files(input.path), value = T ) )
fnRs <- sort(grep( "_R2.*\\.fastq", list.files(input.path), value = T ) )

# Extract sample names, allowing variable filenames; e.g. *_R1[_001].fastq[.gz]

#sample.names <- gsub( "_R1.*\\.fastq(\\.gz)?", "", fnFs, perl = T)
#sample.namesR <- gsub( "_R2.*\\.fastq(\\.gz)?", "", fnRs, perl = T)

sample.names <- gsub( "_R1.*\\.fastq*", "", fnFs, perl = T)
sample.namesR <- gsub( "_R2.*\\.fastq*", "", fnRs, perl = T)

sample.ext <- file_ext(fnFs)
if(identical("gz",sample.ext[1])) sample.ext <- "fastq.gz"

if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")

# Specify the full path to the fnFs and fnRs

fnFs <- file.path(input.path, fnFs)
fnRs <- file.path(input.path, fnRs)
cat(fnFs)
cat(fnRs)

# Create filtered_input/ subdirectory for storing filtered fastq reads
filt_path <- file.path(output.dir, "filtered_input") 
ifelse(!dir.exists(filt_path), dir.create(filt_path, recursive = TRUE), FALSE)

readQC.folder <- file.path(output.dir, "Read_QC")
ifelse(!dir.exists(readQC.folder), dir.create(readQC.folder, recursive = TRUE), FALSE)

# Define filenames for filtered input files
#filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
#filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.", sample.ext))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.", sample.ext))


# Filter the forward and reverse reads:
# Note that:
# 1. Reads are both truncated and then filtered using the maxEE expected errors algorighm from UPARSE.
# 2. Reverse reads are truncated to shorter lengths than forward since they are much lower quality.
# 3. _Both_ reads must pass for the read pair to be output.
# 4. Output files are compressed by default.

rd.counts <- as.data.frame(
  filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,200),
                maxN=0, maxEE=c(1,2), truncQ=2, rm.phix=TRUE,
                compress=TRUE, multithread=TRUE) 
)
# Table of before/after read counts
rd.counts$ratio <- round( rd.counts$reads.out / rd.counts$reads.in, digits = 2 )
rd.counts

# Write rd.counts table to file in readQC.folder
saveRDS(rd.counts, paste0( readQC.folder, "/Read_counts_filt.rds" ))
write.table( rd.counts, paste0( readQC.folder, "/Read_counts_after_filtering.tsv" ), sep = "\t", quote = F, eol = "\n", col.names = NA )