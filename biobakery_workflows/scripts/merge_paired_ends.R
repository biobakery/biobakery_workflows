#!/usr/bin/env Rscript

# load packages
library(dada2); packageVersion("dada2")
library(ggplot2)
library(msa)
library(gridExtra)
library(phangorn)

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

# print contents of folder
cat( grep( "*\\.fastq.gz", list.files(args.list$output_dir), value=T ), sep = "\n" )

# these variables are passed to the workflow
output.path <- normalizePath( args.list$output_dir )
print(output.path)

input.path <- normalizePath( args.list$input_dir )
print(input.path)

# Variable "input.path" containing path to input fastq files directory 
# is inherited from wrapper script dada2_cli.r.

input.file.list <- grep( "*fastq", list.files( input.path ), value = T )
#input.path <- normalizePath("input/")

# List of input files

# Sort ensures forward/reverse reads are in same order
fnFs <- sort(grep( "_R1.*\\.fastq", list.files(input.path), value = T ) )
fnRs <- sort(grep( "_R2.*\\.fastq", list.files(input.path), value = T ) )

# Extract sample names, allowing variable filenames; e.g. *_R1[_001].fastq[.gz]
sample.names <- gsub( "_R1.*\\.fastq(\\.gz)?", "", fnFs, perl = T)
sample.namesR <- gsub( "_R2.*\\.fastq(\\.gz)?", "", fnRs, perl = T)
if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")


# Specify the full path to the fnFs and fnRs

fnFs <- file.path(input.path, fnFs)
fnRs <- file.path(input.path, fnRs)
print(sample.names)
cwd <- getwd()
# Create filtered_input/ subdirectory for storing filtered fastq reads
filt_path <- file.path(output.path, "filtered_input") 
ifelse(!dir.exists(filt_path), dir.create(filt_path, recursive = TRUE), FALSE)

# Define filenames for filtered input files
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))


errF <- readRDS(file.path(output.path,"error_rates_F.rds"))
errR <- readRDS(file.path(output.path,"error_rates_R.rds"))


# Sample inference of dereplicated reads, and merger of paired-end reads
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names
names(filtFs) <- sample.names
names(filtRs) <- sample.names
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  print(filtFs[[sam]])
  derepF <- derepFastq(filtFs[[sam]])
  ddF <- dada(derepF, err=errF, multithread=TRUE)
  derepR <- derepFastq(filtRs[[sam]])
  ddR <- dada(derepR, err=errR, multithread=TRUE)
  merger <- mergePairs(ddF, derepF, ddR, derepR)
  mergers[[sam]] <- merger
}
rm(derepF); rm(derepR)
saveRDS(mergers, paste0( "/Users/anamailyan/dada_output", "/mergers.rds"))