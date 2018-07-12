#!/usr/bin/env Rscript

# load packages
library(dada2); packageVersion("dada2")
library(tools)

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
cat( grep( "*\\.fastq*", list.files(args.list$output_dir), value=T ), sep = "\n" )

# these variables are passed to the workflow
output.path <- normalizePath(args.list$output_dir )
print(output.path)

#Filtered files folder path
filt_path <- file.path(output.path, "filtered_input") 

# Sort ensures forward/reverse reads are in same order
fnFs <- sort(grep( "_F_filt.*\\.fastq", list.files(filt_path), value = T ) )
fnRs <- sort(grep( "_R_filt.*\\.fastq", list.files(filt_path), value = T ) )


# Extract sample names, allowing variable filenames; e.g. *_R1[_001].fastq[.gz]
sample.names <- gsub( "_F_filt.*\\.fastq", "", fnFs, perl = T)
sample.namesR <- gsub( "_R_filt.*\\.fastq", "", fnRs, perl = T)

sample.ext <- tools::file_ext(fnFs)
if(identical("gz",sample.ext[1])){
  sample.ext <- "fastq.gz"
  # Extract sample names, allowing variable filenames
  sample.names <- gsub( "_F_filt.*\\.fastq.gz", "", fnFs, perl = T)
  sample.namesR <- gsub( "_R_filt.*\\.fastq.gz", "", fnRs, perl = T)
}


if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")

cwd <- getwd()

# Define filenames for filtered input files
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.", sample.ext))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.", sample.ext))

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
  derepF <- dada2::derepFastq(filtFs[[sam]])
  ddF <- dada2::dada(derepF, err=errF, multithread=TRUE)
  derepR <- dada2::derepFastq(filtRs[[sam]])
  ddR <- dada2::dada(derepR, err=errR, multithread=TRUE)
  merger <- dada2::mergePairs(ddF, derepF, ddR, derepR)
  mergers[[sam]] <- merger
}
rm(derepF); rm(derepR)
saveRDS(mergers, paste0(output.path, "/mergers.rds"))