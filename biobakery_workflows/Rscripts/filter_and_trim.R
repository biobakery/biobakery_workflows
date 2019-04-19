#!/usr/bin/env Rscript

# Load packages
library(dada2); packageVersion("dada2")
library(gridExtra)
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

# Print content of folder
cat( grep( "*\\.fastq", list.files(args.list$input_dir), value=T ), sep = "\n" )

# These variables are passed to the workflow
input.path <- normalizePath( args.list$input_dir )

output.dir <- ifelse( is.null(args.list$output_dir), "output", args.list$output_dir )

pair_id1 <- args.list$pair_id
pair_id2 <- sub("1","2",pair_id1)

# List of input files
# Sort ensures forward/reverse reads are in same order
fnFs <- sort(grep(paste0(pair_id1,".*\\.fastq"), list.files(input.path), value = T ) )
fnRs <- sort(grep(paste0(pair_id2,".*\\.fastq"), list.files(input.path), value = T ) )

# Extract sample files extension
sample.ext <- tools::file_ext(fnFs)

# Extract sample extension
if(identical("gz",sample.ext[1])){
  sample.ext <- "fastq.gz"
  }
# Extract sample names  allowing variable filenames
sample.names <- gsub( paste0(pair_id1,".*\\.",sample.ext), "", fnFs, perl = T)
sample.namesR <- gsub( paste0(pair_id2,".*\\.",sample.ext), "", fnRs, perl = T)

if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")

# Specify the full path to the fnFs and fnR
fnFs <- file.path(input.path, fnFs)
fnRs <- file.path(input.path, fnRs)

# Create filtered_input/ subdirectory for storing filtered fastq reads
filt_path <- file.path(output.dir, args.list$filtered_dir) 
ifelse(!dir.exists(filt_path), dir.create(filt_path, recursive = TRUE), FALSE)

# Generate plots and save to file
# Forward reads
fwd.qc.plots.list <- list()
for( i in 1 : length(fnFs)) {
  fwd.qc.plots.list[[i]] <- dada2::plotQualityProfile(fnFs[i])
  rm(i)
}
# Save to file
png(args.list$reads_plotF)
gridExtra::marrangeGrob( fwd.qc.plots.list, ncol=2, nrow=3, top = NULL )
dev.off()
rm(fwd.qc.plots.list)

# Reverse reads
rev.qc.plots.list <- list()
for( i in 1 : length(fnRs)) {
  rev.qc.plots.list[[i]] <- dada2::plotQualityProfile(fnRs[i])
  rm(i)
}
# Save to file
png(args.list$reads_plotR)
gridExtra::marrangeGrob( rev.qc.plots.list, ncol=2, nrow=3, top = NULL )
dev.off()
rm(rev.qc.plots.list)


# Define filenames for filtered input files
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.", sample.ext))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.", sample.ext))

# Filter the forward and reverse reads:
# Note that:
# 1. Reads are both truncated and then filtered using the maxEE expected errors algorighm from UPARSE.
# 2. Reverse reads are truncated to shorter lengths than forward since they are much lower quality.
# 3. _Both_ reads must pass for the read pair to be output.
# 4. Output files are compressed by default.
trunc_len_max2 <- strtoi(args.list$trunc_len_max)
if(trunc_len_max2 == 0){
  trunc_len_max1 <- 0
}else{
trunc_len_max1 <- trunc_len_max2 + 40
}
maxee1 <- strtoi(args.list$maxee)
maxee2 <- maxee1 * 2

rd.counts <- as.data.frame(
  dada2::filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(trunc_len_max1,trunc_len_max2),
                       minLen = 50, maxN=0, maxEE=c(maxee1,maxee2), truncQ=2, rm.phix=TRUE,
                compress=TRUE, multithread=as.numeric(args.list$threads)) 
)
# Table of before/after read counts
rd.counts$ratio <- round( rd.counts$reads.out / rd.counts$reads.in, digits = 2 )

# Write rd.counts table to file in output folder
saveRDS(rd.counts,  args.list$readcounts_rds_path )
write.table( rd.counts, args.list$readcounts_tsv_path, sep = "\t", quote = F, eol = "\n", col.names = NA )