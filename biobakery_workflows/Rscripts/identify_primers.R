#!/usr/bin/env Rscript

# Load packages
library(dada2); packageVersion("dada2")
library(ShortRead); packageVersion("ShortRead")
library(Biostrings); packageVersion("Biostrings")
library(tools)

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  # The Biostrings works w/ DNAString objects rather than character vectors
  dna <- Biostrings::DNAString(primer) 
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  # Convert back to character vector
  return(sapply(orients, toString)) 
}

## Collect arguments
args <- commandArgs(TRUE)

## Parse arguments (we expect the form --arg=value)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
args.df <- as.data.frame(do.call("rbind", parseArgs(args)))

args.list <- as.list(as.character(args.df$V2))
names(args.list) <- args.df$V1

# Print args list to STDOUT
for( i in names(args.list) ) {
  cat( i, "\t", args.list[[i]], "\n")
}
# Print content of folder
cat( grep( "*\\.fastq", list.files(args.list$input_dir), value=T ), sep = "\n" )

# These variables are passed to the workflow
input.path <- normalizePath( args.list$input_dir )

pair_id1 <- args.list$pair_id
pair_id2 <- sub("1","2", pair_id1)

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

# Identifay primers and orientation
FWD.orients <- allOrients(args.list$fwd_primer)
REV.orients <- allOrients(args.list$rev_primer)

FWD.orients
REV.orients

# Filter all Ns and put N-filtered files in filtN subdirectory
path.filtN <- args.list$filtn_dir
fnFs.filtN <- file.path(path.filtN, basename(fnFs))
fnRs.filtN <- file.path(path.filtN, basename(fnRs))
if(!dir.exists(path.filtN)) dir.create(path.filtN)

dada2::filterAndTrim(file.path(input.path,fnFs), fnFs.filtN, file.path(input.path,fnRs), fnRs.filtN, maxN = 0, multithread = TRUE)

primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- Biostrings::vcountPattern(primer, ShortRead::sread(ShortRead::readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

primers=rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))

print("Primers identified")
primers

if(primers[1,1] > primers[1,4]){
  FWD <- args.list$fwd_primer
  FWD.RC <- dada2:::rc(FWD)
}else{
  FWD.RC <- args.list$fwd_primer
  FWD <- dada2:::rc(FWD.RC)
}
if(primers[4,1] > primers[4,4]){
  REV <- args.list$rev_primer
  REV.RC <- dada2:::rc(REV)
}else{
  REV.RC <- args.list$rev_primer
  REV <- dada2:::rc(REV.RC)
}
print("Forward primer")
FWD
print("Reverse primer")
REV

if(!dir.exists(args.list$primers_dir)) dir.create(args.list$primers_dir)
fwd_primer_file<-file(args.list$fwd_primer_file)
writeLines(c(FWD,FWD.RC), fwd_primer_file)
close(fwd_primer_file)
rev_primer_file<-file(args.list$rev_primer_file)
writeLines(c(REV,REV.RC), rev_primer_file)
close(rev_primer_file)


