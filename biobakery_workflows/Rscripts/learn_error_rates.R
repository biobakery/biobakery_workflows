#!/usr/bin/env Rscript

# Load packages
library(dada2); packageVersion("dada2")
library(ggplot2)

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

# Print content of folder
cat( grep( "*\\.fastq.gz", list.files(args.list$output_dir), value=T ), sep = "\n" )

# These variables are passed to the workflow
output.path <- normalizePath( args.list$output_dir )
print(output.path)

filt.path = file.path(output.path,args.list$filtered_dir)

filtFs <- file.path(filt.path,sort(grep( "*_F_filt.fastq*", list.files(filt.path), value = T ) ))
filtRs <- file.path(filt.path,sort(grep( "*_R_filt.fastq*", list.files(filt.path), value = T ) ))

set.seed(100)
# Filtered forward read error rates
errF <- dada2::learnErrors(filtFs, nread=1e6, multithread=as.numeric(args.list$threads))
# Filtered reverse read error rates
errR <- dada2::learnErrors(filtRs, nread=1e6, multithread=as.numeric(args.list$threads))


# Visualize the estimated error rates
ggplot2::ggsave(args.list$error_ratesF_png, dada2::plotErrors(errF, nominalQ=TRUE) , device = "png")
ggplot2::ggsave(args.list$error_ratesR_png, dada2::plotErrors(errR, nominalQ=TRUE) , device = "png")

# Save as rds files
saveRDS(errF, args.list$error_ratesF_path) 
saveRDS(errR, args.list$error_ratesR_path)