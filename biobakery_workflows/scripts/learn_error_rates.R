#!/usr/bin/env Rscript

# load packages
library(dada2); packageVersion("dada2")
library(ggplot2)
#library(msa)
#library(gridExtra)
#library(phangorn)


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

# print contents of folder
cat( grep( "*\\.fastq.gz", list.files(args.list$output_dir), value=T ), sep = "\n" )

# these variables are passed to the workflow
output.path <- normalizePath( args.list$output_dir )
print(output.path)

filt.path = file.path(output.path,"filtered_input")

filtFs <- file.path(filt.path,sort(grep( "*_F_filt.fastq.gz", list.files(filt.path), value = T ) ))
filtRs <- file.path(filt.path,sort(grep( "*_R_filt.fastq.gz", list.files(filt.path), value = T ) ))

readQC.folder <- file.path(output.path, "Read_QC")
ifelse(!dir.exists(readQC.folder), dir.create(readQC.folder, recursive = TRUE), FALSE)


set.seed(100)
# Filtered forward reads
errF <- learnErrors(filtFs, nread=1e3, multithread=TRUE)
# Filtered reverse reads
errR <- learnErrors(filtRs, nread=1e3, multithread=TRUE)
# Visualize the estimated error rates
# Forward
#plotErrors(errF, nominalQ=TRUE)
# Save to file 
ggsave(paste0(readQC.folder,"/Error_rates_per_sample_FWD.pdf"), plotErrors(errF, nominalQ=TRUE) , device = "pdf")
# Reverse
#plotErrors(errR, nominalQ=TRUE)
# Save to file 
ggsave(paste0(readQC.folder,"/Error_rates_per_sample_REV.pdf"), plotErrors(errR, nominalQ=TRUE) , device = "pdf")

#save as rds files
saveRDS(errF, paste0(output.path, "/error_rates_F.rds")) 
saveRDS(errR, paste0(output.path, "/error_rates_R.rds"))