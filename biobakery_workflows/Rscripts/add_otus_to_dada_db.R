#!/usr/bin/env Rscript
# This script adds otuids as taxonomy level to dada2 green genes train database
# Uses green genes fasta database to add otuids
# Assumes sequences are in matching order


## Collect arguments
args <- commandArgs(TRUE)

## Parse arguments (we expect the form --arg=value)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
args.df <- as.data.frame(do.call("rbind", parseArgs(args)))

args.list <- as.list(as.character(args.df$V2))
names(args.list) <- args.df$V1

## Args 
if(is.null(args.list$refdb_dada2_path)) {
  stop("Path to dada2 reference db must be supplied (refdb_dada2_path).\n", call.=FALSE)
}
if(is.null(args.list$refdb_usearch_path)) {
  stop("Path to usearch fasta reference db must be supplied (refdb_usearch_path).\n", call.=FALSE)
}
if(is.null(args.list$refdb_usearch_path)) {
  stop("Path to new reformated db must be supplied (refdb_reformat_path).\n", call.=FALSE)
}


# Print args list to STDOUT
for( i in names(args.list) ) {
  cat( i, "\t", args.list[[i]], "\n")
}

con_usearch = file(args.list$refdb_usearch_path, "r")
con_dada2 = file(args.list$refdb_dada2_path, "r")
conw = file(args.list$refdb_reformat_path, "w")
i <- 0
  while ( TRUE ) {
    line_usearch = readLines(con_usearch, n = 1)
    line_dada2 = readLines(con_dada2, n = 1)
    if ( length(line_usearch) == 0 | length(line_dada2) == 0 ) {
      break
    }
    line_reformat <- gsub(";$",paste0(gsub(">",";otuid__",line_usearch),";"),line_dada2)
    writeLines(line_reformat, conw, sep="\n")
    
    line1 = readLines(con_usearch, n = 1)
    line2 = readLines(con_dada2, n = 1)
    if (!identical(line1, line2)){print("Not matching sequences")}
    writeLines(line2, conw, sep="\n")
    i <- i+1
    print(i)
  }
  close(con_usearch)
  close(con_dada2)
  close(conw)