#!/usr/bin/env Rscript

# This script generates an ordination plot from the distance matrix calculated from StrainPhlAn multiple sequence alignment file

# load the optparse library
library(optparse)

# Create command line argument parser
help_description <- "

This script generates an ordination plot from StrainPhlAn output files.

The following positional arguments are required:

clade.distmat.txt (Input file): This is the distmat output file (created by providing the StrainPhlAn MSA file)
metadata.txt (Input file): This is the metadata file
ordination.png (Output file): This is the ordination plot image written"

args <- OptionParser(usage = "%prog clade.distmat.txt metadata.txt ordination.png",
                      add_help_option = TRUE, prog='strainphlan_ordination.R',
                      description=help_description )
args_list <- parse_args(args, positional_arguments=TRUE)

# load libraries after optparse to not load them on help
library(ggplot2)
library(vegan)

# read in the file, skipping the first 8 rows and filling in empty columns, using the tab as sep, and stripping extra white space
data <- read.table( args_list$args[1], skip = 8, fill = TRUE, sep="\t", strip.white = T)

list = data$V33
list = gsub(" .*", "", list)
# remove the first column of the data as it is blank
data[1] <- NULL

# get the header as the last column of the data as a character vector
header <- lapply(data[,ncol(data)], as.character)

# remove the last column from the data as it has been stored as a header
data[ncol(data)] <- NULL

# remove the current last column from the data as it is blank
data[ncol(data)] <- NULL

# add the sample names to the columns and rows of the data matrix
rownames(data) <- list
colnames(data) <- list

# make symmetric, add lower triangle to upper triangle
data[lower.tri(data)] <- t(data)[lower.tri(data)]


# ordinate on the distance matrix
e.sir.pcoa <- cmdscale( data, eig = T )

# variance explained 
variance <- head(eigenvals(e.sir.pcoa)/sum(eigenvals(e.sir.pcoa)))
x_variance <- as.integer(variance[1]*100)
y_variance <- as.integer(variance[2]*100)

# get scores for plotting
e.sir.scores <- as.data.frame( e.sir.pcoa$points )

# read in metadata file
e.sir.meta <- read.delim( args_list$args[2], header = T, sep = "\t", row.names = 1 )
# append to e.sir.scores
e.sir.scores.meta <- merge( e.sir.scores, e.sir.meta, by = 'row.names' )
# set first column as rownames and remove it
rownames( e.sir.scores.meta ) <- e.sir.scores.meta[,1]
e.sir.scores.meta[,1] <- NULL
# change colnames
colnames(e.sir.scores.meta) <- c( "PCo1", "PCo2", "Country" )

# plot ordination
png( args_list$args[3], width = 6, height = 5, res = 300, units = "in" )

ggplot( e.sir.scores.meta, aes(PCo1, PCo2, color=Country) ) + 
  geom_point(size = 4, alpha = 0.75) + theme_classic() + 
  theme(axis.line.x = element_line(colour = 'black', size=0.75, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.75, linetype='solid'),
        axis.ticks = element_blank(), axis.text = element_blank()) + 
  xlab(paste("PCo1 (",x_variance,"% variance explained)")) + ylab(paste("PCo2 (",y_variance,"% variance explained)"))

temp <- dev.off()
