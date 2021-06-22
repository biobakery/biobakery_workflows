#!/usr/bin/env Rscript

# A utility script to compute beta diversity.
# Based on stats visualizations and notes from Kelsey Thompson.

library(vegan)
library(ggplot2)
library(optparse)
library(gridExtra)
library(dplyr)
library(tibble)
library(purrr)
library(readr)

# Add command line arguments #
options <- optparse::OptionParser(
    usage = paste("%prog [options]", " <data.tsv> ", " <metadata.tsv> ", " output_folder ")
)

options <- optparse::add_option(options,
    c("-a", "--min_abundance"),
    type = "double",
    dest = "min_abundance",
    default = 0.0,
    help = paste0("The minimum abundance for each feature",
    " [ Default: %default ]"
    )
)

options <- optparse::add_option(options,
    c("-m", "--min_prevalence"),
    type = "double",
    dest = "min_prevalence",
    default = 0.0,
    help = paste0("The minimum percent of samples for which",
    "a feature is detected at minimum abundance",
    " [ Default: %default ]"
    )
)

options <- optparse::add_option(options,
    c("-f", "--max_missing"),
    type = "double",
    dest = "max_missing",
    default = 20.0,
    help = paste0("The max percent of data for which",
    "a feature can have missing values",
    " [ Default: %default ]"
    )
)

option_not_valid_error <- function(message, valid_options) {
    print(paste(message, ": %s"), toString(valid_options))
    stop("Option not valid", call. = FALSE)
}

# get command line options and positional arguments
parsed_arguments = optparse::parse_args(options, positional_arguments = TRUE)
current_args <- parsed_arguments[["options"]]
positional_args <- parsed_arguments[["args"]]

# check three positional arguments are provided
if (length(positional_args) != 3) {
    optparse::print_help(options)
    stop(paste("Please provide the required","positional arguments","<data.tsv> <metadata.tsv> <output_folder>"))
}

# Input variables
## assume data file has samples as columns
## allow for pounds in header lines
data <- t(read.table(positional_args[1], header = TRUE, row.names = 1, sep="\t", comment.char = ""))
metadata <- data.frame(read.table(positional_args[2], header = TRUE, row.names = 1, sep="\t", comment.char = ""))

# create output folder if it does not exist
if (! dir.exists(positional_args[3]) ) {
  dir.create(positional_args[3])
}

# check for samples as columns or rows
samples_rows <- intersect(rownames(metadata),rownames(data))
if (length(samples_rows) < 1) {
    # allow for special chars in names
    original_metadata_rownames <- rownames(metadata)
    rownames(metadata) <- make.names(rownames(metadata))
    samples_rows <- intersect(rownames(metadata),rownames(data))
    if (length(samples_rows) < 1) {
        samples_col <- intersect(rownames(metadata),colnames(data))
        if (length(samples_col) < 1){
            rownames(metadata) <- original_metadata_rownames
            sample <- colnames(metadata)
            metadata <- type.convert(as.data.frame(t(metadata)))
        } else {
            sample <- colnames(data)
            data <- as.data.frame(t(data))
        }
    }
}

# Filter by abundance using zero as value for NAs
total_samples <- nrow(metadata)
min_samples <- total_samples * current_args$min_prevalence

data_zeros <- data
data_zeros[is.na(data_zeros)] <- 0
filtered_data <- data[,colSums(data_zeros > current_args$min_abundance) > min_samples, drop = FALSE]

# remove empty rows from data
filtered_data <- filtered_data[rowSums(filtered_data != 0, na.rm=TRUE) > 0, , drop = FALSE]

if ((ncol(filtered_data) < 1) || (nrow(filtered_data) < 1)) {
  stop("No data remain in the data after filtering for min abundance and prevalence")
}

# Filter the metadata to remove any samples without data or metadata with more than max values missing
metadata_zeros <- metadata
metadata_zeros[metadata_zeros == "UNK"] <- NA

max_not_missing_values <- (( 100 - current_args$max_missing ) / 100 ) * nrow(metadata)
filtered_metadata <- na.omit(metadata_zeros[ ,colSums(metadata_zeros !='NA', na.rm=TRUE) > max_not_missing_values, drop = FALSE])

if ((ncol(filtered_metadata) < 1) || (nrow(filtered_metadata) < 1)) {
  stop("No metadata remain after filtering for max missing")
}

# filter data and metadata to only include the same samples in the same order
sorted_samples <- sort(intersect(row.names(filtered_metadata), row.names(filtered_data)))
filtered_metadata <- filtered_metadata[sorted_samples, , drop = FALSE]
filtered_data <- filtered_data[sorted_samples, , drop = FALSE]

# remove subject from metadata if present
filtered_metadata <- filtered_metadata[ , !(names(filtered_metadata) %in% c("subject")), drop = FALSE]

# creates alpha diversity results appended to metadata
index <- "invsimpson"
alpha_div <- data.frame(diversity(filtered_data, index = index))
names(alpha_div) <- index
alpha_with_metadata <- inner_join(rownames_to_column(alpha_div), rownames_to_column(filtered_metadata), by = "rowname") %>% column_to_rownames()

# write the data to the file
alpha_with_metadata %>% as_tibble() %>% write_tsv(file.path(positional_args[3], "alpha_diversity_with_metadata.txt"))

# Function credit to Tom (with some modifications for colors and dependencies)
alpha_plot <- function(alpha_meta, category)
{
  index <- names(alpha_meta)[1]
  
  # drop samples with NA for given category
  meta <- alpha_meta[complete.cases(alpha_meta[c(index, category)]), , drop=FALSE]
  
  if (is.character(meta[[category]]) | is.factor(meta[[category]]) ) {
    filename <- file.path(positional_args[3], paste0(category,"_boxplot.png")) 
    print(filename)
    p <- ggplot(data = meta, aes(x = !!sym(category), y = invsimpson, fill = !!sym(category))) + geom_boxplot() + 
      xlab(category) + ylab(paste(index, "Diversity")) + 
      theme_bw(base_size = 14) + theme(axis.text.x = element_text(hjust=0.5, vjust=0.5))
      #scale_fill_manual(values = cat_colors[[category]])
    png(filename, res = 150, height = 800, width = 1100)
    plot(p)
    dev.off()
  } else if (is.numeric(meta[[category]])) {
    filename <- file.path(positional_args[3], paste0(category,"_scatterplot.png"))
    print(filename)
    p <- ggplot(data = meta, aes(x = !!sym(category), invsimpson)) + 
      geom_point(fill = "darkolivegreen4", color = "darkolivegreen4", alpha = 0.5, shape = 21, size = 5, stroke = 0.05) +
      scale_x_continuous(limits = c(min(meta[[category]]), max(meta[[category]]))) + 
      scale_y_continuous(limits = c(min(meta[[index]]), max(meta[[index]]))) + 
      stat_smooth(method = "lm", color = "blue", na.rm = T) + guides(alpha = "none") + 
      labs("") + xlab(category) + ylab(paste(index, "Diversity")) + theme_bw(base_size = 14)
    png(filename, res = 150, height = 800, width = 1100)
    plot(p)
    dev.off()
  } else { stop(paste("Category of", typeof(metadata[[category]]), "not recognized")) }
}

walk(colnames(filtered_metadata), function(x) alpha_plot(alpha_with_metadata, x))
