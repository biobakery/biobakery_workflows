#!/usr/bin/env Rscript

# A utility script to compute beta diversity.
# Based on stats visualizations and notes from Kelsey Thompson.

library(vegan)
library(ggplot2)
library(optparse)
library(gridExtra)

# Add command line arguments #
options <- optparse::OptionParser(
    usage = paste("%prog [options]", " <data.tsv> ", " <metadata.tsv> ", " <output_file.png> ")
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

options <- optparse::add_option(options,
    c("-c", "--covariate_equation"),
    type = "character",
    dest = "covariate_equation",
    default = "",
    help = paste0("The equation for multi-variate studies",
    " [ Default: %default ]"
    )
)

options <- optparse::add_option(options,
    c("-d", "--data_type"),
    type = "character",
    dest = "data_type",
    default = "relab",
    help = paste0("The type of data provided",
    "(examples are 'relab', 'present absent', and 'rna dna norm')",
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
    stop(paste("Please provide the required","positional arguments","<data.tsv> <metadata.tsv> <output_file.png>"))
}

# Input variables
## assume data file has samples as columns
## allow for pounds in header lines
data <- t(read.table(positional_args[1], header = TRUE, row.names = 1, sep="\t", comment.char = ""))
metadata <- data.frame(read.table(positional_args[2], header = TRUE, row.names = 1, sep="\t", comment.char = ""))

# check for samples as columns or rows
samples_rows <- intersect(rownames(metadata),rownames(data))
if (length(samples_rows) < 1) {
    sample <- colnames(metadata)
    metadata <- as.data.frame(t(metadata))
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
filtered_metadata <- na.omit(metadata_zeros[ ,colSums(metadata_zeros != 0, na.rm=TRUE) > max_not_missing_values, drop = FALSE])

if ((ncol(filtered_metadata) < 1) || (nrow(filtered_metadata) < 1)) {
  stop("No metadata remain after filtering for max missing")
}

# filter data and metadata to only include the same samples in the same order
sorted_samples <- sort(intersect(row.names(filtered_metadata), row.names(filtered_data)))
filtered_metadata <- filtered_metadata[sorted_samples, , drop = FALSE]
filtered_data <- filtered_data[sorted_samples, , drop = FALSE]

# remove subject from metadata if present
filtered_metadata <- filtered_metadata[ , !(names(filtered_metadata) %in% c("subject"))]

# compute univariate beta diversity
if (current_args$data_type == "relab") {
  method = "bray"
}
if (current_args$data_type == "present absent") {
  method="jaccard"
}
if (current_args$data_type == "dna rna norm") {
  method="euclidean"
}

# check if valid data type was set
if (! exists("method")) {
  stop(paste("Please provide a data type from the valid set in the help message"))
}

bray = vegdist(filtered_data, method, na.remove = TRUE)

if (current_args$covariate_equation != "") {
  results <- adonis(as.formula(paste("bray ~ ", current_args$covariate_equation)), data = filtered_metadata)
  png(positional_args[3], res=150, height=800, width=1300)
  grid.table(as.data.frame(results$aov.tab))
  dev.off()

} else {
  adonis_pval = vector()
  adonis_rsq = vector()
  for (col in names(filtered_metadata)){
    adonis.univ = adonis(as.formula(paste("bray ~", col)), data = filtered_metadata)
    adonis_pval[col] = adonis.univ$aov.tab[1,]$`Pr(>F)`
    adonis_rsq[col] = adonis.univ$aov.tab[1,]$R2
  }

  univar_tax = rbind(adonis_pval, adonis_rsq)
  univar_tax = as.data.frame(t(univar_tax))
  names(univar_tax) = c("P-Value", "R2")
  univar_tax$p_adj = p.adjust(univar_tax$`P-Value`, "fdr")

  univar_tax$stars = cut(univar_tax$`P-Value`, c(0, 0.001, 0.01, 0.05, 0.1, 1), labels = c("***", "**", "*", "`", ""))

  univar_tax$`P-Value` = round(univar_tax$`P-Value`, 3)
  univar_tax$R2 = univar_tax$R2 *100

  dodge = position_dodge(width = 0.8)

  plot <- ggplot(data = univar_tax, aes(reorder(row.names(univar_tax), univar_tax$R2), y = R2, label = univar_tax$`P-Value`)) + geom_bar(stat = "identity", position = "identity", fill = "#800000") + geom_text(position = dodge, vjust = 0.5, hjust = -0.1, size = 3) + theme_bw(base_size = 12) + ylab("Univarate R-squared") + coord_flip() + ylim(0, as.integer(max(univar_tax$R2))+1) + xlab("") + labs(fill = "")

  png(positional_args[3], res = 150, height = 800, width = 1100)
  print(plot)
  dev.off()
}

