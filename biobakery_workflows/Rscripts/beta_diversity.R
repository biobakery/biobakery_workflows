#!/usr/bin/env Rscript

# A utility script to compute beta diversity.
# Based on stats visualizations and notes from Kelsey Thompson.

library(vegan)
library(ggplot2)
library(optparse)
library(gridExtra)
library(tibble)
library(dplyr)
library(readr)

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
    "(examples are 'relab', 'present_absent', and 'rna_dna_norm')",
    " [ Default: %default ]"
    )
)

options <- optparse::add_option(options,
    c("-o", "--adonis_method"),
    type = "character",
    dest = "adonis_method",
    default = "",
    help = paste0("The method to use for the adnois",
    " [ Default: based on data_type (bray for relab) ]"
    )
)

options <- optparse::add_option(options,
    c("-p", "--pairwise"),
    type = "logical",
    dest = "pairwise",
    default = FALSE,
    help = paste0("Run pairwise comparison on all categorical features ",
    " [ Default: FALSE ]"
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

# compute univariate beta diversity
if (current_args$data_type == "relab") {
  method = "bray"
}
if (current_args$data_type == "present_absent") {
  method="jaccard"
}
if (current_args$data_type == "dna_rna_norm") {
  method="euclidean"
}

# use the user provided method, if set
if (current_args$adonis_method != "") {
  method=current_args$adonis_method
}

# check if valid data type was set
if (! exists("method")) {
  stop(paste("Please provide a data type from the valid set in the help message"))
}

bray = vegdist(filtered_data, method, na.remove = TRUE)

if (current_args$pairwise) {
  used_col <- character()
  tables <- list()
  i <- 1
  png(positional_args[3], res=150, height=800, width=1100)
  theme <- ttheme_default(base_size = 8, padding = unit(c(4, 4), "mm"))
  for (col in names(filtered_metadata)){
      if (is.numeric(filtered_metadata[[col]])) {
        for (col2 in names(filtered_metadata)) {
            if (is.numeric(filtered_metadata[[col2]]) && !(col2 %in% used_col) && (col != col2)) {
                results <- adonis2(as.formula(paste("bray ~", col, " + ", col2)), data = filtered_metadata, by="margin")
                grid <- tableGrob(as.data.frame(results), theme=theme)
                tables[[i]] <- grid
                i <- i + 1
            }
        }
       used_col <- append(used_col, col)
    }
  }
  grid.arrange(grobs=tables, nrow=length(tables), ncol=1, main = "")
  dev.off()
} else if (current_args$covariate_equation != "") {

  # check for filtered covariates
  if (length(names(metadata)) != length(names(filtered_metadata)))
  {
    diff <- setdiff(names(metadata),names(filtered_metadata))
    if (! length(diff) == 1 && diff[1] == "subject") {
      stop("Metadata covariates have been filtered from the data set that are included in the covariate equation: ", paste(diff, collapse=" , "))
    }
    # remove any metadata variables from the equation that might have been filtered
    for (col in names(metadata)) {
      if (! (col %in% names(filtered_metadata))) {
        check_pattern <- paste(col,"+ ")
        if (grepl(check_pattern,current_args$covariate_equation,fixed=TRUE)) {
          current_args$covariate_equation <- sub(check_pattern,"",current_args$covariate_equation,fixed=TRUE)
        }
        check_pattern <- paste0(col,";")
        check_equation <- paste0(current_args$covariate_equation,";")
        if (grepl(check_pattern,check_equation,fixed=TRUE)) {
          current_args$covariate_equation <- sub(check_pattern,"",check_equation,fixed=TRUE)
        }
      }
    }
  }

  results <- adonis2(as.formula(paste("bray ~ ", current_args$covariate_equation)), data = filtered_metadata, by="margin")
  png(positional_args[3], res=150, height=800, width=1100)

  if (length(names(metadata)) > 20) {
    theme <- ttheme_default(base_size = 4, padding = unit(c(2, 2), "mm"))
  } else {
    theme <- ttheme_default(base_size = 8, padding = unit(c(4, 4), "mm"))
  }
  grid.table(as.data.frame(results), theme=theme)
  dev.off()

  # write the equation to a file
  equation_file <- sub(".png","_equation.txt",positional_args[3],fixed=TRUE)
  write(current_args$covariate_equation,equation_file)

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

  plot <- ggplot(data = univar_tax, aes(reorder(row.names(univar_tax), univar_tax$R2), y = R2, label = univar_tax$`P-Value`)) + geom_bar(stat = "identity", position = "identity", fill = "#800000") + geom_text(position = dodge, vjust = 0.5, hjust = -0.1, size = 3) + theme_bw(base_size = 12) + ylab("Univariable R-squared") + coord_flip() + ylim(0, as.integer(max(univar_tax$R2)+1.0)*1.05) + xlab("") + labs(fill = "")

  png(positional_args[3], res = 150, height = 800, width = 1100)
  print(plot)
  dev.off()
  outfile <- gsub(".png",".txt",positional_args[3])
  univar_tax %>% as_tibble() %>% write_tsv(outfile)
}

