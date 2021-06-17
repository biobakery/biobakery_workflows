#!/usr/bin/env Rscript

# A utility script using the functions from the hmp2_analysis repository.

library(vegan)
library(permute)
library(reshape2)
library(RColorBrewer)
library(cowplot)
library(ggplot2)
library(plyr)
library(optparse)
library(tibble)
library(dplyr)
library(readr)

# Stand-alone utility functions (from the hmp2_analysis repository file overview/src/omnibus_tests.r)

#' PERMANOVA with repeat measure-aware permutations. Block sizes are allowed to
#' differ.
#'
#' @param D An N-by-N distance matrix (must be a \code{dist} object).
#' @param permute_within Data frame with N rows containing metadata to test per
#' sample
#' @param blocks A length-N vector containing the block structure.
#' @param block_data Data frame with per-block metadata. If \code{blocks} is
#' numeric, \code{block_data}'s rows must match those indices. If \code{blocks}
#' is a \code{factor}, then the row ordering must match the factor levels. If
#' \code{blocks} is a character vector, \code{block_data} must have rownames
#' matching the contents of \code{blocks}.
#' @param permutations Number of permutations to test
#' @param metadata_order Order of the metadata in the model. If not given, this
#' is assumed to be within-block metadata first, followed by block metadata, in
#' the order given in \code{permute_within} and \code{block_data}.
#' @return Same structure as \code{adonis}, with p-values recalculated based on
#' permutations that are aware of the block structure of the data.
#' Metadata in \code{permute_within} are permuted within blocks, whereas
#' metadata in \code{block_data} are first permuted across blocks, and then
#' assigned to samples according to the block structure.
PERMANOVA_repeat_measures <- function(
    D, permute_within, blocks = NULL, block_data, permutations=999,
    metadata_order = c(names(permute_within), names(block_data)),
    na.rm=F) {

    # Make sure D is a dist object
    if (class(D) != "dist") {
        stop("D must be a dist object")
    }

    # Default to free permutations if blocks is not given
    if (!missing(block_data) && is.null(blocks)) {
        stop("blocks must be given if block_data is present")
    } else if (is.null(blocks)) {
        blocks <- rep(1, nrow(permute_within))
        block_data <- as.data.frame(matrix(0, nrow=1, ncol=0))
    } else if (length(unique(blocks)) == 1) {
        warning("blocks only contains one unique value")
    }

    # Ensure no metadata overlap between permute_within and block_data
    if (length(intersect(names(permute_within), names(block_data))) > 0) {
        stop("metadata is repeated across permute_within and block_data")
    }

    # Ensure that metadata_order only contains stuff in permute_within and block_data
    if(length(setdiff(metadata_order, union(names(permute_within), names(block_data)))) > 0) {
        stop("metadata_order contains metadata not in permute_within and block_data")
    }

    # Ensure that the data in permute_within matches that in dist
    ord <- rownames(as.matrix(D))
    if (length(ord) != nrow(permute_within) || length(blocks) != length(ord)) {
        stop("blocks, permute_within, and D are not the same size")
    }
    if (is.null(rownames(permute_within))) {
        warning("permute_within has no rownames - can't verify sample orders")
    } else if (!all(ord == rownames(permute_within))) {
        stop("rownames do not match between permute_within and D")
    }

    # Ensure matching between blocks and block_data
    if (any(is.na(blocks))) {
        stop("NAs are not allowed in blocks")
    }
    if (is.factor(blocks)) {
        if (any(!(levels(blocks) %in% rownames(block_data)))) {
            stop("not all block levels are contained in block_data")
        }
        # Match blocks with block_data and discard level information
        block_data <- block_data[match(levels(blocks), rownames(block_data)), , drop=F]
        blocks <- as.numeric(blocks)
    } else if (is.numeric(blocks)) {
        if (blocks < 1 || max(blocks) > nrow(block_data)) {
            stop("Numeric blocks has indices out of range")
        }
    } else if (is.character(blocks)) {
        if (is.null(rownames(block_data)) || !all(blocks %in% rownames(block_data))) {
            stop("blocks does not match the rownames of block_data")
        }
        # Transform to numeric
        blocks <- match(blocks, rownames(block_data))
    } else {
        stop("blocks must be a numeric, factor, or character vector")
    }

    # Error out on NA metadata rather than allowing adonis to error out with
    # a totally nonsensical error message
    na.removed <- 0
    if (any(is.na(permute_within)) || any(is.na(block_data))) {
        if (na.rm) {
            n_prerm <- length(blocks)

            # Remove NAs in block_data
            hasna <- (rowSums(is.na(block_data)) > 0) | (sapply(split(rowSums(is.na(permute_within)) > 0, blocks), mean) == 1)
            block_data <- block_data[!hasna,, drop=F]
            keep <- !hasna[blocks]
            blocks <- cumsum(!hasna)[blocks]

            blocks <- blocks[keep]
            permute_within <- permute_within[keep,, drop=F]
            D <- as.matrix(D)[keep, keep]
            # block_data is not subset, as the rows with NAs are no longer referenced in blocks

            # Remove NAs in permute_within
            keep <- rowSums(is.na(permute_within)) == 0
            blocks <- blocks[keep]
            permute_within <- permute_within[keep,, drop=F]
            D <- as.dist(D[keep, keep])

            if (length(blocks) < ncol(permute_within) + ncol(block_data)) {
                stop(sprintf("After omitting samples with NAs, the number of samples (%d) is less than the number of metadata (%d)",
                             length(blocks), ncol(permute_within) + ncol(block_data)))
            } else if (length(blocks) < n_prerm * 0.5) {
                warning(sprintf("Removed %d samples with NA metadata", n_prerm - length(blocks)))
            }
            na.removed <- n_prerm - length(blocks)
        } else {
            stop("Some metadata is NA! adonis does not support any NA in the metadata")
        }
    }
    # Warn on some suspicious input
    persample <- apply(permute_within, 1, function(x)is.factor(x) && !any(duplicated(x)))
    if (any(persample)) {
        warning(sprintf("%s in permute_within has one DOF per sample.", colnames(permute_within)[which(persample)[1]]))
    }
    if (length(unique(blocks)) < nrow(block_data)) {
        warning("Not all blocks have a sample associated with them. Block permutations will still be performed over the full set of blocks - if this is not desired, subset block_data to only the blocks which appear in the data.")
    }
    if (!any(duplicated(blocks))) {
        warning("blocks contains no duplicated elements")
    }

    # Test statistic from non-permuted data
    mtdat <- cbind(permute_within, block_data[blocks,,drop=F])
    ad <- adonis(D ~ ., permutations=0, data=mtdat[, metadata_order, drop=F])
    R2 <- ad$aov.tab$R2
    names(R2) <- rownames(ad$aov.tab)

    # Permutations
    nullsamples <- matrix(NA, nrow=length(R2), ncol=permutations)
    for (i in seq_len(permutations)) {
        within.i <- shuffle(nrow(permute_within), control=how(blocks=blocks))
        block.i <- sample(seq_len(nrow(block_data)))
        mtdat <- cbind(
            permute_within[within.i,,drop=F],
            block_data[block.i,,drop=F][blocks,,drop=F])
        perm.ad <- adonis(D ~ ., permutations=0, data=mtdat[, metadata_order, drop=F])

        nullsamples[,i] <- perm.ad$aov.tab$R2
    }

    # For residuals, test the other direction (i.e. p-value of all covariates)
    n <- length(R2)
    R2[n-1] <- 1 - R2[n-1]
    nullsamples[n-1,] <- 1 - nullsamples[n-1,]

    # P value calculation similar to adonis's
    exceedances <- rowSums(nullsamples > R2)
    P <- (exceedances + 1) / (permutations + 1)

    P[n] <- NA    # No p-values for "Total"
    ad$aov.tab$`Pr(>F)` <- P

    if (na.rm) {
        ad$na.removed <- na.removed
    }

    return (ad)
}

theme_nature <- function() list(
    theme_cowplot(),
    theme(
        text               = element_text(size=6),
        axis.text          = element_text(size=5),
        axis.title.x       = element_text(margin=margin(1, 0, 0.5, 0)),
        axis.title.x.top   = element_text(margin=margin(0, 0, 2, 0)),
        axis.title.y       = element_text(margin=margin(0, 1, 0, 0.5)),
        axis.title.y.right = element_text(margin=margin(0, 0, 0, 2)),
        axis.text.x        = element_text(margin=margin(1, 0, 0, 0)),
        axis.text.x.top    = element_text(margin=margin(0, 0, 1, 0)),
        axis.text.y        = element_text(margin=margin(0, 1, 0, 0)),
        axis.text.y.right  = element_text(margin=margin(0, 0, 0, 1)),
        axis.ticks         = element_line(size=0.3),
        axis.ticks.length  = unit(2, "pt"),
        axis.line          = element_line(size=0.3),
        axis.line.x        = element_line(size=0.3),
        axis.line.y        = element_line(size=0.3),
        line               = element_line(size=0.3),
        legend.margin      = margin(4, 4, 4, 4),
        legend.key.size    = unit(8, "pt"),
        legend.box.spacing = unit(4, "pt"),
        panel.spacing      = unit(1.5, "pt"),
        plot.title         = element_text(size=8),
        plot.margin        = margin(1, 1, 1, 1),
        strip.background   = element_blank(),
        strip.text         = element_text(size=6),
        strip.text.x       = element_text(margin=margin(3, 0, 3, 0)),
        strip.text.y       = element_text(margin=margin(0, 3, 0, 3))
    )
)


#' PERMANOVA results "heatmap"
#'
#' @param R2 A matrix of R2 values from \code{adonis}.
#' @param P A matrix of P-values from \code{adonis}.
#' @param fontsize The desired font size of the % R2 values in the heatmap.
#' @param FDR If \code{T}, P-values are first BH-adjusted, and significance
#' stars are shown as *** < 0.001, ** < 0.01, * < 0.05.
#' @return ggplot object.
PERMANOVA_heatmap <- function(R2, P, fontsize=6, FDR=T, alpha=NA, beta=NA, outfile=NA) {

    df <- melt(R2)
    colnames(df) <- c("Feature", "Dataset", "R2")
    df$P <- melt(P)[,3]
    if (FDR) {
        df$P <- p.adjust(df$P, method="fdr")
    }
    df$varExpPct <- sprintf("%.1f%%", 100*df$R2)
    df$varExpPct[is.na(df$R2)] <- ""
    df$NAtext <- ""
    df$NAtext[is.na(df$R2)] <- "N/A"
    df$stars <- cut(df$P, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), c("***", "**", "*", ""))
    df$stars[is.na(df$R2)] <- ""

    # Try to make a reasonable color scheme that has contrast where needed
    colors <- colorRampPalette(brewer.pal(n = 9, name = "Blues"))(100)
    R2Q <- quantile(R2, c(0.25, 0.75), na.rm=T)
    if (is.na(alpha) || is.na(beta)) {
    labhat <- optim(par=c(0, 0), method="Nelder-Mead",
                    fn=function(lab) sum((pbeta(c(0.25, 0.75), exp(lab[1]), exp(lab[2])) - R2Q)^2))
    abhat <- exp(labhat$par)
    colorvalues <- pbeta(seq(0, 1, length=length(colors)), abhat[1], abhat[2])

    # Label colors flip to white when the color is too dark
        df$lblcolor <- ifelse(qbeta(df$R2, abhat[1], abhat[2]) < 0.8, "black", "white")
        #cat(sprintf("Best-fit alpha = %g, beta = %g\n", abhat[1], abhat[2]))
    } else {
        colorvalues <- pbeta(seq(0, 1, length=length(colors)), alpha, beta)

        # Label colors flip to white when the color is too dark
        df$lblcolor <- ifelse(qbeta(df$R2, alpha, beta) < 0.8, "black", "white")
    }


    ggp <- ggplot(data=df, aes(x=Dataset, y=Feature)) +
        geom_tile(aes(fill=R2)) +
        geom_text(aes(label=varExpPct, color=lblcolor), size=fontsize/(14/5), nudge_y=-0.15) +
        geom_text(aes(label=NAtext), color="grey", size=fontsize/(14/5), nudge_y=-0.15) +
        geom_text(aes(label=stars, color=lblcolor), fontface="bold", size=1.25*fontsize/(14/5), nudge_y=0.12) +
        scale_fill_gradientn(colors=colors, values=colorvalues, limits=c(0, 1), na.value="white") +
        scale_color_manual(values=c(black="black", white="white")) +
        scale_x_discrete(expand=c(0,0)) + xlab(NULL) +
        scale_y_discrete(expand=c(0,0), position = "right", limits = rev(levels(df$Feature))) + ylab(NULL) +
        guides(color="none",
               fill=guide_colourbar(title=NULL, barheight=unit(40,"mm"), label.position = "left")) +
        theme_nature() +
        theme(axis.text.x = element_text(angle=-17, hjust=0),
              panel.border=element_rect(fill=NA),
              legend.position = "left", axis.ticks.y = element_blank())

    df %>% as_tibble() %>% write_tsv(outfile)

    return (ggp)
}

merge_metadata <- function(data, metadata) {
    # modified from original to be generic

    # check for samples without metadata
    extra_feature_samples <-
        setdiff(rownames(data), rownames(metadata))
    if (length(extra_feature_samples) > 0)
        sprintf(
            paste("The following samples were found",
                "to have features but no metadata.",
                "They will be removed. %s"),
            paste(extra_feature_samples, collapse = ",")
        )
        
    # check for metadata samples without features
    extra_metadata_samples <-
        setdiff(rownames(metadata), rownames(data))
    if (length(extra_metadata_samples) > 0)
        sprintf(
            paste("The following samples were found",
                "to have metadata but no features.",
                "They will be removed. %s"),
            paste(extra_metadata_samples, collapse = ",")
        )
        
    # get a set of the samples with both metadata and features
    intersect_samples <- intersect(rownames(data), rownames(metadata))
    sprintf(
        "A total of %s samples were found in both the data and metadata",
        length(intersect_samples)
    )
        
    # now order both data and metadata with the same sample ordering
    sprintf(
        "Reordering data/metadata to use same sample ordering")

    data <- data[intersect_samples, , drop = FALSE]
    metadata <- metadata[intersect_samples, , drop = FALSE]

    pcl <- data.frame(matrix(NA, nrow = nrow(data)))
    pcl$meta <- metadata
    pcl$x <- as.matrix(data)

    return (pcl)
}


bc_omnibus_tests <- function(pcl, meta, covariates, perindividual_covariates, blocks_off=FALSE, method="bray", Nperms=4999) {
    # Main function for performing overall and per-covariate PERMANOVAs for a given dataset

    # Merge in relevant metadata
    pcl <- merge_metadata(pcl, meta)

    # Subset to samples with no NA's (or adonis will puke with a bizarre error message)
    covars.unl <- unique(unlist(covariates))
    ok_covars <-
        ((colMeans(colwise(duplicated)(pcl$meta[,covars.unl])) > 0.75) | (colSums(colwise(is.numeric)(pcl$meta[,covars.unl])) == 1)) # Can't be all unique values
    # Ensure not too many NAs
    ok_covars <- ok_covars & (colMeans(is.na(pcl$meta[,covars.unl])) < 0.5)

    # Get the distance matrix
    D <- vegdist(pcl$x, method, binary=method=="jaccard")

    # Find relevant covariates
    ok_covars[ok_covars] <- apply(pcl$meta[,covars.unl[ok_covars],drop=F], 2, function(x)!all(x==na.omit(x)[1], na.rm=T))
    overall_covars <- covars.unl[ok_covars]
    covars <- union(overall_covars, "subject")

    # Make a per-person covariate table
    pcl$meta$subject <- factor(pcl$meta$subject)
    subj_meta <- pcl$meta[order(pcl$meta$subject),,drop=F]
    subj_meta <- subj_meta[!duplicated(subj_meta$subject), covars[covars %in% perindividual_covariates], drop=F]
    rownames(subj_meta) <- levels(pcl$meta$subject)

    # Do the overall PERMANOVA
    overall.ad <- PERMANOVA_repeat_measures(
        D, permutations=Nperms, na.rm=T,
        permute_within=pcl$meta[,overall_covars, drop=F])

    # Modify the main PERMANOVA to match the expected covariate structure
    ad <- overall.ad
    ad$aov.tab <- ad$aov.tab[seq_len(length(covariates)+1),]
    ad$aov.tab[,] <- NA
    rownames(ad$aov.tab) <- c(names(covariates), "Overall")
    ad$aov.tab["Overall",] <- colSums(overall.ad$aov.tab[seq_len(nrow(overall.ad$aov.tab)-2),])
    ad$aov.tab["Overall","Pr(>F)"] <- overall.ad$aov.tab$`Pr(>F)`[nrow(overall.ad$aov.tab)-1]

    ad$N <- pcl$ns

    # Do per-covariate tests
    blocks <- factor(pcl$meta$subject)
    for (var in covariates) {
        varsub <- covariates[[var]]
        if (any(ok_covars[varsub])) {
            if (any(!ok_covars[varsub])) {
                varsub <- varsub[ok_covars[varsub]]
            }

            if (("subject" %in% varsub) || blocks_off){
                # Subject is among the covariates, so the permutations must be free
                sep.ad <- PERMANOVA_repeat_measures(
                    D, permutations=Nperms, na.rm=T,
                    permute_within=pcl$meta[,varsub, drop=F])
            } else {
                # Per-covariate test with subject-aware permutations
                sep.ad <- PERMANOVA_repeat_measures(
                    D, permutations=Nperms, na.rm=T,
                    permute_within=pcl$meta[,varsub[!(varsub %in% perindividual_covariates)], drop=F],
                    blocks=blocks,
                    block_data=subj_meta[,varsub[varsub %in% perindividual_covariates], drop=F])
            }

            # Store results
            ad$aov.tab[var,] <- colSums(sep.ad$aov.tab[seq_len(nrow(sep.ad$aov.tab)-2),])
            ad$aov.tab[var,"Pr(>F)"] <- sep.ad$aov.tab$`Pr(>F)`[nrow(sep.ad$aov.tab)-1]
        }
    }

    return (ad)
}

# Add command line arguments #
options <- optparse::OptionParser(
    usage = paste("%prog [options]", " <data.tsv> ", " <metadata.tsv> ", " <output_file.png> ")
)

options <- optparse::add_option(options,
    c("-p", "--permutations"),
    type = "integer",
    dest = "nperms",
    default = 4999,
    help = paste0("The minimum number of permutations",
    " [ Default: %default ]"
    )
)

options <- optparse::add_option(options,
    c("-i", "--static_covariates"),
    type = "character",
    dest = "static_covariates",
    default = NULL,
    help = paste0("Treated as per-individual and permuted across individuals instead of longitudinally",
    " [ Default: All covariates ]"
    )
)

options <- optparse::add_option(options,
    c("-s", "--scale"),
    type = "integer",
    dest = "scale",
    default = 1,
    help = paste0("The scale to apply to the data",
    " [ Default: %default ]"
    )
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

## allow for multiple input data files
datatype_list <- c()
for (datafile in unlist(strsplit(positional_args[1], ",", fixed = TRUE))) {

    datatype <- unlist(strsplit(basename(datafile),"_features.txt"))[1]
    datatype_list <- append(datatype_list, datatype)

    data <- t(read.table(datafile, header = TRUE, row.names = 1, sep="\t", comment.char = ""))
    metadata <- data.frame(read.table(positional_args[2], header = TRUE, row.names = 1, sep="\t", comment.char = ""))

    # check for samples as columns or rows
    samples_rows <- intersect(rownames(metadata),rownames(data))
    if (length(samples_rows) < 1) {
        metadata <- type.convert(as.data.frame(t(metadata)))
    }

    covariates <- as.list(colnames(metadata))
    names(covariates) <- covariates

    # Set the individual covariates
    if (is.null(current_args$static_covariates)) {
        blocks_off = TRUE
        current_args$static_covariates <- covariates
    } else {
        blocks_off = FALSE
        current_args$static_covariates <- unlist(strsplit(current_args$static_covariates, ",", fixed = TRUE))
    }

    # add subject based on sample id if subject is not included
    if (!"subject" %in% colnames(metadata)) {
        metadata$subject <- rownames(metadata)
    }

    # apply scale
    data <- data/current_args$scale

    # Filter by abundance using zero as value for NAs
    total_samples <- nrow(metadata)
    min_samples <- total_samples * current_args$min_prevalence

    data_zeros <- data
    data_zeros[is.na(data_zeros)] <- 0
    filtered_data <- data[,colSums(data_zeros > current_args$min_abundance) > min_samples, drop = FALSE]

    # Filter the metadata to remove any samples without data
    metadata_zeros <- metadata
    metadata_zeros[metadata_zeros == "UNK"] <- NA
    filtered_metadata <- metadata[rowSums(metadata_zeros != 0, na.rm=TRUE) > 0, , drop = FALSE]

    ad <- bc_omnibus_tests(filtered_data,filtered_metadata,covariates,current_args$static_covariates,blocks_off,Nperms=current_args$nperms)

    if (! (exists("R2"))) {
        R2 <- matrix(0, nrow=length(covariates) + 1, ncol=0)
        P <- matrix(0, nrow=length(covariates) + 1, ncol=0)
    }

    R2 <<- cbind(R2, ad$aov.tab$R2)
    row.names(R2) <- rownames(ad$aov.tab)
    P <<- cbind(P, ad$aov.tab$`Pr(>F)`)
    row.names(P) <- rownames(ad$aov.tab)

    # generate the stacked barplot
    univar_tax = rbind(ad$aov.tab$`Pr(>F)`, ad$aov.tab$R2)
    univar_tax = as.data.frame(t(univar_tax))
    names(univar_tax) = c("P-Value", "R2")

    univar_tax$stars = cut(univar_tax$`P-Value`, c(0, 0.001, 0.01, 0.05, 0.1, 1), labels = c("***", "**", "*", "`", ""))

    univar_tax$`P-Value` = round(univar_tax$`P-Value`, 3)
    univar_tax$R2 = univar_tax$R2 *100

    row.names(univar_tax) <- rownames(ad$aov.tab)

    dodge = position_dodge(width = 0.8)

    plot <- ggplot(data = univar_tax, aes(reorder(row.names(univar_tax), univar_tax$R2), y = R2, label = univar_tax$`P-Value`)) + geom_bar(stat = "identity", position = "identity", fill = "#800000") + geom_text(position = dodge, vjust = 0.5, hjust = -0.1, size = 3) + theme_bw(base_size = 12) + ylab("Univariable R-squared") + coord_flip() + ylim(0, as.integer(max(univar_tax$R2))*1.05) + xlab("") + labs(fill = "")

    plot_file <- gsub(".png",paste("_",datatype,".png",sep=""),positional_args[3])
    print(paste("Writing to plot file",plot_file))
    png(plot_file, res = 150, height = 800, width = 1100)
    print(plot)
    dev.off()

    text_file <- gsub(".png",paste("_",datatype,".txt",sep=""),positional_args[3]) 
    univar_tax %>% as_tibble() %>% write_tsv(text_file)
}

# update the column names
colnames(R2) <- datatype_list
colnames(P) <- datatype_list

alpha <- 1.9
beta <- 0.1

ggp_all <- PERMANOVA_heatmap(R2, P, fontsize=5, FDR=T, alpha=alpha, beta=beta, outfile=gsub(".png",".txt",positional_args[3]))
png(positional_args[3], res = 150, height = 800, width = 1100)
print(ggp_all)
dev.off()


