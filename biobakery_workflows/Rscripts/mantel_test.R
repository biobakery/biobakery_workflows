#!/usr/bin/env Rscript

# A utility script to compute and display a mantel test.
# Based on stats visualizations and code from the hmp2_analysis repository.

library(vegan)
library(ggplot2)
library(cowplot)
library(optparse)
library(gridExtra)

distance_method <- c(
    "Taxonomy"    = "bray",
    "KOs (DNA)"   = "bray",
    "KOs (RNA)"   = "bray",
    "KOs (Protein)"="bray",
    "Metabolites" = "bray",
    "Viromics"    = "jaccard",
    "Biopsy 16S"  = "bray",
    "Biopsy HTX"  = "bray",
    "Diet"        = "manhattan"
)

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

# Add command line arguments #
options <- optparse::OptionParser(
    usage = paste("%prog [options]", " <data.tsv> ", " <data2.tsv> ", " <output_file.png> ")
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
    stop(paste("Please provide the required","positional arguments","<data1.tsv> <data2.tsv> <output_file.png>"))
}

# Input variables
## assume data file has samples as columns
## allow for pounds in header lines
data <- t(read.table(positional_args[1], header = TRUE, row.names = 1, sep="\t", comment.char = ""))
data2 <- t(read.table(positional_args[2], header = TRUE, row.names = 1, sep="\t", comment.char = ""))

# check for samples as columns or rows
samples_rows <- intersect(rownames(data2),rownames(data))
if (length(samples_rows) < 1) {
    # allow for special chars in names
    original_data2_rownames <- rownames(data2)
    rownames(data2) <- make.names(rownames(data2))
    samples_rows <- intersect(rownames(data2),rownames(data))
    if (length(samples_rows) < 1) {
        samples_col <- intersect(rownames(data2,colnames(data)))
        if (length(samples_col) < 1){
            rownames(data2) <- original_data2_rownames
            sample <- colnames(data2)
            data2 <- as.data.frame(t(data2))
        } else {
            sample <- colnames(data)
            data <- as.data.frame(t(data))
        }
    }
}

# Filter by abundance using zero as value for NAs
total_samples <- nrow(data2)
min_samples <- total_samples * current_args$min_prevalence

data_zeros <- data
data_zeros[is.na(data_zeros)] <- 0
filtered_data <- data[,colSums(data_zeros > current_args$min_abundance) > min_samples, drop = FALSE]

# remove empty rows from data
filtered_data <- filtered_data[rowSums(filtered_data != 0, na.rm=TRUE) > 0, , drop = FALSE]

if ((ncol(filtered_data) < 1) || (nrow(filtered_data) < 1)) {
  stop("No data remain in the data after filtering for min abundance and prevalence")
}

# Filter the data2 to remove any samples without data or data2 with more than max values missing
data2_zeros <- data2
data2_zeros[is.na(data2_zeros)] <- 0
filtered_data2 <- data[,colSums(data2_zeros > current_args$min_abundance) > min_samples, drop = FALSE]

# remove empty rows from data
filtered_data2 <- filtered_data2[rowSums(filtered_data2 != 0, na.rm=TRUE) > 0, , drop = FALSE]

if ((ncol(filtered_data2) < 1) || (nrow(filtered_data2) < 1)) {
  stop("No data2 remain in the data after filtering for min abundance and prevalence")
}

# filter data and data2 to only include the same samples in the same order
sorted_samples <- sort(intersect(row.names(filtered_data2), row.names(filtered_data)))
filtered_data2 <- filtered_data2[sorted_samples, , drop = FALSE]
filtered_data <- filtered_data[sorted_samples, , drop = FALSE]


##############################################################
#                                                            #
# Functions from mantel_test.R from hmp2_analysis repository #
# (modified to make generic)                                 #
##############################################################

# function from mantel_test.R from hmp2_analysis repository
interindividual_mantel_test_doa <- function(pcl1, method1, pcl2, method2, Nperms, Nbs=Nperms) {

    library(ade4)
    library(vegan)
    D1 <- vegdist(pcl1, method=method1, binary=method1=="jaccard")
    D2 <- vegdist(pcl2, method=method2, binary=method2=="jaccard")

    mantelbootstrap <- function(m1, m2, nrepet) {
        s1 <- unclass(m1)
        s2 <- unclass(m2)
        permi <- matrix(1:nrepet, nrow = nrepet, ncol = 1)
        bss <- apply(permi, 1, function(i) {
            bsi <- sample(seq_along(s1), length(s1), replace=T)
            return (cor(s1[bsi], s2[bsi]))
        })
        return (bss)
    }

    mt <- mantel.rtest(D1, D2, nrepet=Nperms)
    mt$bootstraps <- mantelbootstrap(D1, D2, nrepet=Nbs)
    return (mt)
}

mt <- interindividual_mantel_test_doa(
      filtered_data, distance_method["Taxonomy"],
      filtered_data2, distance_method["Taxonomy"],
      Nperms=current_args$nperms)
#mt_inter_doa$C[i,j] <- mt$obs
#mt_inter_doa$Cil[i,j] <- quantile(mt$bootstraps, 0.025, na.rm=T)
#mt_inter_doa$Ciu[i,j] <- quantile(mt$bootstraps, 0.975, na.rm=T)
#mt_inter_doa$P[i,j] <- mt$pvalue

# Plots
library(RColorBrewer)
manteltest_plot <- function(O, P, Ocil, Ociu, title, stars=F) {

    library(reshape2)
    library(viridis)
    df <- melt(O)
    colnames(df) <- c("V1", "V2", "obs")
    df$V2 <- factor(df$V2, levels=rev(levels(df$V2)))
    df$obs <- df$obs^2
    if (!missing(P)) {
        Padj <- matrix(p.adjust(P, method="fdr"), nrow=nrow(P))
        df$padj <- melt(Padj)$value
    }
    if (!missing(Ocil)) {
        df$cil <- melt(Ocil)$value^2
        df$cil[melt(Ocil)$value < 0] <- 0
        df$ciu <- melt(Ociu)$value^2
    }

    # Increase the contrast of the color scale where it matters
    alpha <- 1.9
    beta <- 0.1
    colorvalues <- pbeta(seq(0, 1, length=101), alpha, beta)
    df$lblcolor <- ifelse(qbeta(df$obs, alpha, beta) < 0.8, "black", "white")

    ggp <- ggplot(data=df, aes(x=V1, y=V2)) +
        geom_tile(aes(fill=obs)) +
        scale_fill_gradientn(colours=colorRampPalette(brewer.pal(n = 9, name = "Blues"))(100),
                             values=colorvalues,
                             na.value="white", limits=c(0, 1), name="Variance Explained")
    #ggp
    nudge_y <- 0
    if (!missing(P)) {
        ggp <- ggp + if (stars) {
            nudge_y <- -0.14
            geom_text(aes(label=ifelse(padj<=0.001, "***", ifelse(padj<=0.01, "**", ifelse(padj<=0.05, "*", ""))),
                          color=lblcolor),
                      size=2, nudge_y=0.15, fontface="bold")
        } else {
            geom_text(aes(label=ifelse(is.na(padj), "", sprintf("FDR p\n%.2g", padj)),
                          color=lblcolor), size=1.7)
        }
    }
    ggp <- ggp + if (!missing(Ocil)) {
        geom_text(aes(label=ifelse(is.na(obs), "", sprintf("%.1f%%\n[%.1f%% - %.1f%%]", 100*obs, 100*cil, 100*ciu)),
                      color=lblcolor),
                  size=1.7, nudge_y=nudge_y)
    } else {
        geom_text(aes(label=ifelse(is.na(obs), "", sprintf("%.1f%%", 100*obs)),
                      color=lblcolor),
                  size=1.7, nudge_y=nudge_y)
    }
    ggp <- ggp +
        geom_text(aes(label=ifelse(V1!=V2, ifelse(is.na(obs), "N/A", ""), "")), size=1.7, color="gray") +
        theme_nature() +
        theme(axis.text.x=element_text(angle=-17, hjust=0.05),
              legend.position = "left", axis.ticks.y = element_blank()) +
        scale_x_discrete(expand=c(0, 0)) + scale_y_discrete(expand=c(0, 0), position = "right") +
        scale_color_manual(values=c(white="white", black="black")) +
        guides(fill=guide_colourbar(title=NULL, barheight=unit(0.65,"npc"), label.position = "left"), color="none") +
        xlab(NULL) + ylab(NULL)
    
    return (ggp)
}

#png(positional_args[3], res = 150, height = 800, width = 1100)
#print(manteltest_plot(mt_inter_doa$C, t(mt_inter_doa$P), mt_inter_doa$Cil, mt_inter_doa$Ciu, title="Inter-individual (DOA)"))
#dev.off()

