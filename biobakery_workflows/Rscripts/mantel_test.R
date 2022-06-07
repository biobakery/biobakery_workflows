#!/usr/bin/env Rscript

# A utility script to compute and display a mantel test.
# Based on stats visualizations and code from the hmp2_analysis repository.

library(vegan)
library(ggplot2)
library(cowplot)
library(optparse)
library(gridExtra)
library(tibble)
library(dplyr)
library(readr)

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
        text               = element_text(size=8),
        axis.text          = element_text(size=8),
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
        legend.margin      = margin(4,4,4,4),
        legend.key.size    = unit(8, "pt"),
        legend.box.spacing = unit(4, "pt"),
        panel.spacing      = unit(1.5, "pt"),
        plot.title         = element_text(size=8),
        plot.margin        = margin(5, 5, 5, 5),
        strip.background   = element_blank(),
        strip.text         = element_text(size=8),
        strip.text.x       = element_text(margin=margin(3, 0, 3, 0)),
        strip.text.y       = element_text(margin=margin(0, 3, 0, 3))
    )
)

# Add command line arguments #
options <- optparse::OptionParser(
    usage = paste("%prog [options]"," <metadata.tsv> ", " <data.tsv> ", " <output_file.png> ")
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
    stop(paste("Please provide the required","positional arguments","<metadata.tsv> <data1.tsv> <output_file.png>"))
}

metadata_file <- positional_args[1]
datafiles <- unlist(strsplit(positional_args[2], ",", fixed = TRUE))

M <- matrix(NA, nrow=length(datafiles), ncol=length(datafiles))

mt_intra <- list(C=M, Cil=M, Ciu=M, P=M)
mt_inter_aod <- mt_intra
mt_inter_doa <- mt_intra

datatype_list <- c()
input_list <- vector(mode = "list", length = length(datafiles))
input_index <- 1

### assume all input files have the samples as the column names
for (datafile in datafiles) {

    datatype <- unlist(strsplit(basename(datafile),"_features.txt"))[1]
    datatype_list <- append(datatype_list, datatype)

    data <- t(read.table(datafile, header = TRUE, row.names = 1, sep="\t", comment.char = "", quote = ""))
    samples <- row.names(data)

    # Filter by abundance using zero as value for NAs
    total_samples <- length(samples)
    min_samples <- total_samples * current_args$min_prevalence

    data_zeros <- data
    data_zeros[is.na(data_zeros)] <- 0
    filtered_data <- data[,colSums(data_zeros > current_args$min_abundance) > min_samples, drop = FALSE]

    # remove empty rows from data
    filtered_data <- filtered_data[rowSums(filtered_data != 0, na.rm=TRUE) > 0, , drop = FALSE]

    if ((ncol(filtered_data) < 1) || (nrow(filtered_data) < 1)) {
      stop("No data remain in the data after filtering for min abundance and prevalence")
    }

    input_list[[input_index]] <- filtered_data
    input_index <- input_index + 1
}

rownames(M) <- datatype_list
colnames(M) <- datatype_list

##############################################################
#                                                            #
# Functions from mantel_test.R from hmp2_analysis repository #
# (modified to make generic)                                 #
##############################################################

# function from mantel_test.R from hmp2_analysis repository
intraindividual_mantel_test <- function(pcl1, method1, pcl2, method2, metadata, Nperms=999, Nbs=Nperms) {

    library(vegan)
    D1 <- vegdist(pcl1, method=method1, binary=method1=="jaccard")
    D2 <- vegdist(pcl2, method=method2, binary=method2=="jaccard")

    # Code based on ade4:::mantel.rtest

    library(permute)
    library(ade4)

    permhow <- how(blocks=metadata$subject)
    intrasubjectD <- unclass(as.dist(outer(
        as.character(metadata$subject),
        as.character(metadata$subject), FUN="=="))) > 0

    if (sum(intrasubjectD) < 100) {
        return (list(obs = NA, pvalue = NA, bootstraps = NA))
    }

    nsamples <- nrow(metadata)
    permutedist <- function(m, i) {
        w0 <- permute(i, nsamples, permhow)
        m <- as.matrix(m)
        return(as.dist(m[w0, w0]))
    }
    mantelnoneuclid <- function(m1, m2, nrepet) {
        obs <- cor(unclass(m1)[intrasubjectD], unclass(m2)[intrasubjectD])
        if (nrepet == 0)
            return(obs)
        permi <- matrix(1:nrepet, nrow = nrepet, ncol = 1)
        perm <- apply(permi, 1, function(i)
            cor(unclass(m1)[intrasubjectD],
                unclass(permutedist(m2, i))[intrasubjectD]))
        w <- as.randtest(obs = obs, sim = perm, call = match.call(),
                         subclass = "mantelrtest")
        return(w)
    }

    mantelbootstrap <- function(m1, m2, nrepet) {
        s1 <- unclass(m1)[intrasubjectD]
        s2 <- unclass(m2)[intrasubjectD]
        permi <- matrix(1:nrepet, nrow = nrepet, ncol = 1)
        bss <- apply(permi, 1, function(i) {
            bsi <- sample(seq_along(s1), length(s1), replace=T)
            return (cor(s1[bsi], s2[bsi]))
        })
        return (bss)
    }
    mt <- mantelnoneuclid(D1, D2, nrepet=Nperms)
    mt$bootstraps <- mantelbootstrap(D1, D2, nrepet=Nbs)
    return (mt)
}


# function from mantel_test.R from hmp2_analysis repository
interindividual_mantel_test_aod <- function(pcl1, method1, pcl2, method2, Nperms, Nbs=Nperms) {

    library(ade4)
    library(vegan)
    D1 <- vegdist(pcl1, method=method1, binary=method1=="jaccard")
    D2 <- vegdist(pcl2, method=method2, binary=method2=="jaccard")

    # Build matrices of the averages-of-dissimilarities
    subject <- pcl1$subject
    unqSubj <- unique(subject)
    Da1 <- matrix(0, length(unqSubj), length(unqSubj))
    rownames(Da1) <- unqSubj
    colnames(Da1) <- unqSubj
    Da2 <- Da1
    D1 <- as.matrix(D1)
    D2 <- as.matrix(D2)
    for (i in seq_along(unqSubj)) {
        for (j in seq_along(unqSubj)) {
            if (i > j) {
                mu1 <- mean(D1[subject==unqSubj[i], subject==unqSubj[j]])
                Da1[i,j] <- mu1
                Da1[j,i] <- mu1
                mu2 <- mean(D2[subject==unqSubj[i], subject==unqSubj[j]])
                Da2[i,j] <- mu2
                Da2[j,i] <- mu2
            }
        }
    }
    Da1 <- as.dist(Da1)
    Da2 <- as.dist(Da2)

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

    mt <- mantel.rtest(Da1, Da2, nrepet=Nperms)
    mt$bootstraps <- mantelbootstrap(Da1, Da2, nrepet=Nbs)
    return (mt)
}

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

metadata <- data.frame(read.table(positional_args[1], header = TRUE, row.names = 1, sep="\t", comment.char = ""))

samples_rows <- intersect(rownames(metadata),rownames(input_list[[1]]))
if (length(samples_rows) < 1) {
    metadata <- type.convert(as.data.frame(t(metadata)))
}

for (i in seq_along(datafiles)) {
    for (j in seq_along(datafiles)) {
        if ( i != j && is.na(mt_inter_doa$C[i,j]) && is.na(mt_inter_doa$C[j,i])){

            # filter data and data2 to only include the same samples in the same order
            sorted_samples <- sort(intersect(row.names(input_list[[i]]), row.names(input_list[[j]])))

            filtered_data <- input_list[[i]][sorted_samples, , drop = FALSE]
            filtered_data2 <- input_list[[j]][sorted_samples, , drop = FALSE]
            metadata_sorted <- metadata[sorted_samples, , drop = FALSE]

            #if ("subject" %in% colnames(metadata)) {

            # mt <- intraindividual_mantel_test(
            #      filtered_data, distance_method["Taxonomy"],
            #      filtered_data2, distance_method["Taxonomy"],
            #      metadata_sorted,
            #      Nperms=current_args$nperms)

            #  mt_intra$C[i,j] <- mt$obs
            #  mt_intra$Cil[i,j] <- quantile(mt$bootstraps, 0.025, na.rm=T)
            #  mt_intra$Ciu[i,j] <- quantile(mt$bootstraps, 0.975, na.rm=T)
            #  mt_intra$P[i,j] <- mt$pvalue
            # }

            mtd <- interindividual_mantel_test_doa(
                  filtered_data, distance_method["Taxonomy"],
                  filtered_data2, distance_method["Taxonomy"],
                  Nperms=current_args$nperms)

            mt_inter_doa$C[i,j] <- mtd$obs
            mt_inter_doa$Cil[i,j] <- quantile(mtd$bootstraps, 0.025, na.rm=T)
            mt_inter_doa$Ciu[i,j] <- quantile(mtd$bootstraps, 0.975, na.rm=T)
            mt_inter_doa$P[i,j] <- mtd$pvalue
        }
    }
}

# Plots
library(RColorBrewer)
manteltest_plot <- function(O, P, Ocil, Ociu, datatype_list, title) {

    library(reshape2)
    library(viridis)
    colnames(O) <- datatype_list
    rownames(O) <- datatype_list
    df <- melt(O)
    colnames(df) <- c("V1", "V2", "obs")

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
    alpha <- 1.2
    beta <- 0.1
    colorvalues <- pbeta(seq(0, 1, length=101), alpha, beta)

    ggp <- ggplot(data=df, aes(x=V1, y=V2)) +
        scale_fill_gradientn(colours=colorRampPalette(brewer.pal(n = 9, name = "Blues"))(100),
                             values=colorvalues,
                             na.value="white", limits=c(0, 1), name="Variance Explained") +
        geom_tile(aes(fill=obs)) 

    #ggp
    ggp <- ggp + geom_text(aes(label=ifelse(is.na(padj), "", sprintf("FDR p\n%.2g", padj))), color="black", size=2.5)
    ggp <- ggp + geom_text(aes(label=ifelse(is.na(obs), "", sprintf("%.1f%%\n[%.1f%% - %.1f%%]", 100*obs, 100*cil, 100*ciu))), color="white", size=2.5)
    ggp <- ggp + xlab(NULL) + ylab(NULL)
    ggp <- ggp + theme_nature()
    ggp <- ggp + geom_text(aes(label=ifelse(V1==V2, ifelse(is.na(obs), "N/A", ""), "")), size=2.5, color="gray") 
    ggp <- ggp + theme(axis.text.x=element_text(angle=-17, hjust=0.05), legend.position = "left")  
    ggp <- ggp + guides(fill=guide_colourbar(title=NULL, barheight=unit(0.65,"npc"), label.position = "left"), color="none")
  
    return (ggp)
}

text_file <- sub(".png",".txt", positional_args[3])
png(positional_args[3], res = 150, height = 800, width = 1100)

mt_inter_doa_chars <- data.frame(lapply(mt_inter_doa, as.character))
if ("subject" %in% colnames(metadata)) {
  #print(manteltest_plot(mt_intra$C, t(mt_intra$P), mt_intra$Cil, mt_intra$Ciu, datatype_list, title="Intra-individual"))
  #mt_intra %>% as_tibble() %>% write_tsv(text_file)
  print(manteltest_plot(mt_inter_doa$C, t(mt_inter_doa$P), mt_inter_doa$Cil, mt_inter_doa$Ciu, datatype_list, title="Inter-individual (DOA)"))
  try({ write_tsv(mt_inter_doa_chars,text_file) }, silent = TRUE)
} else {
  print(manteltest_plot(mt_inter_doa$C, t(mt_inter_doa$P), mt_inter_doa$Cil, mt_inter_doa$Ciu, datatype_list, title="Inter-individual (DOA)"))
  try({ write_tsv(mt_inter_doa_chars,text_file) }, silent = TRUE)
}

dev.off()

