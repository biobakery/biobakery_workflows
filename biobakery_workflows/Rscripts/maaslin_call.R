#!/usr/bin/env Rscript

library(tidyverse)
library(vegan)
library(cowplot)
library(readr)
library(reshape2)
library(readxl)
library(ape)
library(gdata)
library(RColorBrewer)
library(viridis)
library(Maaslin2)

###########################################################################################################
#####  Implements prep for Maaslin analysis and calls Maaslin2
#####  Requires:
#####   -metadata.tsv file (sample ids column should be named 'sample')
#####   -closed reference either biom or tsv file 
#####   -closed_reference.tre fasttree file (optional, if available)
#####  Arguments to script are:
#####   -input_dir  (folder where required files are)
#####   -normalize  (if 'no' it will skip relative abundance calculation, by default it calculates)
#############################################################################################################

## Collect arguments
args <- commandArgs(TRUE)

## Parse arguments (we expect the form --arg=value)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
args.df <- as.data.frame(do.call("rbind", parseArgs(args)))

args.list <- as.list(as.character(args.df$V2))
names(args.list) <- args.df$V1

args.list$input_dir <- "~/dada_demo_output_vsearch"
args.list$normalize <- "yes"
args.list$workflow <- "wmgx"
# args.list$workflow <- "16s"
args.list$exclude_metadata <- c("compliance_75", "weight_kg", "cGMPpmolmL", "TotalFreeThiolsuM", "totalCholesterol", "HDLC_", "LDLC_", "Trigl_", "Chol_HDL_Ratio", "HDLP", "LDLP", "ApoA1", "ApoB", "ApoBApoA1", "urinary_syringic_acid", "urinary_chlorogenic_acid", "urinary_isovanillicacid_3_glu", "urinary_m_coumaric_acid", "urinary_DOPAC", "urinary_dihydrocaffeic_acid", "urinary_vanillic_acid", "urinary_hippuric_acid", "urinary_Homovanillic_acid", "urinary_3PHBA_4sul_4PHBA_3sul", "total_serum_derived_metabolites", "total_urine_derived_metabolites", "framingham_CV_score", "CRP_ug_ml_", "IL6_pg_ml_", "TNF__pg_ml_")
args.list$exclude_subjects <- c("671", "1025", "1074", "1098", "1131", "1170", "1355", "1392", "1479", "1513", "1704")
args.list$exclude_samples <- c("VOL1151_0M", "VOL1151_6M", "VOL1280_0M", "VOL1280_6M", "VOL1384_6M", "VOL631_0M", "VOL631_6M", "VOL1267_0M", "VOL1267_6M", "VOL1389_0M", "VOL1389_6M")
args.list$metadata_categorical = c("Age","sex", "group", "timepoint","metabolicsyndromecriteria","statin","bp_med")


#location of files
setwd(args.list$input_dir)

# Taxonomic Profile -----------------------------------------------------------------------------------
metadata <- read.delim(args.list$metadata, header = TRUE, sep="\t", fileEncoding="Latin1")
tax_profile <- read.delim(args.list$abundance, header = TRUE, sep="\t")

rownames(metadata) <- metadata[,1]
metadata <- metadata[,-1]

rownames(tax_profile) <- tax_profile[,1]
subset_ind <- grep( "s__|^unclassified", rownames(tax_profile), invert = F )
unstrat_taxa_t <- tax_profile[subset_ind,]
subset_ind <- grep( "t__", rownames(unstrat_taxa_t), invert = T )
unstrat_taxa <- unstrat_taxa_t[subset_ind,]

# Normalize
sz=dim(unstrat_taxa)
unstrat_taxa <- unstrat_taxa[-sz[2]]
sz=dim(unstrat_taxa)
csums<-colSums(unstrat_taxa[2:sz[1],2:sz[2]])
taxa_profile <- as.data.frame(t(unstrat_taxa[2:sz[1],2:sz[2]])/csums)

rownames(taxa_profile) <- gsub("_taxonomic_profile","",rownames(taxa_profile))

merged_taxa_data <- merge(metadata,taxa_profile, by="row.names", sort=TRUE)

Maaslin2(taxa_profile, metadata, args.list$maaslin_output)