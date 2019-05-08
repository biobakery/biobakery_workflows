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

##  Implements prep for Maaslin analysis and calls Maaslin2
##  Requires:
##  -metadata file 
##  -abundance file 
##  -input_dir  input direcory
##  -maaslin_output output_directory
##  Optional:
##  -normalize
##  -metadata_exclude
##  -subjects_exclude
##  -samples_exclude
##  -metadata_categorical
##  -Maaslin options


# Collect arguments
args <- commandArgs(TRUE)

# Parse arguments (we expect the form --arg=value)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
args.df <- as.data.frame(do.call("rbind", parseArgs(args)))

args.list <- as.list(as.character(args.df$V2))
names(args.list) <- args.df$V1

if(is.null(args.list$normalize)){
 args.list$normalize <- "yes"
}

#args.list$input_dir <- "/Users/anamailyan/dada_demo_output_vsearch"
#args.list$output_dir <- "/Users/anamailyan/dada_demo_output_vsearch_stats"
#args.list$maaslin_output <- "/Users/anamailyan/dada_demo_output_vsearch_stats/maaslin_path"
#args.list$metadata <- "/Users/anamailyan/dada_demo_output_vsearch/metadata.tsv"
#args.list$abundance <- "/Users/anamailyan/dada_demo_output_vsearch/humann2/merged/pathabundance_relab.tsv"
#args.list$core <- 4
#args.list$min_abundance <- 0.01


#args.list$exclude_metadata <- "compliance_75, weight_kg, cGMPpmolmL, TotalFreeThiolsuM, totalCholesterol, HDLC_, LDLC_, Trigl_, Chol_HDL_Ratio, HDLP, LDLP, ApoA1, ApoB, ApoBApoA1, urinary_syringic_acid, urinary_chlorogenic_acid, urinary_isovanillicacid_3_glu, urinary_m_coumaric_acid, urinary_DOPAC, urinary_dihydrocaffeic_acid, urinary_vanillic_acid, urinary_hippuric_acid, urinary_Homovanillic_acid, urinary_3PHBA_4sul_4PHBA_3sul, total_serum_derived_metabolites, total_urine_derived_metabolites, framingham_CV_score, CRP_ug_ml_, IL6_pg_ml_, TNF__pg_ml_"
#args.list$exclude_subjects <-"671, 1025, 1074, 1098, 1131, 1170, 1355, 1392, 1479, 1513, 1704"
#args.list$exclude_samples <- "VOL1151_0M, VOL1151_6M, VOL1280_0M, VOL1280_6M, VOL1384_6M, VOL631_0M, VOL631_6M, VOL1267_0M, VOL1267_6M, VOL1389_0M, VOL1389_6M"
#args.list$metadata_categorical ="age,sex, group, timepoint,metabolicsyndromecriteria,statin,bp_med"


maaslin.args <- args.list
maaslin.args$metadata_exclude <- NULL
maaslin.args$subjects_exclude <- NULL
maaslin.args$samples_exclude <- NULL
maaslin.args$metadata_categorical <- NULL
maaslin.args$abundance <- NULL
maaslin.args$metadata <- NULL
maaslin.args$input_dir <- NULL
maaslin.args$output_dir <- NULL
maaslin.args$maaslin_output <- NULL
maaslin.args$normalize <- NULL


args_string <- c()
if (length(maaslin.args)>0){
 for( arg in names(maaslin.args) ) {
   args_string<- c(args_string,paste0(arg,"='",maaslin.args[[arg]],"'"))
 }
}
print(args_string)

# Set input dir as working
setwd(args.list$input_dir)

# Read metadata
metadata <- read.delim(args.list$metadata, header = TRUE, sep="\t", fileEncoding="Latin1")

# Exclude extra covariates
if(!is.null(args.list$metadata_exclude)){
  exclude <- as.list(strsplit(gsub(" ","",args.list$metadata_exclude), ",")[[1]])
  metadata_selected <- metadata[,-which(names(metadata) %in% exclude)]
  metadata <- metadata_selected
}
if(!is.null(args.list$subject_exclude)){
  subjects_exclude <- as.list(gsub(" ","",strsplit(args.list$subjects_exclude), ",")[[1]])
  metadata_selected <- metadata[,-which(names(metadata) %in% subjects_exclude)]
  metadata <- metadata_selected
}

metadata_cc <- metadata[complete.cases(metadata), ]
metadata <- metadata_cc

# Make sure metadata categorical is char
if(!is.null(args.list$metadata_categorical)){
  categorical <- as.list(strsplit(gsub(" ","",args.list$metadata_categorical), ",")[[1]])
  for(i in 1:length(categorical)){
    covariate <- as.character(categorical[i])
    metadata[,covariate] <- as.character(metadata[,covariate])
   }
 }

# Set rownames
rownames(metadata) <- metadata[,1]
metadata <- metadata[,-1]

# Read taxonomic profile
tax_profile <- read.delim(args.list$abundance, header = TRUE, sep="\t")

# Exclude some samples
if(!is.null(args.list$samples_exclude)){
  samples_exclude <- as.list(strsplit(gsub(" ","",args.list$samples_exclude), ",")[[1]])
  tax_profile_selected <- tax_profile[,-which(names(tax_profile) %in% samples_exclude)]
  tax_profile <- tax_profile_selected
}
# Set row names
rownames(tax_profile) <- tax_profile[,1]
tax_profile <- tax_profile[,-1]

# Unstratify if taxa or path abundance stratified
if (grepl("_taxonomic_profile",colnames(tax_profile)[2])){
  colnames(tax_profile) <- gsub("_taxonomic_profile","",colnames(tax_profile))
  subset_ind <- grep( "s__|^unclassified", rownames(tax_profile), invert = F )
  unstrat_taxa_t <- tax_profile[subset_ind,]
  subset_ind <- grep( "t__", rownames(unstrat_taxa_t), invert = T )
  unstrat_taxa <- unstrat_taxa_t[subset_ind,]
  
}else if(grepl("_Abundance",colnames(tax_profile)[2])){
  colnames(tax_profile) <- gsub("_Abundance","",colnames(tax_profile)) 
  subset_ind <- grep( "\\|", rownames(tax_profile), invert = T )
  unstrat_taxa <- tax_profile[subset_ind,]
  
}else{
  unstrat_taxa <- tax_profile
}

if(args.list$normalize == "yes"){
# Normalize
sz=dim(unstrat_taxa)
csums<-colSums(unstrat_taxa[2:sz[1],2:sz[2]])
taxa_profile <- as.data.frame(t(unstrat_taxa[2:sz[1],2:sz[2]])/csums)
}
merged_taxa_data <- merge(metadata,taxa_profile, by="row.names", sort=TRUE)

abundance_file <- paste0(args.list$maaslin_output,"/","abundance.tsv")
metadata_file <- paste0(args.list$maaslin_output,"/","metadata_processed.tsv")

if(!dir.exists(args.list$maaslin_output)) dir.create(args.list$maaslin_output)

write.table(taxa_profile, abundance_file , sep = "\t", eol = "\n", quote = F, col.names = NA)
write.table(metadata, metadata_file , sep = "\t", eol = "\n", quote = F, col.names = NA)

args_all <- c(list(abundance_file, metadata_file, args.list$maaslin_output),args_string)
if (length(args_string)>0){
  Maaslin2(taxa_profile, metadata, args.list$maaslin_output)
  #do.call(Maaslin2,c(taxa_profile, metadata, args.list$maaslin_output, as.list(args_string)))
  #do.call('Maaslin2',args_all)
}else{
  Maaslin2(taxa_profile, metadata, args.list$maaslin_output)
  #Maaslin2(abundance_file, metadata_file, args.list$maaslin_output,  min_abundance=0.001,cores=4)
}