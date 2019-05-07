#!/usr/bin/env Rscript

library(phyloseq)
library(vegan)
library(RColorBrewer)
library(pheatmap)
library(grDevices)


###########################################################################################################
#####  Implements permanova analysis
#####  Requires:
#####   -metadata.tsv file (sample ids column should be named 'sample')
#####   -closed reference either biom or tsv file 
#####   -closed_reference.tre fasttree file (optional, if available)
#####  Arguments to script are:
#####   -input_dir  (folder where required files are)
#####   -normalize  (if 'no' it will skip relative abundance calculation, by default it calculates)
#####   -fields  (if given, it will create plots with those specific fields in addition to general plot with all fields,
#####    example: source*time or diet+sex) 
#####   -workflow (wmgx or 16s-default)
#####  It outputs:
#####   -permanova_heatmaps.pdf with p-values and r square plots for all fields, and additional plots, 
#####    if specific 'fields' argument is given
#####   -adonis-output.txt file with information about each field
#############################################################################################################

## Collect arguments
args <- commandArgs(TRUE)

## Parse arguments (we expect the form --arg=value)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
args.df <- as.data.frame(do.call("rbind", parseArgs(args)))

args.list <- as.list(as.character(args.df$V2))
names(args.list) <- args.df$V1

#args.list$input_dir <- "/Users/anamailyan/dada_demo_output_vsearch"
#args.list$output_dir <- "/Users/anamailyan/dada_demo_output_vsearch_stats"
#rgs.list$adonis_dir <- "/Users/anamailyan/dada_demo_output_vsearch_stats/permanova_path"
#args.list$metadata <- "/Users/anamailyan/dada_demo_output_vsearch/metadata.tsv"
#args.list$abundance <- "/Users/anamailyan/dada_demo_output_vsearch/humann2/merged/pathabundance_relab.tsv"
args.list$normalize <- "yes"
#args.list$workflow <- "wmgx"
#args.list$workflow <- "16s"
#args.list$fields <- "timepoint*bmi"

## Arg help


## Arg1 default
if(is.null(args.list$input_dir)) {
  #stop("At least one argument must be supplied (input directory).\n", call.=FALSE)
}

if(is.null(args.list$workflow)){args.list$workflow <- "16s"}

# Print args list to STDOUT
for( i in names(args.list) ) {
  cat( i, "\t", args.list[[i]], "\n")
}

if(!dir.exists(args.list$output_dir)) dir.create(args.list$output_dir)
if(!dir.exists(args.list$adonis_dir)) dir.create(args.list$adonis_dir)

# These variables are passed to the workflow
input.path <- normalizePath( args.list$input_dir )
output.path <- normalizePath( args.list$output_dir )
adonis.path <- normalizePath( args.list$adonis_dir )
setwd(args.list$input_dir)

metadata_file = args.list$metadata
if(file.exists(metadata_file)){
    meta.dat <- read.table(metadata_file,header = TRUE)
    meta.dat <- meta.dat[complete.cases(meta.dat),]
}else{
    stop("Exiting. Metadata file  not found", call.=FALSE)
}

closed_ref_biom = paste0(args.list$input_dir, "/all_samples_taxonomy_closed_reference.biom")
closed_ref_tsv = args.list$abundance
if(file.exists(closed_ref_biom)){
      otu.dat <-  phyloseq::import_biom( 
      BIOMfilename = closed_ref_biom,
      parseFunction = parse_taxonomy_greengenes
      )
      D <- as.data.frame(otu_table(otu.dat))
}else if(file.exists(closed_ref_tsv)){
        otu.dat <- read.delim(closed_ref_tsv, sep = '\t', header=TRUE )
     
       if(grepl("_taxonomic_profile",colnames(otu.dat)[2])){
         colnames(otu.dat) <- gsub("_taxonomic_profile","",colnames(otu.dat))
         rownames(otu.dat) <- otu.dat[,1]
         subset_ind <- grep( "s__|unclassified", rownames(otu.dat), invert = F )
         unstrat_otu_st.dat <- otu.dat[subset_ind,]
         subset_ind <- grep( "t__", rownames(unstrat_otu_st.dat), invert = T )
         otu.dat <- unstrat_otu_st.dat[subset_ind,]
         
         }else if(grepl("_Abundance",colnames(otu.dat)[2])){
           colnames(otu.dat) <- gsub("_Abundance","",colnames(otu.dat)) 
           rownames(otu.dat) <- otu.dat[,1]
           subset_ind <- grep( "\\|", rownames(otu.dat), invert = T )
           otu.dat <- otu.dat[subset_ind,]
           
         }else{
           colnames(otu.dat) <- gsub( "^X", "", colnames(otu.dat))
         }
      
     otu.dat <- otu.dat[2:length(otu.dat)] 
     otu.dat.filt <- otu.dat[,which(colnames(otu.dat) %in% meta.dat$sample)]
     D <- otu.dat.filt
}else{
    stop("Exiting. Closed reference file not found in ", call.=FALSE)
}

phylotree_file = paste0(args.list$input_dir, "/closed_reference.tre")
if (file.exists(phylotree_file)){
    phylotree.data <- phyloseq::read_tree(phylotree_file)
    meta <- phyloseq::sample_data(meta.dat)
    rownames(meta) <- phyloseq::sample_names(otu.dat)
    otu.dat <- phyloseq::merge_phyloseq(otu.dat, meta, phylotree.data)
    D <- as.data.frame(otu_table(otu.dat))
}

D <- D[, which(colSums(D)>0)]

if(is.null(args.list$normalize)){args.list$normalize <- "yes"}
if (args.list$normalize == "yes"){
    D <- D/colSums(D)
    D.filt <- D[ apply(D, 1, function(x) sum(x > 0.0001) > 0.1 * ncol(D) ), ]
}else{
    D.filt <- D
}
 
D.filt <- D.filt[, which(colSums(D.filt)>0)]
meta.dat.filt <- meta.dat[meta.dat$sample %in% colnames(D.filt) ,,drop=T]
meta.dat <- meta.dat.filt
 
#write.table(D, "D.txt", sep = "\t", eol = "\n", quote = F, row.names=TRUE,col.names = TRUE)
#write.table(D.filt, "D-filt.txt", sep = "\t", eol = "\n", quote = F, row.names=TRUE, col.names = TRUE)

D.dist <- vegan::vegdist(t(D.filt), "bray")

sink(file = paste0(adonis.path,"/adonis-output.txt"), append = TRUE, type = "output", split = FALSE)
for (metacol in colnames(meta.dat)[2:length(colnames(meta.dat))]){
  adonis.col <- vegan::adonis(as.formula(paste("D.dist ~ ", metacol)), data = meta.dat[,-1]) 
  print(adonis.col)
}
sink()

# hmap.plots <-list()
plots.list <- list()
# Figure 1. The PERMANOVA
adonis.full <- vegan::adonis(D.dist ~. , data = meta.dat[,-1])
prflength <- length(adonis.full$aov.tab$`Pr(>F)`)-2

pdf(file=paste0(adonis.path,"/permanova_heatmaps.pdf"))

# Figure 1(A). P-value Summary with 999 permutations
# hmap.plots[[1]] <- "hmap-pv.pdf"
pheatmap::pheatmap(adonis.full$aov.tab$`Pr(>F)`[1:prflength], display_numbers = TRUE, number_color= "black",cluster_rows = FALSE,
         cluster_cols = FALSE, number_format = "%.3f", 
         color = colorRampPalette(brewer.pal(n = 10, name = "RdYlGn"))(501),fontsize = 12, fontsize_number = 12,
         breaks=pbeta(seq(0, 1, len=501), 1.1, 0.15), main = "P-value Summary", labels_row = row.names(adonis.full$aov.tab)[1:prflength])

plots.list[[1]] <- recordPlot(load="pheatmap,RColorBrewer", attach="pheatmap,RColorBrewer")
#graphics.off()

# Figure 1(B). R-square Summary 
# hmap.plots[[2]] <- "hmap-r2.pdf"
#png(file=paste0(adonis.path,"/all_r2.png"))
pheatmap::pheatmap(100*adonis.full$aov.tab$R2[1:prflength], display_numbers = TRUE, number_color= "black",cluster_rows = FALSE,
         cluster_cols = FALSE, number_format = "%.3f%%",
         color = colorRampPalette(brewer.pal(n = 9, name = "Oranges"))(501),fontsize = 12, fontsize_number = 12,
         breaks=100*pbeta(seq(0, 1, len=501), 1.1, 0.1), main = "R-square Summary", labels_row = row.names(adonis.full$aov.tab)[1:prflength])

plots.list[[2]] <- recordPlot(load="pheatmap,RColorBrewer", attach="pheatmap,RColorBrewer")
#graphics.off()

if (!is.null(args.list$fields)){
   
    if (grepl("\\*", args.list$fields)){
      splitchar <- "\\*"
    }else if(grepl("\\+", args.list$fields)){
      splitchar <- "\\+"
    }else{
      stop("At least 2 fields must be given with an operator * or + ")
    }
    
    meta_cols <- strsplit(args.list$fields, splitchar)
    for (mcol in meta_cols){
      if (ncol(meta.dat[mcol %in% colnames(meta.dat)]) == 0){
        stop("Field not found in metadata")
      }
    }
   
    fmla <- as.formula(paste0("D.dist ~",args.list$fields ))
    adonis.full <- adonis(fmla, data = meta.dat[,-1])

    prflength <- length(adonis.full$aov.tab$`Pr(>F)`)-2
    png(file=paste0(adonis.path,"/fields_heatmap.png"))
   # hmap.plots[[3]] <- "hmap-pv_fields.pdf"
    pheatmap::pheatmap(adonis.full$aov.tab$`Pr(>F)`[1:prflength], display_numbers = TRUE, number_color= "black",cluster_rows = FALSE,
           cluster_cols = FALSE,  number_format = "%.3f",
           color = colorRampPalette(brewer.pal(n = 10, name = "RdYlGn"))(501),fontsize = 12, fontsize_number = 12,
           breaks=pbeta(seq(0, 1, len=501), 1.1, 0.15), main = "P-value Summary", 
           labels_row = row.names(adonis.full$aov.tab)[1:prflength])
  
    plots.list[[3]] <- recordPlot(load="pheatmap,RColorBrewer", attach="pheatmap,RColorBrewer")
#    graphics.off()
   
#    png(file=paste0(adonis.path,"/fields_r2.png")) 
   #  hmap.plots[[4]] <- "hmap-r2_fields.pdf"
     pheatmap::pheatmap(100*adonis.full$aov.tab$R2[1:prflength], display_numbers = TRUE, number_color= "black",cluster_rows = FALSE,
           cluster_cols = FALSE,  number_format = "%.3f%%",
           color = colorRampPalette(brewer.pal(n = 9, name = "Oranges"))(501),fontsize = 12, fontsize_number = 12,
           breaks=100*pbeta(seq(0, 1, len=501), 1.1, 0.1), main = "R-square Summary", 
           labels_row = row.names(adonis.full$aov.tab)[1:prflength])
     plots.list[[4]] <- recordPlot(load="pheatmap,RColorBrewer", attach="pheatmap,RColorBrewer")
#     graphics.off()
  
}

for (plot in plots.list) {
  replayPlot(plot, reloadPkgs=TRUE)
 }
 graphics.off()


