#!/usr/bin/env Rscript

library(docopt)
library(dplyr)
library(data.table)
library(doParallel)
library(tibble)
library(tidyr)
library(stringr)

'Usage:
   mash_clusters.R [--mash <mash> --checkm <checkm> --phylo <phylophlan_relab> --mag_dir <mag_dir> --out_dir <out_dir> --threads <threads>]

Options:
   --mash Mash table
   --checkm CheckM QA table
   --phylo Phylophlan Relabeled table
   --mag_dir MAGs directory
   --out_dir Output directory
   --threads Number of threads to use [default: 1]
   
' -> doc 

opts <- docopt(doc)

out_dir <- paste0(gsub("/$", "", opts$out_dir), "/")

mag_dir <- paste0(gsub("/$", "", opts$mag_dir), "/")

threads <- as.numeric(opts$threads)
registerDoParallel(cores=threads)

############################################
# function to tidy the Mash distance table #
############################################

boil_the_spuds <- function(input, genomes) {
  
  mash <- input %>%
    rename(genome = 1) %>%
    mutate(genome = gsub(".*\\/", "", genome)) %>%
    filter(genome %in% genomes) %>%
    setNames(gsub(".*\\/", "", names(.))) %>%
    select(genome, all_of(genomes)) %>%
    arrange(genome) %>%
    column_to_rownames("genome") %>%
    select(sort(tidyselect::peek_vars())) %>%
    as.dist()
  
}

#################################
# function to identify clusters #
#################################

get_clusters <- function(dist_matrix, cutoff) {
  
  hclust_out <- hclust(dist_matrix, "average")
  
  hclust_out$height <- round(hclust_out$height, 6)
  
  # cut tree
  as.data.frame(cutree(hclust_out, h=cutoff)) %>%
    rownames_to_column("genome") %>%
    rename(cluster=2) %>%
    mutate(cluster = paste0("cluster_", cluster)) 
  
}

#################################################
# function to pick representative from clusters #
#################################################

get_reps <- function(clusters) {
  
  set.seed(123)
  
  clusters %>%
    inner_join(., checkm, by=c("genome"="bin_id")) %>%
    group_by(cluster) %>%
    # rank completeness (decreasing)
    arrange(desc(completeness)) %>%
    mutate(rank_completeness = row_number()) %>%
    # rank contamination (increasing)
    arrange(contamination) %>%
    mutate(rank_contamination = row_number()) %>%
    # rank n50 (decreasing)
    arrange(desc(n50)) %>%
    mutate(rank_n50 = row_number()) %>%
    # primary rank
    mutate(rank_primary = rank_completeness + rank_contamination) %>%
    top_n(-5, rank_primary) %>%
    # secondary rank
    mutate(rank_secondary = rank_contamination) %>%
    top_n(-1, rank_secondary) %>%
    # pick 1 genome if multiple genomes rank equally
    sample_n(1) %>%
    ungroup() %>%
    select(cluster, genome)
  
}

####################################
# function to list cluster members #
####################################

get_members <- function(x) {
  
  x %>%
    group_by(cluster) %>%
    summarise(cluster_members = paste0(genome, collapse=",")) %>%
    ungroup()
  
}

###################
# CheckM QA table #
###################

phylophlan_relab <- fread(opts$phylo)

checkm <- fread(opts$checkm) %>%
  filter(bin_id %in% phylophlan_relab$mag[phylophlan_relab$taxon == "UNKNOWN"]) %>%
  filter(keep == "keep") %>%
  mutate(bin_id = paste0(bin_id, ".fa"))

#######################
# Mash distance table #
#######################

mash_raw <- fread(opts$mash)

if (nrow(mash_raw)[1] < 2) {
  colnames(mash_raw)[2] <- gsub(".*\\/", "", colnames(mash_raw))[2]
  sgb_dir <- paste0(out_dir, "sgbs/")
  dir.create(sgb_dir, recursive = TRUE)
  
  sgbs <- data.frame(cbind("cluster"="cluster_1", "genome"=colnames(mash_raw)[2], "cluster_members" = colnames(mash_raw)[2], n_gemones=1), check.names = F)
  sgbs <- left_join(sgbs, checkm, by=c("genome"="bin_id"))
  sgbs$cluster_name <- "cluster_1"
  sgbs$sgb <- "sgb_01"
  
  write.table(sgbs, paste0(sgb_dir, "SGB_info.tsv"), quote=FALSE, sep="\t", row.names=FALSE)
  
  # remove pre-existing bins from the output directory
  files <- list.files(sgb_dir)
  if (any(grepl("\\.fa", files))) {
    system(paste0("rm ", sgb_dir, "*.fa"))
  }
  
  commands <- sgbs %>%
    select(genome, sgb) %>%
    mutate(command = paste0("cp ", mag_dir, genome, " ", sgb_dir, sgb, ".fa")) %>%
    select(command) %>%
    .$command
  
  foreach(i=commands) %dopar% {
    system(i)
  }
  
  quit()
}

# list of MAGs to include
mags <- levels(as.factor(checkm$bin_id))

# make distance matrix
mash <- boil_the_spuds(input=mash_raw, genomes=mags)

# cluster + cut at 5%
mash_clusters <- get_clusters(dist_matrix=mash, cutoff=0.05)

# map
cluster_members <- get_members(mash_clusters)

######################################################
# merge with CheckM and pick cluster representatives #
######################################################

results <- get_reps(clusters=mash_clusters) %>%
  inner_join(., cluster_members, by="cluster")

mash_res <- results %>%
  select(genome, cluster_members)

#######################
# run fastANI on sgbs #
#######################

ani_dir <- paste0(out_dir, "fastANI/")
dir.create(ani_dir, recursive=T)

# fastANI input

ani_in <- paste0(ani_dir, "SGB_list.txt")

ani_list <- results %>%
  select(genome) %>%
  mutate(genome = paste0(mag_dir, genome))

write.table(ani_list, ani_in, quote=F, col.names=F, row.names=F)

# fastANI output

ani_out <- paste0(ani_dir, "fastANI_sgbs.txt")

# run fastANI

fastANI_command <- paste0(
  "fastANI --ql ",
  ani_in,
  " --rl ",
  ani_in,
  " -o ",
  ani_out,
  " -t ",
  threads,
  " --matrix"
)

system(fastANI_command)

# read fastANI outputs

ani_raw <- fread(ani_out) %>%
  mutate(V3 = 1 - V3 / 100)

# 0.25 chosen based on reporting threshold from https://github.com/ParBLiSS/FastANI
ani_dist <- ani_raw %>%
  select(V1, V2, V3) %>%
  pivot_wider(names_from=V2, values_from=V3, values_fill=0.25) %>%
  arrange(V1) %>%
  column_to_rownames("V1") %>%
  select(sort(tidyselect::peek_vars())) %>%
  as.dist()

ani_clusters <- get_clusters(dist_matrix=ani_dist, cutoff=0.05) %>%
  mutate(genome = gsub(".*\\/", "", genome)) %>%
  inner_join(., mash_res, by="genome") %>%
  select(cluster, cluster_members) %>%
  separate_rows(cluster_members, sep=",") %>%
  rename(genome=2)

ani_results <- get_reps(clusters=ani_clusters)

ani_cluster_members <- get_members(ani_clusters)

ani_res <- ani_results %>%
  inner_join(., ani_cluster_members, by="cluster") %>%
  mutate(n_genomes = str_count(cluster_members, ",") + 1) %>%
  inner_join(., checkm, by=c("genome"="bin_id"))

sgbs <- ani_res %>%
  mutate(cluster_name = paste0("cluster_", row_number())) %>%
  mutate(sgb = sprintf(paste0("sgb_%0", nchar(nrow(.)), "d"), as.numeric(gsub(".*_", "", cluster_name))))

#################################
# copy SGBs to output directory #
#################################

sgb_dir <- paste0(out_dir, "sgbs/")
dir.create(sgb_dir, recursive = TRUE)

write.table(sgbs, paste0(sgb_dir, "SGB_info.tsv"), quote=FALSE, sep="\t", row.names=FALSE)

# remove pre-existing bins from the output directory
files <- list.files(sgb_dir)
if (any(grepl("\\.fa", files))) {
  system(paste0("rm ", sgb_dir, "*.fa"))
}

commands <- sgbs %>%
  select(genome, sgb) %>%
  mutate(command = paste0("cp ", mag_dir, genome, " ", sgb_dir, sgb, ".fa")) %>%
  select(command) %>%
  .$command

foreach(i=commands) %dopar% {
  system(i)
}

#