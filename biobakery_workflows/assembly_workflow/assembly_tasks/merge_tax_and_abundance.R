#! usr/bin/env Rscript

library(docopt)
library(dplyr)
library(data.table)
library(stringr)
rm(list = ls())

'Usage:
   merge_tax_and_abundance.R [-i <input> --tax <tax> --qa <combined> --sgbs <sgbs> -o <output>]

Options:
   -i abundance_by_X directory
   --tax phylophlan_relab.tsv
   --qa checkm_qa_and_n50.tsv
   --sgbs SGB_info.tsv
   -o output

' -> doc 

opts <- docopt(doc)

# Read folder of abundances
in_dir <- gsub("/$", "", opts$i)
files <- list.files(path=in_dir, full.names=TRUE)

read_nums <- files[grep("read_num",files)]

abundance_files <- files[grepl("\\.abundance\\.",files) & !grepl("done", files)]

# Combine abundance files
abundances <- fread(abundance_files[1])[,c(1,6)]
colnames(abundances) <- c("ID", (str_split(abundance_files[1], "/")[[1]] %>% tail(1) %>%
                             str_split(pattern = "\\.abundance"))[[1]][1])

for (abundance_file in abundance_files[-1]) {
  abundances$ID = as.character(abundances$ID)
  data <- fread(abundance_file)[,c(1,6)]
  colnames(data) <- c("ID", (str_split(abundance_file, "/")[[1]] %>% tail(1) %>%
                        str_split(pattern = "\\.abundance"))[[1]][1])
  if (nrow(data) > 0) {
    abundances <- full_join(abundances, data, by="ID")
  } else {
    abundances[,(str_split(abundance_file, "/")[[1]] %>% tail(1) %>%
                  str_split(pattern = "\\.abundance"))[[1]][1]] <- NA
  }
}

# Merge in bin taxonomies
map <- fread(gsub("/$", "", opts$tax))
map <- map[,c(2,7)]
colnames(map) <- c("ID", "Taxonomy")
abundance <- left_join(abundances, map, by = "ID")

# Keep only medium and high quality
qa <- fread(gsub("/$", "", opts$qa))
qa <- qa[,c(1,5)]

abundance <- left_join(abundance, qa, by=c("ID"="bin_id"))

sgbs <- fread(gsub("/$", "", opts$sgbs))
sgbs <- sgbs[,c(3,11)]
if (nrow(sgbs) > 0) {
  for (i in 1:nrow(sgbs)) {
    if (str_count(sgbs$cluster_members[i], ",") > 0) {
      samples = str_split(sgbs$cluster_members[i], ",")[[1]]
      sgbs[i,] <- list(samples[1], sgbs$sgb[i])
      sgbs <- rbind(sgbs, cbind(samples[-1], rep(sgbs$sgb[i], length(samples) - 1)), use.names=FALSE)
    }
  }
}
colnames(sgbs) <- c("ID", "SGB")
sgbs$ID <- gsub("\\.fa", "", sgbs$ID)

abundance <- left_join(abundance, sgbs, by="ID")
abundance$Taxonomy <- ifelse(abundance$Taxonomy=="UNKNOWN", abundance$SGB, abundance$Taxonomy)

# Reformat
abundance$Taxonomy <- gsub(pattern=" \\(root\\)", replacement="", abundance$Taxonomy)
abundance <- mutate(abundance, Taxonomy = case_when(
  grepl("unbinned", ID) ~ "UNKNOWN",
  is.na(Taxonomy) ~ "UNKNOWN",
  grepl("reject", keep) ~ "UNKNOWN",
  TRUE ~ Taxonomy
))

abundance$ID <- gsub("\\.bin.*", "", abundance$ID)
abundance <- data.frame(abundance, check.names = F)

drops = c("ID", "keep", "SGB")
profile <- abundance[, !(colnames(abundance) %in% drops)]

profile[, colnames(profile) != "Taxonomy"][is.na(profile[, colnames(profile) != "Taxonomy"])] <- 0

# Combine same taxonomic groups
profile <- aggregate(.~Taxonomy,data=profile,FUN=sum)

# Upcycle taxonomic groups 
for (name in profile$Taxonomy) {
  abun_tmp <- profile[profile$Taxonomy==name,-which(colnames(profile) == "Taxonomy")]
  newTax <- name
  while (grepl("\\|", newTax)) {
    newTax <- gsub("(\\|[^\\|]*)$", "", newTax)
    profile[nrow(profile) + 1,] <- c(newTax, abun_tmp)
  }
}

profile[,-1] <- data.frame(lapply(data.frame(profile[,-1]), as.numeric))

# Combine same taxonomic groups again
total_profile <- aggregate(.~Taxonomy,data=profile,FUN=sum)

# Alphabetize
total_profile <- total_profile[order(total_profile$Taxonomy),]

# Replace NA with 0
total_profile[is.na(total_profile)] <- 0

mapped_props <- c()
for (read_num in read_nums) {
  text <- readChar(read_num, file.info(read_num)$size)
  temp <- strsplit(text, "\n")[[1]] %>% as.numeric()
  mapped <- temp[1]/temp[2]
  sample <- gsub(".*\\/", "", read_num) %>% gsub(pattern="\\.mapped_read_num.*", replacement="")
  mapped_props <- rbind(mapped_props, c(sample, mapped))
}

colnames(mapped_props) <- c("ID", "prop")
mapped_props <- as.data.frame(mapped_props)
mapped_props$prop <- as.numeric(as.character(mapped_props$prop))

for (ID in mapped_props$ID) {
  if (mapped_props$prop[mapped_props$ID==ID] == 0 ) {
    total_profile[,'new'] <- rep(0, nrow(total_profile))
    colnames(total_profile)[colnames(total_profile)=="new"] <- ID
  }
  total_profile[,ID] <- total_profile[,ID] * mapped_props$prop[mapped_props$ID==ID]
  total_profile[total_profile$Taxonomy=="UNKNOWN", ID] <- total_profile[total_profile$Taxonomy=="UNKNOWN", ID] + 
    100 * (1-mapped_props$prop[mapped_props$ID==ID])
}

write.table(total_profile, opts$o, quote=FALSE, sep="\t", row.names=FALSE)

#
