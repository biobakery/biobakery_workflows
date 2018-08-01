#!/usr/bin/env python

from rpy2 import robjects
from rpy2.robjects.packages import importr

import sys
import argparse

base = importr('base')
dada2 = importr('dada2')


def parse_arguments(args):
    """ Parse the arguments"""

    parser = argparse.ArgumentParser(
        description="Construct sequence table arguments\n",
        formatter_class=argparse.RawTextHelpFormatter)


    parser.add_argument(
        "--readcounts_steps_path",
        help="file path to write read counts on each step\n[REQUIRED]",
    )
    parser.add_argument(
        "--all_counts_tsv",
        help="path to all counts tsv file\n[REQUIRED]",
    )
    parser.add_argument(
        "--mergers_rds",
        help="path to mergers rds file\n[REQUIRED]",
    )
    parser.add_argument(
        "--sample_names_rds",
        help="path to sample_names\n[REQUIRED]",
    )
    parser.add_argument(
        "--seqtab_rds",
        help="path to sample_names\n[REQUIRED]",
    )
    parser.add_argument(
        "--readcounts_rds",
        help="path to read counts rds file\n[REQUIRED]",
    )
    parser.add_argument(
        "--threads",
        help="number of threads\n[REQUIRED]",
    )

    return parser.parse_args()


def main():

    # parse arguments from the user
    args = parse_arguments(sys.argv)

    robjects.globalenv["mergers_rds"] = args.mergers_rds
    robjects.globalenv["sample_names_rds"] = args.sample_names_rds
    robjects.globalenv["all_samples_counts_tsv"] = args.all_counts_tsv
    robjects.globalenv["readcounts_steps_path"] = args.readcounts_steps_path
    robjects.globalenv["seqtab_rds"] = args.seqtab_rds
    robjects.globalenv["readcounts_rds"] = args.readcounts_rds
    robjects.globalenv["threads"] = args.threads

    robjects.r('''
    
        mergers <- readRDS(mergers_rds)
        sample_names <- readRDS(sample_names_rds)
        sample.names <- unlist(sample_names, use.names=FALSE)

        seqtab <- dada2::makeSequenceTable(mergers)
        dim(seqtab)

        # Inspect distribution of sequence lengths
        table(nchar(getSequences(seqtab)))

        # The sequence table is a matrix with rows corresponding to (and named by) the samples,
        # and columns corresponding to (and named by) the sequence variants.

        # Remove chimeric sequences:
        seqtab.nochim <- dada2::removeBimeraDenovo(seqtab, method="consensus", multithread=as.numeric(threads), verbose=TRUE)
        dim(seqtab.nochim)

        # ratio of chimeric sequence reads
        1 - sum(seqtab.nochim)/sum(seqtab)

        # write sequence variants count table to file
        write.table(t(seqtab.nochim), all_samples_counts_tsv, sep="\t", eol="\n", quote=F, col.names = NA )
        saveRDS(seqtab.nochim, seqtab_rds)
        
        rd.counts <- readRDS(readcounts_rds)
        
        # remove rows with 0 reads after filtering and trimming
        rdf.counts <- rd.counts[which(rd.counts$reads.out != 0), ]

        getN <- function(x)
        sum(getUniques(x))
        track <- cbind(rdf.counts, sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
        colnames(track) <- c("input", "filtered", "ratio", "merged", "tabled", "nonchim")
        rownames(track) <- sample.names
        # print table
        track
        # save to file
        write.table(track, readcounts_steps_path, sep = "\t", quote = F, eol = "\n", col.names = NA )

        ''')


if __name__ == "__main__":
        main()
