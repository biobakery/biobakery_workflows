#!/usr/bin/env python

from rpy2 import robjects
from rpy2.robjects.packages import importr
import os.path
import sys
import re
import argparse

base = importr('base')
dada2 = importr('dada2')
gridExtra = importr('gridExtra')
grdevices = importr('grDevices')
ggplot2 = importr('ggplot2')


def parse_arguments(args):
    """ Parse the arguments"""

    parser = argparse.ArgumentParser(
        description="Filter and trim arguments\n",
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        "--input",
        help="input directory where fastq(.gz) files are \n[REQUIRED]",
    )
    parser.add_argument(
        "--output",
        help="output directory \n[REQUIRED]",
    )
    parser.add_argument(
        "--maxee",
        help="maxee \n[REQUIRED]",
    )
    parser.add_argument(
        "--trunc_len_max",
        help="truncate at max length\n[REQUIRED]",
    )
    parser.add_argument(
        "--readcounts_tsv_path",
        help="file path to write read counts tsv format\n[REQUIRED]",
    )
    parser.add_argument(
        "--readcounts_rds_path",
        help = "file path to write read counts rds file\n[REQUIRED]",
    )
    parser.add_argument(
        "--plotF",
        help="file to plot forward reads quality\n[REQUIRED]",
    )
    parser.add_argument(
        "--plotR",
        help="file to plot referse reads quality\n[REQUIRED]",
    )
    parser.add_argument(
        "--pair_id",
        help="pair_identifier\n[REQUIRED]",
    )
    parser.add_argument(
        "--threads",
        help="number of threads\n[REQUIRED]",
    )

    return parser.parse_args()


def main():

    # parse arguments from the user
    args = parse_arguments(sys.argv)

    pair_id0 = args.pair_id.split("1_")[0]
    pair_id1 = pair_id0 + "1"
    pair_id2 = pair_id0 + "2"

    fnFs = sorted([f for f in os.listdir(args.input) if re.match(r'(.*)'+ pair_id1 + '(.*)', f)])
    fnRs = sorted([f for f in os.listdir(args.input) if re.match(r'(.*)'+ pair_id2 + '(.*)', f)])

    sample_names = [f.split(pair_id1)[0] for f in fnFs]
    sample_namesR = [f.split(pair_id2)[0] for f in fnRs]

    print(sample_names)

    if (sample_names != sample_namesR):
        sys.exit("ERROR: Forward and reverse files do not match.")

    fnFs = [os.path.join(args.input, f) for f in fnFs]
    fnRs = [os.path.join(args.input, f) for f in fnRs]

    filt_path = os.path.join(args.output, "filtered_input")

    # Create filtered_input/ subdirectory for storing filtered fastq reads
    if not os.path.isdir(filt_path):
        try:
            os.makedirs(filt_path)
        except EnvironmentError:
            sys.exit("ERROR: Unable to create filtered_input directory: ")

    fileext = os.path.splitext(fnFs[1])[-1]
    if fileext == ".gz": fileext = ".fastq.gz"
    print(fileext)

    filtFs = [os.path.join(args.output, filt_path, name + "_F_filt" + fileext) for name in sample_names]
    filtRs = [os.path.join(args.output, filt_path, name + "_R_filt" + fileext) for name in sample_names]

    robjects.globalenv["output"] = args.output
    robjects.globalenv["fnFs"] = fnFs
    robjects.globalenv["fnRs"] = fnRs
    robjects.globalenv["maxee"] = args.maxee
    robjects.globalenv["trunc_len_max"] = args.trunc_len_max
    robjects.globalenv["readcounts_tsv_path"] = args.readcounts_tsv_path
    robjects.globalenv["readcounts_rds_path"] = args.readcounts_rds_path
    robjects.globalenv["reads_plotF"] = args.plotF
    robjects.globalenv["reads_plotR"] = args.plotR
    robjects.globalenv["ext"] = fileext
    robjects.globalenv["filtFs"] = filtFs
    robjects.globalenv["filtRs"] = filtRs
    robjects.globalenv["threads"] = args.threads

    robjects.r('''

    output.path <- normalizePath(output)

    fnF <- unlist(fnFs, use.names=FALSE)
    fnR <- unlist(fnRs, use.names=FALSE)
    filtF <- unlist(filtFs, use.names=FALSE)
    filtR <- unlist(filtRs, use.names=FALSE)

    # Generate plots and save to file
    # Forward reads
    fwd.qc.plots.list <- list()
    for( i in 1 : length(fnF)) {
      fwd.qc.plots.list[[i]] <- dada2::plotQualityProfile(fnF[i])
      rm(i)
    }
    # Save to file
    png(reads_plotF)
    gridExtra::marrangeGrob( fwd.qc.plots.list, ncol=2, nrow=3, top = NULL )
    dev.off()
    rm(fwd.qc.plots.list)

    # Reverse reads
    rev.qc.plots.list <- list()
    for( i in 1 : length(fnR)) {
      rev.qc.plots.list[[i]] <- dada2::plotQualityProfile(fnR[i])
      rm(i)
    }
    # Save to file
    png(reads_plotR)
    gridExtra::marrangeGrob( rev.qc.plots.list, ncol=2, nrow=3, top = NULL )
    dev.off()
    rm(rev.qc.plots.list)


    # Filter the forward and reverse reads:
    # Note that:
    # 1. Reads are both truncated and then filtered using the maxEE expected errors algorighm from UPARSE.
    # 2. Reverse reads are truncated to shorter lengths than forward since they are much lower quality.
    # 3. _Both_ reads must pass for the read pair to be output.
    # 4. Output files are compressed by default.
    trunc_len_max2 <- strtoi(trunc_len_max)
    trunc_len_max1 <- trunc_len_max2 + 40
    maxee1 <- strtoi(maxee)
    maxee2 <- maxee1 * 2

    rd.counts <- as.data.frame(
      filterAndTrim(fnF, filtF, fnR, filtR, truncLen=c(trunc_len_max1,trunc_len_max2),
                    maxN=0, maxEE=c(maxee1,maxee2), truncQ=2, rm.phix=TRUE,
                    compress=TRUE, multithread=as.numeric(threads)) 
    )
    
    # Table of before/after read counts
    rd.counts$ratio <- round( rd.counts$reads.out / rd.counts$reads.in, digits = 2 )
    
    # Write rd.counts table to file in output folder
    saveRDS(rd.counts, readcounts_rds_path )
    write.table( rd.counts, readcounts_tsv_path, sep = "\t", quote = F, eol = "\n", col.names = NA )
    
    ''')

    rd_counts = robjects.globalenv["rd.counts"]

if __name__ == "__main__":
        main()