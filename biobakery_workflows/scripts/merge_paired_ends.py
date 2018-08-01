#!/usr/bin/env python

from rpy2 import robjects
from rpy2.robjects.packages import importr
import os.path
import sys
import re
import argparse

base = importr('base')
dada2 = importr('dada2')

def parse_arguments(args):
    """ Parse the arguments"""

    parser = argparse.ArgumentParser(
        description="Merge paired ends arguments\n",
        formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument(
        "--filt_path",
        help="path to filtered files \n[REQUIRED]",
    )
    parser.add_argument(
        "--error_ratesF_rds",
        help="error rates forward reads data\n[REQUIRED]",
    )
    parser.add_argument(
        "--error_ratesR_rds",
        help="error rates reverse reads data\n[REQUIRED]",
    )
    parser.add_argument(
        "--mergers_rds",
        help="merged reads rds file\n[REQUIRED]",
    )
    parser.add_argument(
        "--sample_names_rds",
        help="sample names rds file\n[REQUIRED]",
    )
    parser.add_argument(
        "--threads",
        help="number of threads\n[REQUIRED]",
    )

    return parser.parse_args()


def main():

    # parse arguments from the user
    args = parse_arguments(sys.argv)

    filtFs = sorted([f for f in os.listdir(args.filt_path) if re.match(r'(.*)_F_filt(.*)', f)])
    filtRs = sorted([f for f in os.listdir(args.filt_path) if re.match(r'(.*)_R_filt(.*)', f)])

    sample_names = [f.split("_F_filt")[0] for f in filtFs]
    sample_namesR = [f.split("_R_filt")[0] for f in filtRs]

    print(sample_names)

    if (sample_names != sample_namesR):
        sys.exit("ERROR: Forward and reverse files do not match.")

    filtFs = [os.path.join(args.filt_path, f) for f in filtFs]
    filtRs = [os.path.join(args.filt_path, f) for f in filtRs]

    robjects.globalenv["error_ratesF_rds"] = args.error_ratesF_rds
    robjects.globalenv["error_ratesR_rds"] = args.error_ratesR_rds
    robjects.globalenv["mergers_rds"] = args.mergers_rds
    robjects.globalenv["filtFs"] = filtFs
    robjects.globalenv["filtRs"] = filtRs
    robjects.globalenv["sample.names"] = sample_names
    robjects.globalenv["sample_names_rds"] = args.sample_names_rds
    robjects.globalenv["threads"] = args.threads

    robjects.r('''
      
        errF <- readRDS(error_ratesF_rds)
        errR <- readRDS(error_ratesR_rds)
        
        filtF <- unlist(filtFs, use.names=FALSE)
        filtR <- unlist(filtRs, use.names=FALSE)
        sample.names <- unlist(sample.names, use.names=FALSE)
        saveRDS(sample.names, sample_names_rds)
        mergers <- vector("list", length(sample.names))
        names(mergers) <- sample.names
        names(filtF) <- sample.names
        names(filtR) <- sample.names
        for(sam in sample.names) {
          cat("Processing:", sam, "\n")
          print(filtF[[sam]])
          derepF <- dada2::derepFastq(filtF[[sam]])
          ddF <- dada2::dada(derepF, err=errF, multithread=as.numeric(threads))
          derepR <- dada2::derepFastq(filtR[[sam]])
          ddR <- dada2::dada(derepR, err=errR, multithread=as.numeric(threads))
          merger <- dada2::mergePairs(ddF, derepF, ddR, derepR)
          mergers[[sam]] <- merger
        }
        rm(derepF); rm(derepR)
        
        saveRDS(mergers, mergers_rds)
    
        ''')

if __name__ == "__main__":
        main()
