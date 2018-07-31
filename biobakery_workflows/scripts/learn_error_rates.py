#!/usr/bin/env python

from rpy2 import robjects
from rpy2.robjects.packages import importr
import os.path
import sys
import re
import argparse

base = importr('base')
dada2 = importr('dada2')
ggplot2 = importr('ggplot2')

def parse_arguments(args):
    """ Parse the arguments"""

    parser = argparse.ArgumentParser(
        description="Learn error rates for each sample\n",
        formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument(
        "--filt_path",
        help="path to filtered files \n[REQUIRED]",
    )
    parser.add_argument(
        "--error_ratesF_png",
        help="error rates forward reads plot\n[REQUIRED]",
    )
    parser.add_argument(
        "--error_ratesR_png",
        help="error rates reverse reads plot\n[REQUIRED]",
    )
    parser.add_argument(
        "--error_ratesF_rds",
        help="error rates forward reads data\n[REQUIRED]",
    )
    parser.add_argument(
        "--error_ratesR_rds",
        help = "error rates reverse reads data\n[REQUIRED]",
    )
    parser.add_argument(
        "--threads",
        help="number of threads\n[REQUIRED]",
    )

    return parser.parse_args()


def main():

    # parse arguments from the user
    args = parse_arguments(sys.argv)

    filtFs = sorted(
        [os.path.join(args.filt_path, f) for f in os.listdir(args.filt_path) if re.match(r'(.*)_F_filt.fastq(.*)', f)])
    filtRs = sorted(
        [os.path.join(args.filt_path, f) for f in os.listdir(args.filt_path) if re.match(r'(.*)_R_filt.fastq(.*)', f)])

    robjects.globalenv["filtFs"] = filtFs
    robjects.globalenv["filtRs"] = filtRs
    robjects.globalenv["error_ratesF_png"] = args.error_ratesF_png
    robjects.globalenv["error_ratesR_png"] = args.error_ratesR_png
    robjects.globalenv["error_ratesF_rds"] = args.error_ratesF_rds
    robjects.globalenv["error_ratesR_rds"] = args.error_ratesR_rds
    robjects.globalenv["threads"] = args.threads


    robjects.r('''

        filtF <- unlist(filtFs, use.names=FALSE)
        filtR <- unlist(filtRs, use.names=FALSE)

        set.seed(100)
        # Filtered forward reads
        errF <- dada2::learnErrors(filtF, nread=1e3, multithread=as.numeric(threads))
        # Filtered reverse reads
        errR <- dada2::learnErrors(filtR, nread=1e3, multithread=as.numeric(threads))

        # Visualize the estimated error rates
         ggplot2::ggsave(error_ratesF_png, dada2::plotErrors(errF, nominalQ=TRUE), device = "png")
         ggplot2::ggsave(error_ratesR_png, dada2::plotErrors(errR, nominalQ=TRUE), device = "png")
         
         saveRDS(errF, error_ratesF_rds)
         saveRDS(errR, error_ratesR_rds)

        ''')

    errF = robjects.globalenv["errF"]
    errR = robjects.globalenv["errR"]

if __name__ == "__main__":
        main()
