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
        description="Phylogeny\n",
        formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument(
        "--seqtab_rds",
        help="sequence tabe rds file\n[REQUIRED]",
    )
    parser.add_argument(
        "--msa_fasta_path",
        help="msa fasta file\n[REQUIRED]",
    )
    parser.add_argument(
        "--threads",
        help="number of threads\n[REQUIRED]",
    )

    return parser.parse_args()


def main():

    # parse arguments from the user
    args = parse_arguments(sys.argv)


    robjects.globalenv["msa_fasta_path"] = args.msa_fasta_path
    robjects.globalenv["seqtab_rds"] = args.seqtab_rds

    robjects.r('''

    library(phangorn)
    library(msa)

    seqtab.nochim <- readRDS(seqtab_rds)

    # Get sequences
    seqs <- dada2::getSequences(seqtab.nochim)
    names(seqs) <- seqs
    # This propagates to the tip labels of the tree

    # Multiple seqeuence alignment
    mult <- msa::msa(seqs, method="ClustalOmega", type="dna", order="input")
    # Save msa to file; convert first to phangorn object
    phang.align <- phangorn::as.phyDat(mult, type="DNA", names=dada2::getSequences(seqtab.nochim))
    write.phyDat(phang.align, format='fasta', file=msa_fasta_path)

    detach("package:phangorn", unload=TRUE)
    detach("package:msa", unload=TRUE)

    ''')

if __name__ == "__main__":
        main()