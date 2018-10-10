#!/usr/bin/env python

# Given a raw interleaved fastq file and a balanced interleaved fastq (containing no orpahns)
# generate the read ID lists necessary to extract all orphans from the original file.
#
# Makes use of the seqtk utility

import argparse
import os
import subprocess


def parse_cli_arguments():
    """
    """
    parser = argparse.ArgumentParser('Extracts orphan reads from a interleaved sequence file '
                                     'and produces an orphan sequence file.')
    parser.add_argument('-r', '--raw-sequence', required=True, 
                        help='The raw interleaved sequence file.')
    parser.add_argument('-b', '--balanced-sequence', required=True,
                        help='Balanced sequence file with no orphan sequences.')
    parser.add_argument('-o', '--output-dir', required=True,
                        help='Output directory to write orphan sequence files too.')    

    return parser.parse_args()


def get_ids_from_sequences(sample_name, raw_seq, balanced_seq, out_dir):
    """
    Extracts raw, balanced and orphan sequence IDs from the provided 
    sequence files.
    """
    raw_ids = os.path.join(out_dir, "%s.raw_ids.txt" % sample_name)
    balanced_ids = os.path.join(out_dir, "%s.matched_ids.txt" % sample_name)
    orphan_ids = os.path.join(out_dir, "%s.orphan_ids.txt" % sample_name)

    for (input_seq, output_ids) in [(raw_seq, raw_ids), (balanced_seq, balanced_ids)]:
        with open(output_ids, 'wb') as out_ids:
            ps_grep = subprocess.Popen(['grep', '-e', '^@.*/[1|2]$', input_seq], stdout=subprocess.PIPE)
            ps_sed = subprocess.Popen(['sed', '-e', 's/^@//'], stdin=ps_grep.stdout, stdout=subprocess.PIPE)
            ps_grep.stdout.close()

            ps_sort = subprocess.Popen(['sort'], stdin=ps_sed.stdout, stdout=out_ids)
            ps_sed.stdout.close()

            ps_sort.communicate()

    with open(orphan_ids, 'wb') as orphan_ids_out:
        p = subprocess.Popen(['comm', '-23', raw_ids, balanced_ids], stdout=orphan_ids_out)
        p.communicate()

    return (raw_ids, balanced_ids, orphan_ids)


def generate_orphan_sequences(sample_name, raw_seqs, orphan_ids, out_dir):
    """
    Generates an orphan sequence file from the supplied interleaved sequence file.
    """
    orphan_seqs_file = os.path.join(out_dir, "%s_orphans.fastq" % sample_name)

    with open(orphan_seqs_file, 'wb') as orphan_seqs:
        p = subprocess.Popen(['seqtk', 'subseq', raw_seqs, orphan_ids], stdout=orphan_seqs)
        p.communicate()


def main(args):
    sample_name = os.path.basename(args.raw_sequence).split(os.extsep, 1)[0]

    (raw_ids, balanced_ids, orphan_ids) = get_ids_from_sequences(sample_name,
                                                                 args.raw_sequence, 
                                                                 args.balanced_sequence,
                                                                 args.output_dir)

    generate_orphan_sequences(sample_name, args.raw_sequence, orphan_ids, args.output_dir)

    os.remove(raw_ids)
    os.remove(balanced_ids)
    os.remove(orphan_ids)


if __name__ == "__main__":
    main(parse_cli_arguments())
