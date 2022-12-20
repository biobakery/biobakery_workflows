#!/usr/bin/env python3

###############################################################################
#
# checkm - main program entry point. See checkm/main.py for internals.
#
###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

__author__ = "Donovan Parks, Connor Skennerton, Michael Imelfort"
__copyright__ = "Copyright 2014"
__credits__ = ["Donovan Parks", "Connor Skennerton", "Michael Imelfort"]
__license__ = "GPL3"
__maintainer__ = "Donovan Parks"
__email__ = "donovan.parks@gmail.com"
__status__ = "Development"

import os
import sys
import logging
import tempfile
import argparse

from checkm import main
from checkm.defaultValues import DefaultValues
from checkm.customHelpFormatter import CustomHelpFormatter
from checkm.logger import logger_setup


class ChangeTempAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        if os.path.isdir(values):
            tempfile.tempdir = values
        else:
            raise argparse.ArgumentTypeError(
                'The value of %s must be a valid directory' % option_string)


def version():
    import checkm
    versionFile = open(os.path.join(checkm.__path__[0], 'VERSION'))
    return versionFile.readline().strip()


def printHelp():
    print('')
    print('                ...::: CheckM v' + version() + ' :::...''')
    print('''\

  Utility functions:
    coverage     -> Calculate coverage of sequences
    profile      -> Calculate percentage of reads mapped to each bin

  Use 'checkm data setRoot <checkm_data_dir>' to specify the location of CheckM database files.

  Usage: checkm <command> -h for command specific help
    ''')


if __name__ == '__main__':
    # initialize the options parser
    parser = argparse.ArgumentParser('checkm', add_help=False)
    subparsers = parser.add_subparsers(help="--", dest='subparser_name')

    data_parser = subparsers.add_parser('data',
                                        formatter_class=CustomHelpFormatter,
                                        description='Set path to the CheckM database files.',
                                        epilog='Example: checkm data setRoot')
    data_parser.add_argument('action', nargs="+",
                             help='''b
  setRoot <DB_PATH> -> set the the checkm data directory to <DB_PATH> (requires permissions)
  or
  setRoot <DB_PATH> <CONFIG_PATH> -> set the the checkm data directory to <DB_PATH> and the checkM DATA_CONFIG FILE to CONFIG_PATH (requires permissions)
            ''')

    # Calculate coverage
    coverage_parser = subparsers.add_parser('coverage',
                                            formatter_class=CustomHelpFormatter,
                                            description='Calculate coverage of sequences.',
                                            epilog='Example: checkm coverage ./bins coverage.tsv example_1.bam example_2.bam')
    coverage_parser.add_argument(
        'bin_input', help="directory containing bins (fasta format) or "
                          "path to file describing genomes/genes - tab "
                          "separated in 2 or 3 columns [genome ID, genome fna, genome translation file (pep)]")
    coverage_parser.add_argument('output_file', help="print results to file")
    coverage_parser.add_argument(
        'bam_files', nargs='+', help="BAM files to parse")
    coverage_parser.add_argument('-x', '--extension', default='fna',
                                 help="extension of bins (other files in directory are ignored)")
    coverage_parser.add_argument('-r', '--all_reads', action='store_true',
                                 help="use all reads to estimate coverage instead of just those in proper pairs")
    coverage_parser.add_argument(
        '-a', '--min_align', help='minimum alignment length as percentage of read length', type=float, default=0.98)
    coverage_parser.add_argument(
        '-e', '--max_edit_dist', help='maximum edit distance as percentage of read length', type=float, default=0.02)
    coverage_parser.add_argument(
        '-m', '--min_qc', help='minimum quality score (in phred)', type=int, default=15)
    coverage_parser.add_argument(
        '-t', '--threads', type=int, default=1, help="number of threads")
    coverage_parser.add_argument('-q', '--quiet', dest='bQuiet',
                                 action="store_true", default=False, help="suppress console output")

    # Calculate community profile
    profile_parser = subparsers.add_parser('profile',
                                           formatter_class=CustomHelpFormatter,
                                           description='Calculate percentage of reads mapped to each bin.',
                                           epilog='Example: checkm profile coverage.tsv')
    profile_parser.add_argument(
        'coverage_file', help="file indicating coverage of each sequence (see coverage command)")
    profile_parser.add_argument(
        '-f', '--file', default='stdout', help="print results to file")
    profile_parser.add_argument('--tab_table', dest='bTabTable', action="store_true",
                                default=False, help="print tab-separated values table")
    profile_parser.add_argument('-q', '--quiet', dest='bQuiet',
                                action="store_true", default=False, help="suppress console output")

    # get and check options
    args = None
    if(len(sys.argv) == 1 or sys.argv[1] == '-h' or sys.argv == '--help'):
        printHelp()
        sys.exit(0)
    else:
        args = parser.parse_args()

    # setup logging
    silent = False
    if hasattr(args, 'bQuiet') and args.bQuiet:
        silent = True

    try:
        logger_setup(args.output_dir, "checkm.log",
                     "CheckM", version(), silent)
    except:
        logger_setup(None, "checkm.log", "CheckM", version(), silent)

    # do what we came here to do
    try:
        logger = logging.getLogger('timestamp')
        logger.info(f"CheckM data: {DefaultValues.CHECKM_DATA_DIR}")

        checkmParser = main.OptionsParser()
        if(False):
            # import pstats
            # p = pstats.Stats('prof')
            # p.sort_stats('cumulative').print_stats(10)
            # p.sort_stats('time').print_stats(10)
            import cProfile
            cProfile.run('checkmParser.parseOptions(args)', 'prof')
        elif False:
            import pdb
            pdb.run(checkmParser.parseOptions(args))
        else:
            checkmParser.parseOptions(args)
    except SystemExit:
        print("\n  Controlled exit resulting from an unrecoverable error or warning.")
        raise
    except:
        print("\nUnexpected error: %s" % (sys.exc_info()[0]))
        raise
