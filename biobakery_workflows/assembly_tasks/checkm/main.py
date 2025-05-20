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

import os
import sys
import logging

from checkm.defaultValues import DefaultValues
from checkm.timeKeeper import TimeKeeper
from checkm.coverage import Coverage
from checkm.profile import Profile
from checkm.common import (checkBinInputExists,
                           makeSurePathExists,
                           checkFileExists,
                           binIdFromFilename,
                           getBinIdsFromOutDir,
                           checkDirExists)

from checkm.checkmData import DBManager

class OptionsParser():
    def __init__(self):
        self.logger = logging.getLogger('timestamp')
        self.timeKeeper = TimeKeeper()

    def updateCheckM_DB(self, options):
        self.logger.info(
            '[CheckM - data] Check for database updates. [%s]' % options.action[0])

        action = options.action
        if action and action[0] == 'setRoot':
            if len(action) > 2:
                DBM = DBManager(set_path=action[1], configFile=action[2])
            elif len(action) > 1:
                DBM = DBManager(set_path=action[1])
        else:
            self.logger.error(
                'Path to the CheckM reference data must be specified.')

    def binFiles(self, binInput, binExtension, bCalledGenes):
        binFiles = []
        binIDs = set()
        isInputDir = True
        if binInput is not None:
            if os.path.isdir(binInput):
                if binExtension[0] != '.':
                    binExtension = '.' + binExtension

                all_files = os.listdir(binInput)
                for f in all_files:
                    if f.endswith(binExtension):
                        binFile = os.path.join(binInput, f)
                        if os.stat(binFile).st_size == 0:
                            self.logger.warning(
                                "Skipping bin %s as it has a size of 0 bytes." % f)
                        else:
                            binFiles.append(binFile)
                            binIDs.add(os.path.basename(binFile))
            else:
                with open(binInput, "r") as oh:
                    for line in oh:
                        files = line.strip().split("\t")
                        binFile = files[1]
                        if bCalledGenes:
                            binFile = files[2]
                        if not os.path.exists(binFile):
                            self.logger.warning(
                                "Skipping bin %s as it doesn't exists." % binFile)
                        elif os.stat(binFile).st_size == 0:
                            self.logger.warning(
                                "Skipping bin %s as it has a size of 0 bytes." % binFile)
                        else:
                            binFiles.append(binFile)
                            binIDs.add(os.path.basename(binFile))

        if not binFiles:
            if isInputDir:
                self.logger.error(
                    "No bins found. Check the extension (-x) used to identify bins.")
            else:
                self.logger.error(
                    "No bins found. Check the bins input table to verify bins exists.")
            sys.exit(1)

        if len(binIDs) != len(binFiles):
            self.logger.error(
                "There are redundant bin IDs, please check and update.")
            sys.exit(1)

        return sorted(binFiles)

    def coverage(self, options):
        """Coverage command"""

        self.logger.info(
            '[CheckM - coverage] Calculating coverage of sequences.')

        checkBinInputExists(options.bin_input, False)
        makeSurePathExists(os.path.dirname(options.output_file))

        binFiles = self.binFiles(
            options.bin_input, options.extension, False)

        coverage = Coverage(options.threads)
        coverage.run(binFiles, options.bam_files, options.output_file, options.all_reads,
                     options.min_align, options.max_edit_dist, options.min_qc)

        self.logger.info(
            'Coverage information written to: ' + options.output_file)

        self.timeKeeper.printTimeStamp()

    def profile(self, options):
        """Profile command"""

        self.logger.info(
            '[CheckM - profile] Calculating percentage of reads mapped to each bin.')

        checkFileExists(options.coverage_file)

        profile = Profile()
        profile.run(options.coverage_file, options.file, options.bTabTable)

        if options.file != '':
            self.logger.info('Profile information written to: ' + options.file)

        self.timeKeeper.printTimeStamp()

    def parseOptions(self, options):
        """Parse user options and call the correct pipeline(s)"""

        try:
            if options.file == "stdout":
                options.file = ''
        except:
            pass

        if(options.subparser_name == "data"):
            self.updateCheckM_DB(options)
        elif(options.subparser_name == 'coverage'):
            self.coverage(options)
        elif(options.subparser_name == 'profile'):
            self.profile(options)
        else:
            self.logger.error('Unknown CheckM command: ' +
                              options.subparser_name + '\n')
            sys.exit(1)

        return 0
