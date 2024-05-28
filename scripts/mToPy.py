#!/usr/bin/env python
import sys
import os
import shutil
import translateToPy
import mToPy  # touch mToPy to assertain location of installation
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
program = 'mToPy.py'
version = '1.0'
versionReleaseDate = '09/24/12'
origAuthor = 'Mike Pellegrin'
desc = '\nMain control unit for converting an matlab script file to python'
#
#   Note: Output will be put in the same (absolute) location as the input.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                History
#    Date        Developer           Modification
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#    09 / 24 / 12    Michael Pellegrin    Initial Release         V1.0
#
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


def convert(inputFile):
    try:
        if os.path.exists(inputFile + '.m') and os.path.isfile(inputFile + '.m'):
            checkBackupOutputFile(inputFile)
            convertMToPy(inputFile)
        else:
            print('Specified input file: ' + inputFile + '.m doesn\'t appear to exist')
    finally:
        print('')


def convertMToPy(inputFileName):
    translateToPy.convertToPython(inputFileName + '.m', inputFileName + '.py')


def checkBackupOutputFile(inputFile):
    mFile = inputFile + '.m'
    pyFile = inputFile + '.py'
    if os.path.exists(pyFile):
        i = 1
        bkupName = pyFile + str(i)
        while os.path.exists(bkupName):
            i += 1
            bkupName = pyFile + str(i)
        os.rename(pyFile, bkupName)

    shutil.copyfile(mFile, pyFile)


if __name__ == "__main__":
    convert(sys.argv[1])
