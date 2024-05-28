
import codecs
import unicodedata
import sys
import datetime
import os

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
program = 'translateToPy.py'
version = '1.0'
versionReleaseDate = '09/24/12'
origAuthor = 'Mike Pellegrin'
desc = '\nMatlab script conversion into python'
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                History
#    Date        Developer           Modification
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#    09 / 24 / 12 Michael Pellegrin    Initial Release     V1.0
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


outputLocation = sys.stdout
inputFile = ""

# other global vars
indentLevel = 0


def setupOutputLocation(outFile):
    if outFile != sys.stdout:
        globals()['outputLocation'] = open(outFile, 'w')  # clobber


def translateFile(inputFile):
    f = codecs.open(inputFile, encoding='utf-8')
    try:
        for line in f:
            # print "in: " + line

            asciiLine = unicodedata.normalize('NFKD', line).encode('ascii', 'ignore')
            line = asciiLine
            translateLine(line)

    finally:
        f.close()


def translateLine(line):

    if len(line) == 1:    # blank line
        output(line)

    elif line.split()[0][0] == '%':        # comment line
        output("# " + line.replace('%', ''))

    else:        # needs cleanup.  this is a real - quick - n - dirty implimentation
        #print line
        res = line.replace('{', '[')
        res = res.replace('}', ']')
        res = res.replace('model', 'model()')
        res = res.replace('SolutionEnum', 'SolutionEnum()')
        res = res.replace('StressTensorEnum', 'StressTensorEnum()')
        res = res.replace('.par', '.py')
        res = res.replace('=extrude(md, ', '.extrude(')

        res = res.replace('thickness(pos)', 'thickness[pos]')
        res = res.replace('find(md.', 'numpy.nonzero(md.')

        res = res.replace('\n', '')

        # handle inline comments
        res = res.replace('%', '#')

        res = res.replace('...', '\\')

        # determine if the m file has mult. line cmd (real quick solution)
        multCmds = res.split(';')
        numLines = len(multCmds) - 2
        allParts = ''
        for part in multCmds:
            allParts += part
            #allParts += re.sub('^\s + ', '', part)
            #allParts += part.strip()
            if numLines > 0:
                allParts += '\n'
                numLines -= 1

        res = allParts
        res = res.replace(';', '')
        res = convertFieldValues(res)
        #print 'resulting line:' + str(res) + '\n'
        output(res)


def convertFieldValues(currentLine):
    # before utilizing regex's {starting w / eg. \([0 - 9]\) } for special case: ...(#)...
    # i noticed what i'm looking for is only TransientSolution(* ). So, ...

    res = currentLine
    if 'md.results' in currentLine:
        res = res.replace('(md.results.', 'md.results[\'')

        if 'TransientSolution(' in currentLine:        # got a TransientSolution([0 - 9..]) case
            res = res.replace('TransientSolution(', 'TransientSolution\'][')
            parts = res.split(')')
            res = parts[0] + '][\'' + parts[1].replace('.', '') + '\']' + parts[2]

        else:                # handle the other cases for md.results

            res = res.replace('Solution.Vx)', 'Solution\'][1][\'Vx\']')
            res = res.replace('Solution.Vy)', 'Solution\'][1][\'Vy\']')
            res = res.replace('Solution.Vz)', 'Solution\'][1][\'Vz\']')
            res = res.replace('Solution.Vel)', 'Solution\'][1][\'Vel\']')

            res = res.replace('Solution.Pressure)', 'Solution\'][1][\'Pressure\']')

            res = res.replace('Solution.StressTensorxx)', 'Solution\'][1][\'StressTensorxx\']')
            res = res.replace('Solution.StressTensorxy)', 'Solution\'][1][\'StressTensorxy\']')
            res = res.replace('Solution.StressTensoryy)', 'Solution\'][1][\'StressTensoryy\']')
            res = res.replace('Solution.StressTensorzz)', 'Solution\'][1][\'StressTensorzz\']')
            res = res.replace('Solution.StressTensorxz)', 'Solution\'][1][\'StressTensorxz\']')
            res = res.replace('Solution.StressTensoryz)', 'Solution\'][1][\'StressTensoryz\']')

            res = res.replace('Solution.FrictionCoefficient)', 'Solution\'][1][\'FrictionCoefficient\']')
            res = res.replace('Solution.SurfaceforcingsMasBalance)', 'Solution\'][1][\'SurfaceforcingsMasBalance\']')
            res = res.replace('Solution.MaskElementonfloatingice)', 'Solution\'][1][\'MaskElementonfloatingice\']')
            res = res.replace('Solution.J)', 'Solution\'][1][\'J\']')
            res = res.replace('Solution.BalancethicknessThickeningRate)', 'Solution\'][1][\'BalancethicknessThickeningRate\']')

            res = res.replace('Solution.Gradient1)', 'Solution\'][1][\'Gradient1\']')
            res = res.replace('Solution.Gradient2)', 'Solution\'][1][\'Gradient2\']')

            res = res.replace('Solution.MaterialsRheologyZbar)', 'Solution\'][1][\'MaterialsRheologyZbar\']')
            res = res.replace('Solution.MaterialsRheologyBbar)', 'Solution\'][1][\'MaterialsRheologyBbar\']')
            res = res.replace('Solution.MaterialsRheologyB)', 'Solution\'][1][\'MaterialsRheologyB\']')

            res = res.replace('Solution.Thickness)', 'Solution\'][1][\'Thickness\']')

            res = res.replace('Solution.Temperature)', 'Solution\'][1][\'Temperature\']')

            res = res.replace('Solution.BasalforcingsMeltingRate)', 'Solution\'][1][\'BasalforcingsMeltingRate\']')

            res = res.replace('Solution.SurfaceSlopeX)', 'Solution\'][1][\'SurfaceSlopeX\']')
            res = res.replace('Solution.SurfaceSlopeY)', 'Solution\'][1][\'SurfaceSlopeY\']')
            res = res.replace('Solution.SurfaceSlopeZ)', 'Solution\'][1][\'SurfaceSlopeZ\']')

            res = res.replace('Solution.BedSlopeX)', 'Solution\'][1][\'BedSlopeX\']')
            res = res.replace('Solution.BedSlopeY)', 'Solution\'][1][\'BedSlopeY\']')
            res = res.replace('Solution.BedSlopeZ)', 'Solution\'][1][\'BedSlopeZ\']')

            res = res.replace('Solution.Enthalpy)', 'Solution\'][1][\'Enthalpy\']')
            res = res.replace('Solution.Waterfraction)', 'Solution\'][1][\'Waterfraction\']')
            res = res.replace('Solution.Temperature)', 'Solution\'][1][\'Temperature\']')

            # special case
            res = res.replace('.DiagnosticSolution.J', '[\'DiagnosticSolution\'][1][\'J\']')

    return res


def output(line):
    numTabs = indentLevel
    while numTabs:
        numTabs -= 1
        print('\t', end="", file=outputLocation)

    print(line, end="", file=outputLocation)


def outputTopOfSript(inputFile):

    global indentLevel

    output("\"\"\"")
    output("====================================== ")
    output("Auto generated python script for ISSM: {}".format(inputFile))
    output("Created on {} via {} Ver {} by {}".format(datetime.date.today(), program, version, os.getlogin()))
    output("====================================== ")
    #output("")
    output(desc)
    output("%s Author: Michael Pellegrin" % (program))
    output("%s Date: %s" % (program, versionReleaseDate))
    output("====================================== ")
    output("\"\"\"")
    output("")


def outputBottomOfScript():

    global indentLevel

    output("")


def genericImports():
    output("from MatlabFuncs import * ")
    output("from model import * ")
    output("from EnumDefinitions import * ")
    output("from numpy import * ")


def createImports(inputFile):
    genericImports()

    # cycle through eachline to assertain import needs
    f = codecs.open(inputFile, encoding='utf-8')
    try:
        for line in f:
            # print "in: " + line

            # toss blank lines
            if len(line) == 1:
                continue

            asciiLine = unicodedata.normalize('NFKD', line).encode('ascii', 'ignore')
            line = asciiLine

            for il in importList:
                if line.find(il) != -1:
                    output("from %s import * " % (il))
                    importList.remove(il)    # already got it

    finally:
        output("")
        f.close()


def initImportList():
    global importList

    importList = ['triangle',
                  'setmask',
                  'parameterize',
                  'setflowequation',
                  'meshconvert',
                  'solve']


def convertToPython(inFile, outFile=sys.stdout):
    #print ' in cnvrt to python w / file:' + inFile
    initImportList()
    setupOutputLocation(outFile)
    outputTopOfSript(inFile)
    createImports(inFile)
    translateFile(inFile)
    #    outputBottomOfScript()


if __name__ == "__main__":
    #print ' in main w / arg:' + sys.argv[1] + ' ' + sys.argv[2]
    if len(sys.argv) == 2:
        convertToPython(sys.argv[1], sys.argv[2])
    else:
        convertToPython(sys.argv[1])
