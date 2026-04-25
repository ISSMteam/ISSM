import numpy as np
from fielddisplay import fielddisplay
from checkfield import checkfield
from WriteData import WriteData
from isnans import isnans
import MatlabFuncs as m


class rifts(object):
    """
    RIFTS class definition

       Usage:
          rifts = rifts()
    """

    def __init__(self):  # {{{
        self.riftstruct = []
        self.riftproperties = []

    #set defaults
        self.setdefaultparameters()

    # }}}

    def __repr__(self):  # {{{
        string = '   rifts parameters:'
        string += '\n{}'.format(fielddisplay(self, 'riftstruct', 'structure containing all rift information (vertices coordinates, segments, type of melange, ...)'))
        string += '\n{}'.format(fielddisplay(self, 'riftproperties', ''))
        return string
    # }}}

    def setdefaultparameters(self):  # {{{
        return self
    # }}}

    def checkconsistency(self, md, solution, analyses):  # {{{
        if (not self.riftstruct) or np.any(isnans(self.riftstruct)):
            numrifts = 0
        else:
            numrifts = len(self.riftstruct)

        if numrifts:
            if not m.strcmp(md.mesh.domaintype(), '2Dhorizontal'):
                md.checkmessage('models with rifts are only supported in 2d for now!')
            if not isinstance(self.riftstruct, list):
                md.checkmessage('rifts.riftstruct should be a list of dict!')
            if np.any(md.mesh.segmentmarkers >= 2):
                #We have segments with rift markers, but no rift structure!
                md.checkmessage('model should be processed for rifts (run meshprocessrifts)!')
            for i, rift in enumerate(self.riftstruct):
                md = checkfield(md, 'fieldname', 'rifts.riftstruct[{}][\'fill\']'.format(i), 'values', ['Water', 'Air', 'Ice', 'Melange', 0, 1, 2, 3])
        else:
            if self.riftstruct and np.any(np.logical_not(isnans(self.riftstruct))):
                md.checkmessage('riftstruct should be NaN since numrifts is 0!')

        return md
    # }}}

    def marshall(self, prefix, md, fid):  # {{{
        #Process rift info
        if (not self.riftstruct) or np.any(isnans(self.riftstruct)):
            numrifts = 0
        else:
            numrifts = len(self.riftstruct)

        numpairs = 0
        if numrifts > 0:
            for rift in self.riftstruct:
                numpairs += np.size(rift['penaltypairs'], axis=0)

            # Convert strings in riftstruct to hard coded numbers
            FillDict = {'Air': 0,
                        'Ice': 1,
                        'Melange': 2,
                        'Water': 3}
            for rift in self.riftstruct:
                if rift['fill'] in ['Air', 'Ice', 'Melange', 'Water']:
                    rift['fill'] = FillDict[rift['fill']]

            # 2 for nodes + 2 for elements + 2 for  normals + 1 for length + 1 for fill + 1 for friction + 1 for fraction + 1 for fractionincrement + 1 for state.
            data = np.zeros((numpairs, 12))
            count = 0
            for rift in self.riftstruct:
                numpairsforthisrift = np.size(rift['penaltypairs'], 0)
                data[count:count + numpairsforthisrift, 0:7] = rift['penaltypairs']
                data[count:count + numpairsforthisrift, 7] = rift['fill']
                data[count:count + numpairsforthisrift, 8] = rift['friction']
                data[count:count + numpairsforthisrift, 9] = rift['fraction']
                data[count:count + numpairsforthisrift, 10] = rift['fractionincrement']
                data[count:count + numpairsforthisrift, 11] = rift['state'].reshape(-1)
                count += numpairsforthisrift
        else:
            data = np.zeros((numpairs, 12))
        WriteData(fid, prefix, 'data', numrifts, 'name', 'md.rifts.numrifts', 'format', 'Integer')
        WriteData(fid, prefix, 'data', data, 'name', 'md.rifts.riftstruct', 'format', 'DoubleMat', 'mattype', 3)
    # }}}
