import numpy as np

from checkfield import checkfield
from fielddisplay import fielddisplay
from WriteData import WriteData


class surfaceload(object):
    """SURFACELOAD class definition

    Usage:
        surfaceload = surfaceload()
    """

    def __init__(self, *args):  #{{{
        self.icethicknesschange = np.nan
        self.waterheightchange = np.nan
        self.other = np.nan

        nargin = len(args)

        if nargin == 0:
            self.setdefaultparameters()
        else:
            raise Exception('constructor not supported')
    # }}}

    def __repr__(self):  #{{{
        s = '   surfaceload:\n'
        s += '{}\n'.format(fielddisplay(self, 'icethicknesschange', 'thickness change: ice height equivalent [mIce/yr]'))
        s += '{}\n'.format(fielddisplay(self, 'waterheightchange', 'water height change: water height equivalent [mWater/yr]'))
        s += '{}\n'.format(fielddisplay(self, 'other', 'other loads (sediments) [kg/m^2/yr]'))
        return s
    # }}}

    def setdefaultparameters(self):  # {{{
        return self
    # }}}

    def checkconsistency(self, md, solution, analyses):  # {{{
        if ('SealevelchangeAnalysis' not in analyses) or (solution == 'TransientSolution' and not md.transient.isslc):
            return md
        if type(self.icethicknesschange) == np.ndarray:
            md = checkfield(md, 'fieldname', 'solidearth.surfaceload.icethicknesschange', 'timeseries', 1, 'NaN', 1, 'Inf', 1)
        if type(self.waterheightchange) == np.ndarray:
            md = checkfield(md, 'fieldname', 'solidearth.surfaceload.waterheightchange', 'timeseries', 1, 'NaN', 1, 'Inf', 1)
        if type(self.other) == np.ndarray:
            md = checkfield(md, 'fieldname', 'solidearth.surfaceload.other', 'timeseries', 1, 'NaN', 1, 'Inf', 1)
        return md
    # }}}

    def marshall(self, prefix, md, fid):  #{{{
        # TODO
        # - When surfaceload class is used eventually, check premarshall 
        # modifications against those done in MATLAB
        #

        # Deal with ice thickness change
        if np.isempty(self.icethicknesschange):
            self.icethicknesschange = np.zeros((md.mesh.numberofelements + 1, ))

        yts = md.constants.yts

        WriteData(fid, prefix, 'object', self, 'fieldname', 'icethicknesschange', 'name', 'md.solidearth.surfaceload.icethicknesschange', 'format', 'MatArray', 'timeserieslength', md.mesh.numberofelements + 1, 'yts', yts, 'scale', 1 / yts)

        # Deal with ice thickness change
        if np.isempty(self.waterheightchange):
            self.waterheightchange = np.zeros((md.mesh.numberofelements + 1, ))

        WriteData(fid, prefix, 'object', self, 'fieldname', 'waterheightchange', 'name', 'md.solidearth.surfaceload.waterheightchange', 'format', 'MatArray', 'timeserieslength', md.mesh.numberofelements + 1, 'yts', yts, 'scale', 1 / yts)

        # Deal with ice thickness change
        if np.isempty(self.otherchange):
            self.otherchange = np.zeros((md.mesh.numberofelements + 1, ))

        WriteData(fid, prefix, 'object', self, 'fieldname', 'otherchange', 'name', 'md.solidearth.surfaceload.otherchange', 'format', 'MatArray', 'timeserieslength', md.mesh.numberofelements + 1, 'yts', yts, 'scale', 1 / yts)
    # }}}

    def extrude(self, md):  #{{{
        return self
    # }}}
