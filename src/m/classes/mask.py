import numpy as np

from checkfield import checkfield
from fielddisplay import fielddisplay
from project3d import project3d
from WriteData import WriteData


class mask(object):
    """MASK class definition

    Usage:
        mask = mask()
    """
    
    def __init__(self):  # {{{
        self.ice_levelset = float('NaN')
        self.ocean_levelset = float('NaN')

        #set defaults
        self.setdefaultparameters()
    # }}}

    def __repr__(self):  # {{{
        string = "   masks:"

        string = "%s\n%s" % (string, fielddisplay(self, "ocean_levelset", "presence of ocean if < 0, coastline/grounding line if = 0, no ocean if > 0"))
        string = "%s\n%s" % (string, fielddisplay(self, "ice_levelset", "presence of ice if < 0, icefront position if = 0, no ice if > 0"))
        return string
    # }}}

    def setdefaultparameters(self):  # {{{
        return self
    # }}}

    def extrude(self, md):  # {{{
        self.ice_levelset = project3d(md, 'vector', self.ice_levelset, 'type', 'node')
        self.ocean_levelset = project3d(md, 'vector', self.ocean_levelset, 'type', 'node')

        return self
    # }}}

    def mask(self, *args):  # {{{
        if len(args) == 0:
            self.setdefaultparameters()
        else:
            raise RuntimeError("constructor not supported")

        return self
    # }}}

    def oceanset(self, *args):  # {{{
        if len(args) == 2:
            ocean = args[0]
            index = args[1]
            self.ocean_levelset[index] = -ocean
        elif len(args) == 1:
            ocean = args[0]
            self.ocean_levelset = -ocean
        else:
            raise RuntimeError("oceanset error message: not supported yet")

        return self
    # }}}

    def iceset(self, *args):  # {{{
        if len(args) == 2:
            ice = args[0]
            index = args[1]
            self.ice_levelset[index] = -ice
        elif len(args) == 1:
            ice = args[0]
            self.ocean_levelset = -ice
        else:
            raise RuntimeError("iceset error message: not supported yet")

        return self
    # }}}

    def checkconsistency(self, md, solution, analyses):  # {{{
        if solution == 'LoveSolution':
            return

        md = checkfield(md, 'fieldname', 'mask.ice_levelset', 'size', [md.mesh.numberofvertices])
        isice = np.array(md.mask.ice_levelset <= 0, int)
        if np.sum(isice) == 0:
            raise TypeError("no ice present in the domain")

        return md
    # }}}

    def marshall(self, prefix, md, fid):  # {{{
        WriteData(fid, prefix, 'object', self, 'fieldname', 'ocean_levelset', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', md.constants.yts)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'ice_levelset', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', md.constants.yts)
    # }}}
