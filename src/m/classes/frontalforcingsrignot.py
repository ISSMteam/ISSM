# -*- coding: utf-8 -*-

import numpy as np

from checkfield import checkfield
from WriteData import WriteData
from fielddisplay import fielddisplay


class frontalforcingsrignot(object):
    """FRONTALFORCINGSRIGNOT class definition

    Usage:
        frontalforcingsrignot = frontalforcingsrignot()
    """

    def __init__(self, *args):  # {{{
        self.basin_id = np.nan
        self.num_basins = 0
        self.subglacial_discharge = np.nan
        self.thermalforcing = np.nan

        if len(args) == 0:
            self.setdefaultparameters()
        else:
            error('constructor not supported')

    # }}}

    def __repr__(self):  # {{{
        s = '   Frontalforcings parameters:\n'
        s += '{}\n'.format(fielddisplay(self, 'basin_id', 'basin ID for elements'))
        s += '{}\n'.format(fielddisplay(self, 'num_basins', 'number of basins'))
        s += '{}\n'.format(fielddisplay(self, 'subglacial_discharge', 'sum of subglacial discharge for each basin [m/d]'))
        s += '{}\n'.format(fielddisplay(self, 'thermalforcing', 'thermal forcing [∘C]'))
        return s
    # }}}

    def setdefaultparameters(self):  # {{{
        self.basin_id = np.nan
        self.num_basins = 0
        self.subglacial_discharge = np.nan
        self.thermalforcing = np.nan

        return self
    # }}}

    def extrude(self, md):  # {{{
        return self
    # }}}

    def checkconsistency(self, md, solution, analyses):  # {{{
        # Early return
        if (solution != 'TransientSolution') or (not md.transient.ismovingfront):
            return md

        md = checkfield(md, 'fieldname', 'frontalforcings.num_basins', 'numel', [1], 'NaN', 1, 'Inf', 1, '>', 0)
        md = checkfield(md, 'fieldname', 'frontalforcings.basin_id', 'Inf', 1, '>=', 0, '<=', md.frontalforcings.num_basins, 'size', [md.mesh.numberofelements])
        md = checkfield(md, 'fieldname', 'frontalforcings.subglacial_discharge', '>=', 0, 'NaN', 1, 'Inf', 1, 'timeseries', 1)
        md = checkfield(md, 'fieldname', 'frontalforcings.thermalforcing', 'NaN', 1, 'Inf', 1, 'timeseries', 1)
        return md
    # }}}

    def marshall(self, prefix, md, fid):  # {{{
        WriteData(fid, prefix, 'name', 'md.frontalforcings.parameterization', 'data', 2, 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'basin_id', 'data', self.basin_id - 1, 'name', 'md.frontalforcings.basin_id', 'format', 'IntMat', 'mattype', 2)  # 0-indexed
        WriteData(fid, prefix, 'object', self, 'fieldname', 'num_basins', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'subglacial_discharge', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', md.constants.yts)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'thermalforcing', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', md.constants.yts)
    # }}}
