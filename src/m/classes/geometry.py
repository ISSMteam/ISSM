import numpy as np

from checkfield import checkfield
from fielddisplay import fielddisplay
from project3d import project3d
from WriteData import WriteData


class geometry(object):
    """GEOMETRY class definition

    Usage:
        geometry = geometry()
    """

    def __init__(self, *args):  # {{{
        self.surface = np.nan
        self.thickness = np.nan
        self.base = np.nan
        self.bed = np.nan
        self.hydrostatic_ratio = np.nan

        if len(args) == 0:
            self.setdefaultparameters()
        else:
            raise Exception('constructor not supported')
    # }}}

    def __repr__(self):  # {{{
        s = '   geometry parameters:\n'
        s += '{}\n'.format(fielddisplay(self, 'surface', 'ice upper surface elevation [m]'))
        s += '{}\n'.format(fielddisplay(self, 'thickness', 'ice thickness [m]'))
        s += '{}\n'.format(fielddisplay(self, 'base', 'ice base elevation [m]'))
        s += '{}\n'.format(fielddisplay(self, 'bed', 'bed elevation [m]'))
        return s
    # }}}

    def setdefaultparameters(self):  # {{{
        return
    # }}}

    def checkconsistency(self, md, solution, analyses):  # {{{
        if solution == 'LoveSolution':
            return md
        else:
            md = checkfield(md, 'fieldname', 'geometry.surface', 'NaN', 1, 'Inf', 1, 'size', [md.mesh.numberofvertices])
            md = checkfield(md, 'fieldname', 'geometry.base', 'NaN', 1, 'Inf', 1, 'size', [md.mesh.numberofvertices])
            md = checkfield(md, 'fieldname', 'geometry.thickness', 'NaN', 1, 'Inf', 1, 'size', [md.mesh.numberofvertices], '>=', 0)
            if any(abs(self.thickness - self.surface + self.base) > 1e-9):
                md.checkmessage('equality thickness = surface-base violated')
            if solution == 'TransientSolution' and md.transient.isgroundingline:
                md = checkfield(md, 'fieldname', 'geometry.bed', 'NaN', 1, 'Inf', 1, 'size', [md.mesh.numberofvertices])
                if np.any(self.bed - self.base > 1e-12):
                    md.checkmessage('base < bed on one or more vertices')
                pos = np.where(md.mask.ocean_levelset > 0)
                if np.any(np.abs(self.bed[pos] - self.base[pos]) > 1e-9):
                    md.checkmessage('equality base = bed on grounded ice violated')
                md = checkfield(md, 'fieldname', 'geometry.bed', 'NaN', 1, 'Inf', 1, 'size', [md.mesh.numberofvertices])
        return md
    # }}}

    def marshall(self, prefix, md, fid):  # {{{
        if isinstance(self.thickness, (list, np.ndarray)):
            length_thickness = len(self.thickness)
        else:
            length_thickness = 1

        if (length_thickness == md.mesh.numberofvertices) or (length_thickness == md.mesh.numberofvertices + 1):
            WriteData(fid, prefix, 'object', self, 'fieldname', 'thickness', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', md.constants.yts)
        elif (length_thickness == md.mesh.numberofelements) or (length_thickness == md.mesh.numberofelements + 1):
            WriteData(fid, prefix, 'object', self, 'fieldname', 'thickness', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', md.mesh.numberofelements + 1, 'yts', md.constants.yts)
        else:
            raise RuntimeError('geometry thickness time series should be a vertex or element time series')

        WriteData(fid, prefix, 'object', self, 'fieldname', 'surface', 'format', 'DoubleMat', 'mattype', 1)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'base', 'format', 'DoubleMat', 'mattype', 1)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'bed', 'format', 'DoubleMat', 'mattype', 1)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'hydrostatic_ratio', 'format', 'DoubleMat', 'mattype', 1)
    # }}}

    def extrude(self, md):  # {{{
        self.surface = project3d(md, 'vector', self.surface, 'type', 'node')
        self.thickness = project3d(md, 'vector', self.thickness, 'type', 'node')
        self.hydrostatic_ratio = project3d(md, 'vector', self.hydrostatic_ratio, 'type', 'node')
        self.base = project3d(md, 'vector', self.base, 'type', 'node')
        self.bed = project3d(md, 'vector', self.bed, 'type', 'node')
        return self
    # }}}
