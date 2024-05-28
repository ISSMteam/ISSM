import numpy as np

from checkfield import checkfield
from fielddisplay import fielddisplay
from project3d import project3d
from structtoobj import structtoobj
from WriteData import WriteData


class frictionwaterlayer(object):
    """FRICTIONWATERLAYER class definition

    Usage:
        friction = frictionwaterlayer(md)
    """

    def __init__(self, *args):  # {{{
        self.coefficient = np.nan
        self.f = np.nan
        self.p = np.nan
        self.q = np.nan
        self.water_layer = np.nan

        nargs = len(args)
        if nargs == 0:
            self.setdefaultparameters()
        elif nargs == 1:
            self = structtoobj(self, args[0])
    # }}}

    def __repr__(self):  # {{{
        s = 'Basal shear stress parameters: tau_b = coefficient^2 * Neff ^r * |u_b|^(s - 1) * u_b * 1 / f(T)\n(effective stress Neff = rho_ice * g * thickness + rho_water * g * (bed + water_layer), r = q / p and s = 1 / p)\n'
        s = "{}\n".format(fielddisplay(self, 'coefficient', 'frictiontemp coefficient [SI]'))
        s = "{}\n".format(fielddisplay(self, 'f', 'f variable for effective pressure'))
        s = "{}\n".format(fielddisplay(self, 'p', 'p exponent'))
        s = "{}\n".format(fielddisplay(self, 'q', 'q exponent'))
        s = "{}\n".format(fielddisplay(self, 'water_layer', 'water thickness at the base of the ice (m)'))
        return s
    # }}}

    def setdefaultparameters(self):  # {{{
        return self
    # }}}

    def checkconsistency(self, md, solution, analyses):  #{{{
        # Early return
        if 'StressbalanceAnalysis' not in analyses and 'ThermalAnalysis' not in analyses:
            return
        md = checkfield(md, 'fieldname', 'friction.coefficient', 'timeseries', 1, 'NaN', 1, 'Inf', 1)
        md = checkfield(md, 'fieldname', 'friction.f', 'size', [1], 'NaN', 1, 'Inf', 1)
        md = checkfield(md, 'fieldname', 'friction.q', 'NaN', 1, 'Inf', 1, 'size', [md.mesh.numberofelements])
        md = checkfield(md, 'fieldname', 'friction.p', 'NaN', 1, 'Inf', 1, 'size', [md.mesh.numberofelements])
        md = checkfield(md, 'fieldname', 'thermal.spctemperature', 'Inf', 1, 'timeseries', 1, '>=', 0.)
        return md
    # }}}

    def extrude(self, md):  # {{{
        self.coefficient = project3d(md, 'vector', self.coefficient, 'type', 'node', 'layer', 1)
        self.p = project3d(md, 'vector', self.p, 'type', 'element')
        self.q = project3d(md, 'vector', self.q, 'type', 'element')
        self.water_layer = project3d(md, 'vector', self.water_layer, 'type', 'node', 'layer', 1)
        return self
    # }}}

    def marshall(self, prefix, md, fid):  #{{{
        WriteData(fid, prefix, 'name', 'md.friction.law', 'data', 5, 'format', 'Integer')
        WriteData(fid, prefix, 'class', 'friction', 'object', self, 'fieldname', 'coefficient', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', md.constants.yts)
        WriteData(fid, prefix, 'class', 'friction', 'object', self, 'fieldname', 'f', 'format', 'Double')
        WriteData(fid, prefix, 'class', 'friction', 'object', self, 'fieldname', 'p', 'format', 'DoubleMat', 'mattype', 2)
        WriteData(fid, prefix, 'class', 'friction', 'object', self, 'fieldname', 'q', 'format', 'DoubleMat', 'mattype', 2)
        WriteData(fid, prefix, 'class', 'friction', 'object', self, 'fieldname', 'water_layer', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', md.constants.yts)
    # }}}
