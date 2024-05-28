import numpy as np

from checkfield import checkfield
from fielddisplay import fielddisplay
from project3d import project3d
from structtoobj import structtoobj
from WriteData import WriteData


class frictionshakti(object):
    """FRICTIONSHAKTI class definition

    Usage:
        friction = frictionshakti()
    """

    def __init__(self, *args):  # {{{
        self.coefficient = np.nan
        
        nargs = len(args)
        if nargs == 0:
            self.setdefaultparameters()
        elif nargs == 1:
            self = structtoobj(self, args[0])
    # }}}

    def __repr__(self):  # {{{
        s = 'Basal shear stress parameters: Sigma_b = coefficient^2 * Neff * u_b\n'
        s += '(effective stress Neff = rho_ice * g * thickness + rho_water * g * (head - b))\n'
        s += '{}\n'.format(fielddisplay(self, 'coefficient', 'friction coefficient [SI]'))
        return s
    # }}}

    def setdefaultparameters(self):  # {{{
        return self
    # }}}

    def extrude(self, md):  # {{{
        self.coefficient = project3d(md, 'vector', self.coefficient, 'type', 'node', 'layer', 1)
        return self
    # }}}

    def checkconsistency(self, md, solution, analyses):  # {{{
        # Early return
        if 'StressbalanceAnalysis' not in analyses and 'ThermalAnalysis' not in analyses:
            return md
        md = checkfield(md, 'fieldname', 'friction.coefficient', 'timeseries', 1, 'NaN', 1, 'Inf', 1)
        return md
    # }}}

    def marshall(self, prefix, md, fid):  # {{{
        WriteData(fid, prefix, 'name', 'md.friction.law', 'data', 8, 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'coefficient', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', md.constants.yts)
    # }}}
