import numpy as np

from checkfield import checkfield
from fielddisplay import fielddisplay
from project3d import project3d
from WriteData import WriteData


class frictioncoulomb2(object):
    """frictioncoulomb2 class definition

    Usage:
        frictioncoulomb2 = frictioncoulomb2()
    """

    def __init__(self):  # {{{
        self.coefficient = np.nan
        self.coefficientcoulomb = np.nan
        self.p = np.nan
        self.q = np.nan
        self.coupling = 0
        self.effective_pressure = np.nan
        self.effective_pressure_limit = 0

        self.setdefaultparameters()
    # }}}

    def __repr__(self):  # {{{
        s = 'Basal shear stress parameters: Sigma_b = min(coefficient^2 * Neff ^r * |u_b|^(s - 1) * u_b,\n'
        s += 'coefficientcoulomb^2 * Neff), (effective stress Neff = rho_ice * g * thickness + rho_water * g * bed, r = q / p and s = 1 / p).\n'
        s += '{}\n'.format(fielddisplay(self, "coefficient", "power law (Weertman) friction coefficient [SI]"))
        s += '{}\n'.format(fielddisplay(self, "coefficientcoulomb", "Coulomb friction coefficient [SI]"))
        s += '{}\n'.format(fielddisplay(self, "p", "p exponent"))
        s += '{}\n'.format(fielddisplay(self, "q", "q exponent"))
        s += '{}\n'.format(fielddisplay(self, 'coupling', 'Coupling flag: 0 for default, 1 for forcing(provide md.friction.effective_pressure)  and 2 for coupled(not implemented yet)'))
        s += '{}\n'.format(fielddisplay(self, 'effective_pressure', 'Effective Pressure for the forcing if not coupled [Pa]'))
        s += '{}\n'.format(fielddisplay(self, 'effective_pressure_limit', 'Neff do not allow to fall below a certain limit: effective_pressure_limit * rho_ice * g * thickness (default 0)'))
        return s
    # }}}

    def setdefaultparameters(self):  # {{{
        self.effective_pressure_limit = 0
        return self
    # }}}

    def extrude(self, md):  # {{{
        self.coefficient = project3d(md, 'vector', self.coefficient, 'type', 'node', 'layer', 1)
        self.coefficientcoulomb = project3d(md, 'vector', self.coefficientcoulomb, 'type', 'node', 'layer', 1)
        self.p = project3d(md, 'vector', self.p, 'type', 'element')
        self.q = project3d(md, 'vector', self.q, 'type', 'element')
        if self.coupling == 1:
            self.effective_pressure = project3d(md, 'vector', self.effective_pressure, 'type', 'node', 'layer', 1)
        elif self.coupling == 2:
            raise ValueError('not implemented yet')
        elif self.coupling > 2:
            raise ValueError('not supported yet')
        return self
    # }}}
    def checkconsistency(self, md, solution, analyses):  # {{{
        # Early return
        if 'StressbalanceAnalysis' not in analyses and 'ThermalAnalysis' not in analyses:
            return md
        md = checkfield(md, 'fieldname', 'friction.coefficient', 'timeseries', 1, 'NaN', 1, 'Inf', 1)
        md = checkfield(md, 'fieldname', 'friction.coefficientcoulomb', 'timeseries', 1, 'NaN', 1, 'Inf', 1)
        md = checkfield(md, 'fieldname', 'friction.q', 'NaN', 1, 'Inf', 1, 'size', [md.mesh.numberofelements])
        md = checkfield(md, 'fieldname', 'friction.p', 'NaN', 1, 'Inf', 1, 'size', [md.mesh.numberofelements])
        md = checkfield(md, 'fieldname', 'friction.coupling', 'numel', [1], 'values', [0, 1, 2])
        md = checkfield(md, 'fieldname', 'friction.effective_pressure_limit', 'numel', [1], '>=', 0)
        if self.coupling == 1:
            md = checkfield(md, 'fieldname', 'friction.effective_pressure', 'NaN', 1, 'Inf', 1, 'timeseries', 1)
        elif self.coupling == 2:
            raise ValueError('not implemented yet')
        elif self.coupling > 2:
            raise ValueError('not supported yet')
        return md
    # }}}

    def marshall(self, prefix, md, fid):  # {{{
        WriteData(fid, prefix, 'name', 'md.friction.law', 'data', 7, 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'coefficient', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', md.constants.yts)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'coefficientcoulomb', 'format', 'DoubleMat', 'mattype', 1)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'p', 'format', 'DoubleMat', 'mattype', 2)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'q', 'format', 'DoubleMat', 'mattype', 2)
        WriteData(fid, prefix, 'class', 'friction', 'object', self, 'fieldname', 'coupling', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'class', 'friction', 'fieldname', 'effective_pressure_limit', 'format', 'Double')
        if self.coupling == 1:
            WriteData(fid, prefix, 'class', 'friction', 'object', self, 'fieldname', 'effective_pressure', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', md.constants.yts)
        elif self.coupling == 2:
            raise ValueError('not implemented yet')
        elif self.coupling > 2:
            raise ValueError('not supported yet')
    # }}}
