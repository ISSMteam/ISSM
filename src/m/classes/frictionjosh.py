import numpy as np

from checkfield import checkfield
from fielddisplay import fielddisplay
from project3d import project3d
from WriteData import WriteData


class frictionjosh(object):
    """FRICTIONJOSH class definition

    Usage:
        frictionjosh = frictionjosh()
    """

    def __init__(self):  # {{{
        self.coefficient                   = np.nan
        self.pressure_adjusted_temperature = np.nan
        self.gamma                         = 0.
        self.effective_pressure_limit      = 0.
        self.coefficient_max               = 0.

        self.setdefaultparameters()
    # }}}

    def __repr__(self):  # {{{
        s = 'Basal shear stress parameters: Sigma_b = coefficient^2 * Neff ^r * |u_b|^(s - 1) * u_b,\n'
        s += '(effective stress Neff = rho_ice * g * thickness + rho_water * g * base, r = q / p and s = 1 / p)\n'
        s += '{}\n'.format(fielddisplay(self, "coefficient", "friction coefficient [SI]"))
        s += '{}\n'.format(fielddisplay(self, 'pressure_adjusted_temperature', 'friction pressure_adjusted_temperature (T - Tpmp) [K]'))
        s += '{}\n'.format(fielddisplay(self, 'gamma', '(T - Tpmp)/gamma [K]'))
        s += '{}\n'.format(fielddisplay(self, 'effective_pressure_limit', 'Neff do not allow to fall below a certain limit: effective_pressure_limit * rho_ice * g * thickness (default 0)'))
        s += '{}\n'.format(fielddisplay(self, 'coefficient_max', 'effective friction C = min(coefficient_max, sqrt(exp(T_b(modern) - T_b(t))/gamma) * coefficient)'))
        return s
    # }}}

    def extrude(self, md):  # {{{
        self.coefficient = project3d(md, 'vector', self.coefficient, 'type', 'node', 'layer', 1)
        self.pressure_adjusted_temperature = project3d(md, 'vector', self.pressure_adjusted_temperature, 'type', 'node', 'layer', 1)
        return self
    # }}}

    def setdefaultparameters(self):  # {{{

        #Default gamma: 1
        self.gamma = 1.

        #Default 0.
        self.effective_pressure_limit = 0.0

        #Default max friction coefficient: 300
        self.coefficient_max = 300.

        return self
    # }}}

    def checkconsistency(self, md, solution, analyses):  # {{{

        # Early return
        if 'StressbalanceAnalysis' not in analyses and 'ThermalAnalysis' not in analyses:
            return md

        md = checkfield(md, 'fieldname', 'friction.coefficient', 'timeseries', 1, 'NaN', 1, 'Inf', 1)
        md = checkfield(md, 'fieldname', 'friction.pressure_adjusted_temperature','NaN',1,'Inf',1)
        md = checkfield(md, 'fieldname', 'friction.gamma','numel',1,'NaN',1,'Inf',1,'>',0.)
        md = checkfield(md, 'fieldname', 'friction.effective_pressure_limit', 'numel', [1], '>=', 0)
        md = checkfield(md,'fieldname', 'friction.coefficient_max', 'numel', 1, 'NaN', 1, 'Inf', 1, '>', 0.)

        # Check that temperature is provided
        md = checkfield(md,'fieldname','initialization.temperature','NaN',1,'Inf',1,'size','universal')
        return md
    # }}}

    def marshall(self, prefix, md, fid):  # {{{
        WriteData(fid,prefix, 'name', 'md.friction.law', 'data',9, 'format', 'Integer')
        WriteData(fid,prefix, 'class', 'friction', 'object',self, 'fieldname', 'coefficient', 'format', 'DoubleMat', 'mattype',1, 'timeserieslength',md.mesh.numberofvertices+1, 'yts',md.constants.yts)
        WriteData(fid,prefix, 'class', 'friction', 'object',self, 'fieldname', 'pressure_adjusted_temperature', 'format', 'DoubleMat', 'mattype',1, 'timeserieslength',md.mesh.numberofvertices+1, 'yts',md.constants.yts)
        WriteData(fid,prefix, 'class', 'friction', 'object',self, 'fieldname', 'gamma', 'format', 'Double')
        WriteData(fid,prefix, 'object',self, 'class', 'friction', 'fieldname', 'effective_pressure_limit', 'format', 'Double')
        WriteData(fid,prefix, 'class', 'friction', 'object',self, 'fieldname', 'coefficient_max', 'format', 'Double')
    # }}}
