import numpy as np
from project3d import project3d
from fielddisplay import fielddisplay
from checkfield import checkfield
from WriteData import WriteData


class frictionhydro(object):
    """
    friction hydro is the friction law from Schoof 2005 or Gagliardini2007

    Usage:
        friction = frictionhydro()
    """
    def __init__(self):  # {{{
        self.coupling = 0
        self.q = np.nan
        self.C = np.nan
        self.As = np.nan
        self.effective_pressure = np.nan
        self.effective_pressure_limit = 0
        #set defaults
        self.setdefaultparameters()
    # }}}

    def __repr__(self):  # {{{
        string = 'Effective Pressure based friction law described in Gagliardini 2007'
        string = "%s\n%s" % (string, fielddisplay(self, 'coupling', 'Coupling flag 0: uniform sheet (negative pressure ok, default), 1: ice pressure only, 2: water pressure assuming uniform sheet (no negative pressure), 3: use provided effective_pressure, 4: used coupled model (not implemented yet)'))
        string = "%s\n%s" % (string, fielddisplay(self, 'q', 'friction law exponent q >= 1'))
        string = "%s\n%s" % (string, fielddisplay(self, 'C', 'friction law max value (Iken bound)'))
        string = "%s\n%s" % (string, fielddisplay(self, 'As', 'Sliding Parameter without cavitation [m Pa^ - n s^ - 1]'))
        string = "%s\n%s" % (string, fielddisplay(self, 'effective_pressure', 'Effective Pressure for the forcing if not coupled [Pa]'))
        string = "%s\n%s" % (string, fielddisplay(self, 'effective_pressure_limit', 'Neff do not allow to fall below a certain limit: effective_pressure_limit * rho_ice * g * thickness (default 0)'))

        return string
    # }}}

    def extrude(self, md):  # {{{
        self.q = project3d(md, 'vector', self.q, 'type', 'element')
        self.C = project3d(md, 'vector', self.C, 'type', 'element')
        self.As = project3d(md, 'vector', self.As, 'type', 'element')
        if self.coupling in[3, 4]:
            self.effective_pressure = project3d(md, 'vector', self.effective_pressure, 'type', 'node', 'layer', 1)
        elif self.coupling > 4:
            raise ValueError('md.friction.coupling larger than 4, not supported yet')
        return self
    # }}}

    def setdefaultparameters(self):  # {{{
        self.coupling = 0
        self.effective_pressure = np.nan
        self.effective_pressure_limit = 0
        return self
    # }}}

    def checkconsistency(self, md, solution, analyses):  #{{{
        #Early return
        if 'StressbalanceAnalysis' in analyses and 'ThermalAnalysis' in analyses:
            return md

        md = checkfield(md, 'fieldname', 'friction.coupling', 'numel', [1], 'values', [0, 1, 2, 3, 4])
        md = checkfield(md, 'fieldname', 'friction.q', 'NaN', 1, 'Inf', 1, 'size', [md.mesh.numberofelements])
        md = checkfield(md, 'fieldname', 'friction.C', 'NaN', 1, 'Inf', 1, 'size', [md.mesh.numberofelements])
        md = checkfield(md, 'fieldname', 'friction.As', 'NaN', 1, 'Inf', 1, 'size', [md.mesh.numberofelements])
        md = checkfield(md, 'fieldname', 'friction.effective_pressure_limit', 'numel', [1], '>=', 0)
        if self.coupling == 3:
            md = checkfield(md, 'fieldname', 'friction.effective_pressure', 'NaN', 1, 'Inf', 1, 'timeseries', 1)
        elif self.coupling > 4:
            raise ValueError('md.friction.coupling larger than 4, not supported yet')
    # }}}

    def marshall(self, prefix, md, fid):  #{{{
        WriteData(fid, prefix, 'name', 'md.friction.law', 'data', 3, 'format', 'Integer')
        WriteData(fid, prefix, 'class', 'friction', 'object', self, 'fieldname', 'coupling', 'format', 'Integer')
        WriteData(fid, prefix, 'class', 'friction', 'object', self, 'fieldname', 'q', 'format', 'DoubleMat', 'mattype', 2)
        WriteData(fid, prefix, 'class', 'friction', 'object', self, 'fieldname', 'C', 'format', 'DoubleMat', 'mattype', 2)
        WriteData(fid, prefix, 'class', 'friction', 'object', self, 'fieldname', 'As', 'format', 'DoubleMat', 'mattype', 2)
        WriteData(fid, prefix, 'object', self, 'class', 'friction', 'fieldname', 'effective_pressure_limit', 'format', 'Double')
        if self.coupling in[3, 4]:
            WriteData(fid, prefix, 'class', 'friction', 'object', self, 'fieldname', 'effective_pressure', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', md.constants.yts)
        elif self.coupling > 4:
            raise ValueError('md.friction.coupling larger than 4, not supported yet')
    # }}}
