from fielddisplay import fielddisplay
from checkfield import checkfield
from WriteData import WriteData
from project3d import project3d


class frictionweertman(object):
    """
    FRICTIONWEERTMAN class definition

       Usage:
          frictionweertman = frictionweertman()
    """

    def __init__(self):  # {{{
        self.C = float('NaN')
        self.m = float('NaN')
        self.linearize = 0

    #set defaults
        self.setdefaultparameters()

    # }}}

    def __repr__(self):  # {{{
        string = "Weertman sliding law parameters: Sigma_b = C^(- 1 / m) * |u_b|^(1 / m - 1) * u_b"

        string = "%s\n%s" % (string, fielddisplay(self, "C", "friction coefficient [SI]"))
        string = "%s\n%s" % (string, fielddisplay(self, "m", "m exponent"))
        string = "%s\n%s" % (string, fielddisplay(self, "linearize", "0: not linearized, 1: interpolated linearly, 2: constant per element (default is 0)"))
        return string
    # }}}

    def extrude(self,md): # {{{
        print('-------------- file: frictionweertman.m line: 35')
        self.C=project3d(md,'vector',self.C,'type','node','layer',1)
        self.m=project3d(md,'vector',self.m,'type','element')
        return self
        # }}}

    def setdefaultparameters(self):  # {{{
        self.linearize = 0
        return self
    # }}}

    def checkconsistency(self, md, solution, analyses):  # {{{

        #Early return
        if 'StressbalanceAnalysis' not in analyses and 'ThermalAnalysis' not in analyses:
            return md

        md = checkfield(md, 'fieldname', 'friction.C', 'timeseries', 1, 'NaN', 1, 'Inf', 1)
        md = checkfield(md, 'fieldname', 'friction.m', 'NaN', 1, 'Inf', 1, 'size', [md.mesh.numberofelements])
        md = checkfield(md, 'fieldname', 'friction.linearize', 'numel', [1], 'values', [0, 1, 2])

        return md
    # }}}

    def marshall(self, prefix, md, fid):  # {{{
        WriteData(fid, prefix, 'name', 'md.friction.law', 'data', 2, 'format', 'Integer')
        WriteData(fid, prefix, 'class', 'friction', 'object', self, 'fieldname', 'C', 'format', 'DoubleMat', 'mattype', 1)
        WriteData(fid, prefix, 'class', 'friction', 'object', self, 'fieldname', 'm', 'format', 'DoubleMat', 'mattype', 2)
        WriteData(fid, prefix, 'class', 'friction', 'object', self, 'fieldname', 'linearize', 'format', 'Integer')
    # }}}
