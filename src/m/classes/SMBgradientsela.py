from fielddisplay import fielddisplay
from checkfield import checkfield
from WriteData import WriteData


class SMBgradientsela(object):
    """SMBGRADIENTSELA class definition

    Usage:
        SMBgradientsela = SMBgradientsela()
    """

    def __init__(self, *args):  # {{{
        self.ela = float('NaN')
        self.b_pos = float('NaN')
        self.b_neg = float('NaN')
        self.b_max = float('NaN')
        self.b_min = float('NaN')
        self.steps_per_step = 1
        self.averaging = 0
        self.requested_outputs = []

        if len(args) == 0:
            self.setdefaultparameters()
        else:
            error('constructor not supported')
    # }}}

    def __repr__(self):  # {{{
        string = "   surface forcings parameters:"
        string += '\n   SMB gradients ela parameters:'

        string = "%s\n%s" % (string, fielddisplay(self, 'ela', ' equilibrium line altitude from which deviation is used to calculate smb using the smb gradients ela method [m a.s.l.]'))
        string = "%s\n%s" % (string, fielddisplay(self, 'b_pos', ' vertical smb gradient (dB/dz) above ela'))
        string = "%s\n%s" % (string, fielddisplay(self, 'b_neg', ' vertical smb gradient (dB/dz) below ela'))
        string = "%s\n%s" % (string, fielddisplay(self, 'b_max', ' upper cap on smb rate, default: 9999 (no cap) [m ice eq./yr]'))
        string = "%s\n%s" % (string, fielddisplay(self, 'b_min', ' lower cap on smb rate, default: -9999 (no cap) [m ice eq./yr]'))
        string = "%s\n%s" % (string, fielddisplay(self, 'steps_per_step', 'number of smb steps per time step'))
        string = "%s\n%s" % (string, fielddisplay(self, 'averaging', 'averaging methods from short to long steps'))
        string = "%s\n\t\t%s" % (string, '0: Arithmetic (default)')
        string = "%s\n\t\t%s" % (string, '1: Geometric')
        string = "%s\n\t\t%s" % (string, '2: Harmonic')

        string = "%s\n%s" % (string, fielddisplay(self, 'requested_outputs', 'additional outputs requested'))
        return string
    # }}}

    def extrude(self, md):  # {{{
        #Nothing for now
        return self
    # }}}

    def defaultoutputs(self, md):  # {{{
        return ['SmbMassBalance']
    # }}}

    def initialize(self, md):  # {{{
        #Nothing for now
        return self
    # }}}

    def setdefaultparameters(self):  # {{{
        self.b_max = 9999
        self.b_min = -9999
        self.requested_outputs = ['default']
        return self
    # }}}

    def checkconsistency(self, md, solution, analyses):    # {{{
        if 'MasstransportAnalysis' in analyses:
            md = checkfield(md, 'fieldname', 'smb.ela', 'timeseries', 1, 'NaN', 1, 'Inf', 1)
            md = checkfield(md, 'fieldname', 'smb.b_pos', 'timeseries', 1, 'NaN', 1, 'Inf', 1)
            md = checkfield(md, 'fieldname', 'smb.b_neg', 'timeseries', 1, 'NaN', 1, 'Inf', 1)
            md = checkfield(md, 'fieldname', 'smb.b_max', 'timeseries', 1, 'NaN', 1, 'Inf', 1)
            md = checkfield(md, 'fieldname', 'smb.b_min', 'timeseries', 1, 'NaN', 1, 'Inf', 1)

        md = checkfield(md, 'fieldname', 'smb.steps_per_step', '>=', 1, 'numel', [1])
        md = checkfield(md, 'fieldname', 'smb.averaging', 'numel', [1], 'values', [0, 1, 2])
        md = checkfield(md, 'fieldname', 'smb.requested_outputs', 'stringrow', 1)
        return md
    # }}}
    def marshall(self, prefix, md, fid):    # {{{

        yts = md.constants.yts

        WriteData(fid, prefix, 'name', 'md.smb.model', 'data', 9, 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'ela', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', yts)
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'b_pos', 'format', 'DoubleMat', 'mattype', 1, 'scale', 1. / yts, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', yts)
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'b_neg', 'format', 'DoubleMat', 'mattype', 1, 'scale', 1. / yts, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', yts)
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'b_max', 'format', 'DoubleMat', 'mattype', 1, 'scale', 1. / yts, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', yts)
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'b_min', 'format', 'DoubleMat', 'mattype', 1, 'scale', 1. / yts, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', yts)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'steps_per_step', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'averaging', 'format', 'Integer')

        #process requested outputs
        outputs = self.requested_outputs
        indices = [i for i, x in enumerate(outputs) if x == 'default']
        if len(indices) > 0:
            outputscopy = outputs[0:max(0, indices[0] - 1)] + self.defaultoutputs(md) + outputs[indices[0] + 1:]
            outputs = outputscopy
        WriteData(fid, prefix, 'data', outputs, 'name', 'md.smb.requested_outputs', 'format', 'StringArray')

    # }}}
