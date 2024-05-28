import numpy as np
from fielddisplay import fielddisplay
from checkfield import checkfield
from WriteData import WriteData


class SMBgradients(object):
    """SMBgradients Class definition

    Usage:
        SMBgradients = SMBgradients();
    """

    def __init__(self):  # {{{
        self.href = np.nan
        self.smbref = np.nan
        self.b_pos = np.nan
        self.b_neg = np.nan
        self.steps_per_step = 1
        self.averaging = 0
        self.requested_outputs = []

        # Set defaults
        self.setdefaultparameters()
    # }}}

    def __repr__(self):  # {{{
        s = '   surface forcings parameters:\n'
        s += '{}\n'.format(fielddisplay(self, 'issmbgradients', 'is smb gradients method activated (0 or 1, default is 0)'))
        s += '{}\n'.format(fielddisplay(self, 'href', 'reference elevation from which deviation is used to calculate SMB adjustment in smb gradients method'))
        s += '{}\n'.format(fielddisplay(self, 'smbref', 'reference smb from which deviation is calculated in smb gradients method [m/yr ice equiv]'))
        s += '{}\n'.format(fielddisplay(self, 'b_pos', 'slope of hs - smb regression line for accumulation regime required if smb gradients is activated'))
        s += '{}\n'.format(fielddisplay(self, 'b_neg', 'slope of hs - smb regression line for ablation regime required if smb gradients is activated'))
        s += '{}\n'.format(fielddisplay(self, 'steps_per_step', 'number of smb steps per time step'))
        s += '{}\n'.format(fielddisplay(self, 'averaging', 'averaging methods from short to long steps'))
        s += '\t\t{}\n'.format('0: Arithmetic (default)')
        s += '\t\t{}\n'.format('1: Geometric')
        s += '\t\t{}\n'.format('2: Harmonic')
        s += '{}\n'.format(fielddisplay(self, 'requested_outputs', 'additional outputs requested'))

        return s
    # }}}

    def extrude(self, md):  # {{{
        # Nothing for now
        return self
    # }}}

    def defaultoutputs(self, md):  # {{{
        return ['SmbMassBalance']
    # }}}

    def setdefaultparameters(self):  # {{{
        # Output default
        self.requested_outputs = ['default']
        return self
    # }}}

    def initialize(self, md):  # {{{
        # Nothing for now
        return self
    # }}}

    def checkconsistency(self, md, solution, analyses):  # {{{
        if 'MasstransportAnalysis' in analyses:
            md = checkfield(md, 'fieldname', 'smb.href', 'timeseries', 1, 'NaN', 1, 'Inf', 1)
            md = checkfield(md, 'fieldname', 'smb.smbref', 'timeseries', 1, 'NaN', 1, 'Inf', 1)
            if np.max(np.max(np.abs(md.smb.smbref[0:-1,]))) < 1:
                print('!!! Warning: SMBgradients now expects smbref to be in m/yr ice eq. instead of mm/yr water eq.')
            md = checkfield(md, 'fieldname', 'smb.b_pos', 'timeseries', 1, 'NaN', 1, 'Inf', 1)
            md = checkfield(md, 'fieldname', 'smb.b_neg', 'timeseries', 1, 'NaN', 1, 'Inf', 1)

        md = checkfield(md, 'fieldname', 'smb.steps_per_step', '>=', 1, 'numel', [1])
        md = checkfield(md, 'fieldname', 'smb.averaging', 'numel', [1], 'values', [0, 1, 2])
        md = checkfield(md, 'fieldname', 'masstransport.requested_outputs', 'stringrow', 1)
        return md
    # }}}

    def marshall(self, prefix, md, fid):  # {{{
        yts = md.constants.yts

        WriteData(fid, prefix, 'name', 'md.smb.model', 'data', 6, 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'href', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', yts)
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'smbref', 'format', 'DoubleMat', 'mattype', 1, 'scale', 1. / yts, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', yts)
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'b_pos', 'format', 'DoubleMat', 'mattype', 1, 'scale', 1. / yts, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', yts)
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'b_neg', 'format', 'DoubleMat', 'mattype', 1, 'scale', 1. / yts, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', yts)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'steps_per_step', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'averaging', 'format', 'Integer')

        # Process requested outputs
        outputs = self.requested_outputs
        indices = [i for i, x in enumerate(outputs) if x == 'default']
        if len(indices) > 0:
            outputscopy = outputs[0:max(0, indices[0] - 1)] + self.defaultoutputs(md) + outputs[indices[0] + 1:]
            outputs = outputscopy
        WriteData(fid, prefix, 'data', outputs, 'name', 'md.smb.requested_outputs', 'format', 'StringArray')
    # }}}
