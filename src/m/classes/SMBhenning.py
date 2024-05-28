import numpy as np

from checkfield import checkfield
from fielddisplay import fielddisplay
from project3d import project3d
from WriteData import WriteData


class SMBforcing(object):
    """SMBhenning class definition

    Usage:
        SMBhenning = SMBhenning()
    """

    def __init__(self, *args):  # {{{
        self.smbref = np.nan
        self.steps_per_step = 1
        self.averaging = 0
        self.requested_outputs = []

        nargs = len(args)
        if nargs == 0:
            self.setdefaultparameters()
        elif nargs == 1:
            # TODO: Replace the following with constructor
            self.setdefaultparameters()
        else:
            raise Exception('constructor not supported')
    # }}}

    def __repr__(self):  # {{{
        s = '   surface forcings parameters:\n'
        s += '{}\n'.format(fielddisplay(self, 'mass_balance', 'surface mass balance [m/yr ice eq]'))
        s += '{}\n'.format(fielddisplay(self, 'steps_per_step', 'number of smb steps per time step'))
        s += '{}\n'.format(fielddisplay(self, 'averaging', 'averaging methods from short to long steps'))
        s += '\t\t{}\n'.format('0: Arithmetic (default)')
        s += '\t\t{}\n'.format('1: Geometric')
        s += '\t\t{}\n'.format('2: Harmonic')
        s += '{}\n'.format(fielddisplay(self, 'requested_outputs', 'additional outputs requested'))
        return s
    # }}}

    def extrude(self, md):  # {{{
        self.smbref = project3d(md, 'vector', self.smbref, 'type', 'node')
        return self
    # }}}

    def defaultoutputs(self, md):  # {{{
        return ['SmbMassBalance']
    # }}}

    def initialize(self, md):  # {{{
        if np.all(np.isnan(self.smbref)):
            self.smbref = np.zeros((md.mesh.numberofvertices))
            print("      no smb.smbref specified: values set as zero")
        return self
    # }}}

    def checkconsistency(self, md, solution, analyses):  # {{{
        if solution == 'TransientSolution' and not md.transient.issmb:
            return
        if 'MasstransportAnalysis' in analyses:
            md = checkfield(md, 'fieldname', 'smb.mass_balance', 'timeseries', 1, 'NaN', 1, 'Inf', 1)
        if 'BalancethicknessAnalysis' in analyses:
            md = checkfield(md, 'fieldname', 'smb.mass_balance', 'size', [md.mesh.numberofvertices], 'NaN', 1, 'Inf', 1)
        md = checkfield(md, 'fieldname', 'smb.steps_per_step', '>=', 1, 'numel', [1])
        md = checkfield(md, 'fieldname', 'smb.averaging', 'numel', [1], 'values', [0, 1, 2])
        md = checkfield(md, 'fieldname', 'smb.requested_outputs', 'stringrow', 1)
        return md
    # }}}

    def marshall(self, prefix, md, fid):  # {{{
        yts = md.constants.yts
        WriteData(fid, prefix, 'name', 'md.smb.model', 'data', 1, 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'mass_balance', 'format', 'DoubleMat', 'mattype', 1, 'scale', 1 / yts, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', yts)
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

    def setdefaultparameters(self):  #{{{
        self.requested_outputs = ['default']
        return self
    # }}}
