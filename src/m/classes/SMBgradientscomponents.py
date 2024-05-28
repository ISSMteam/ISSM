from fielddisplay import fielddisplay
from checkfield import checkfield
from WriteData import WriteData


class SMBgradientscomponents(object):
    """
    SMBgradients Class definition

       Usage:
          SMBgradients = SMBgradientscomponents();
    For now it has accumulation, runoff ans retention which could be aither refreezing and/or evaporation
    """

    def __init__(self):  # {{{
        self.accuref = float('NaN')
        self.accualti = float('NaN')
        self.accugrad = float('NaN')
        self.runoffref = float('NaN')
        self.runoffalti = float('NaN')
        self.runoffgrad = float('NaN')
        self.steps_per_step = 1
        self.averaging = 0
        self.requested_outputs = ['default']
    # }}}

    def __repr__(self):  # {{{
        string = "   surface forcings parameters:"
        string = "%s\n%s" % (string, fielddisplay(self, 'accuref', ' reference value of the accumulation'))
        string = "%s\n%s" % (string, fielddisplay(self, 'accualti', ' Altitude at which the accumulation is equal to the reference value'))
        string = "%s\n%s" % (string, fielddisplay(self, 'accugrad', ' Gradient of the variation of the accumulation (0 for uniform accumulation)'))
        string = "%s\n%s" % (string, fielddisplay(self, 'runoffref', ' reference value of the runoff m w.e. y-1 (temperature times ddf)'))
        string = "%s\n%s" % (string, fielddisplay(self, 'runoffalti', ' Altitude at which the runoff is equal to the reference value'))
        string = "%s\n%s" % (string, fielddisplay(self, 'runoffgrad', ' Gradient of the variation of the runoff (0 for uniform runoff) m w.e. m-1 y-1 (lapse rate times ddf)'))
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
        list = ['SmbMassBalance']
        if self.steps_per_step > 1:
            list.extend(['SmbMassBalanceSubstep'])
        return list
    # }}}

    def initialize(self, md):  # {{{
        #Nothing for now
        return self
    # }}}

    def checkconsistency(self, md, solution, analyses):  # {{{
        if 'MasstransportAnalysis' in analyses:
            md = checkfield(md, 'fieldname', 'smb.accualti', 'numel', [1], 'NaN', 1, 'Inf', 1)
            md = checkfield(md, 'fieldname', 'smb.accuref', 'singletimeseries', 1, 'NaN', 1, 'Inf', 1)
            md = checkfield(md, 'fieldname', 'smb.accugrad', 'singletimeseries', 1, 'NaN', 1, 'Inf', 1)
            md = checkfield(md, 'fieldname', 'smb.runoffalti', 'numel', [1], 'NaN', 1, 'Inf', 1)
            md = checkfield(md, 'fieldname', 'smb.runoffref', 'singletimeseries', 1, 'NaN', 1, 'Inf', 1)
            md = checkfield(md, 'fieldname', 'smb.runoffgrad', 'singletimeseries', 1, 'NaN', 1, 'Inf', 1)

        md = checkfield(md, 'fieldname', 'smb.steps_per_step', '>=', 1, 'numel', [1])
        md = checkfield(md, 'fieldname', 'smb.averaging', 'numel', [1], 'values', [0, 1, 2])
        md = checkfield(md, 'fieldname', 'masstransport.requested_outputs', 'stringrow', 1)
        return md
    # }}}

    def marshall(self, prefix, md, fid):    # {{{
        yts = md.constants.yts
        WriteData(fid, prefix, 'name', 'md.smb.model', 'data', 11, 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'accuref', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', 2, 'yts', yts, 'scale', 1. / yts)
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'accugrad', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', 2, 'yts', yts, 'scale', 1. / yts)
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'accualti', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'runoffref', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', 2, 'yts', yts, 'scale', 1. / yts)
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'runoffgrad', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', 2, 'yts', yts, 'scale', 1. / yts)
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'runoffalti', 'format', 'Double')
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
