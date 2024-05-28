from fielddisplay import fielddisplay
from MatlabFuncs import *
from model import *
from checkfield import checkfield
from WriteData import WriteData


class esa(object):
    """
    ESA class definition

        Usage:
          esa = esa();
    """

    def __init__(self):  # {{{
        self.deltathickness = float('NaN')
        self.love_h = 0  #provided by PREM model()
        self.love_l = 0  #ideam
        self.hemisphere = 0
        self.degacc = 0
        self.requested_outputs = []
        self.transitions = []

    #set defaults
        self.setdefaultparameters()
    # }}}

    def __repr__(self):  # {{{
        string = '   esa parameters:'
        string = "%s\n%s" % (string, fielddisplay(self, 'deltathickness', 'thickness change: ice height equivalent [m]'))
        string = "%s\n%s" % (string, fielddisplay(self, 'love_h', 'load Love number for radial displacement'))
        string = "%s\n%s" % (string, fielddisplay(self, 'love_l', 'load Love number for horizontal displaements'))
        string = "%s\n%s" % (string, fielddisplay(self, 'hemisphere', 'North-south, East-west components of 2-D horiz displacement vector:-1 south, 1 north'))
        string = "%s\n%s" % (string, fielddisplay(self, 'degacc', 'accuracy (default .01 deg) for numerical discretization of the Green''s functions'))
        string = "%s\n%s" % (string, fielddisplay(self, 'transitions', 'indices into parts of the mesh that will be icecaps'))
        string = "%s\n%s" % (string, fielddisplay(self, 'requested_outputs', 'additional outputs requested (default: EsaUmotion)'))

        return string
    # }}}

    def setdefaultparameters(self):  # {{{
        #numerical discretization accuracy
        self.degacc = 0.01
        #computational flags:
        self.hemisphere = 0
        #output default:
        self.requested_outputs = ['default']
        #transitions should be a cell array of vectors:
        self.transitions = []
        #default output
        self.requested_outputs = ['default']
        return self
    # }}}

    def checkconsistency(self, md, solution, analyses):  # {{{
        #Early return
        if (solution != 'EsaAnalysis'):
            return md

        md = checkfield(md, 'fieldname', 'esa.deltathickness', 'NaN', 1, 'Inf', 1, 'size', [md.mesh.numberofelements, 1])
        md = checkfield(md, 'fieldname', 'esa.love_h', 'NaN', 1, 'Inf', 1)
        md = checkfield(md, 'fieldname', 'esa.love_l', 'NaN', 1, 'Inf', 1)
        md = checkfield(md, 'fieldname', 'esa.hemisphere', 'NaN', 1, 'Inf', 1)
        md = checkfield(md, 'fieldname', 'esa.degacc', 'size', [1, 1], '>=', 1e-10)
        md = checkfield(md, 'fieldname', 'esa.requested_outputs', 'stringrow', 1)

    #check that love numbers are provided at the same level of accuracy:
        if (size(self.love_h, 0) != size(self.love_l, 0)):
            error('esa error message: love numbers should be provided at the same level of accuracy')
        return md
    # }}}

    def defaultoutputs(self, md):  # {{{
        return ['EsaUmotion']
    # }}}

    def marshall(self, prefix, md, fid):  # {{{
        WriteData(fid, prefix, 'object', self, 'fieldname', 'deltathickness', 'format', 'DoubleMat', 'mattype', 2)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'love_h', 'format', 'DoubleMat', 'mattype', 1)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'love_l', 'format', 'DoubleMat', 'mattype', 1)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'hemisphere', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'degacc', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'transitions', 'format', 'MatArray')

    #process requested outputs
        outputs = self.requested_outputs
        indices = [i for i, x in enumerate(outputs) if x == 'default']
        if len(indices) > 0:
            outputscopy = outputs[0:max(0, indices[0] - 1)] + self.defaultoutputs(md) + outputs[indices[0] + 1:]
            outputs = outputscopy
        WriteData(fid, prefix, 'data', outputs, 'name', 'md.esa.requested_outputs', 'format', 'StringArray')
    # }}}
