from fielddisplay import fielddisplay
from checkfield import checkfield
from WriteData import WriteData
from project3d import project3d


class calvinglevermann(object):
    """
    CALVINGLEVERMANN class definition

       Usage:
          calvinglevermann = calvinglevermann()
    """

    def __init__(self):  # {{{

        self.coeff = float('NaN')

    #set defaults
        self.setdefaultparameters()

    # }}}

    def __repr__(self):  # {{{
        string = '   Calving Levermann parameters:'
        string = "%s\n%s" % (string, fielddisplay(self, 'coeff', 'proportionality coefficient in Levermann model'))

        return string
    # }}}

    def extrude(self, md):  # {{{
        self.coeff = project3d(md, 'vector', self.coeff, 'type', 'node')
        return self
    # }}}

    def setdefaultparameters(self):  # {{{
        #Proportionality coefficient in Levermann model
        self.coeff = 2e13
    # }}}

    def checkconsistency(self, md, solution, analyses):  # {{{
        #Early return
        if (solution != 'TransientSolution') or (not md.transient.ismovingfront):
            return md

        md = checkfield(md, 'fieldname', 'calving.coeff', 'size', [md.mesh.numberofvertices], '>', 0)
        return md
    # }}}

    def marshall(self, prefix, md, fid):  # {{{
        WriteData(fid, prefix, 'name', 'md.calving.law', 'data', 3, 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'coeff', 'format', 'DoubleMat', 'mattype', 1)
    # }}}
