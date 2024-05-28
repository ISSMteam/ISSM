import numpy as np
from fielddisplay import fielddisplay
from checkfield import checkfield
from WriteData import WriteData


class hydrologypism(object):
    """
    HydrologyPism class definition

    Usage:
      hydrologypism = hydrologypism()
    """

    def __init__(self):  # {{{
        self.drainage_rate = np.nan
        self.watercolumn_max = np.nan

        # set defaults
        self.setdefaultparameters()
    # }}}

    def extrude(self, md):  # {{{
        return self
    # }}}

    def defaultoutputs(self, md):  # {{{
        list = ['Watercolumn']
        return list

    def setdefaultparameters(self):  # {{{
        return self

    def checkconsistency(self, md, solution, analyses):  # {{{

        #Early return
        if 'HydrologyPismAnalysis' not in analyses:
            return

        md = checkfield(md, 'fieldname', 'hydrology.drainage_rate', 'Inf', 1, 'NaN', 1, '>=', 0, 'size', [md.mesh.numberofvertices])
        md = checkfield(md, 'fieldname', 'friction.watercolumn_max', 'NaN', 1, 'Inf', 1, '>', 0., 'size', [md.mesh.numberofvertices])

    def __repr__(self):  # {{{
        string = '   hydrologypism solution parameters:'
        string = "%s\n%s" % (string, fielddisplay(self, 'drainage_rate', 'fixed drainage rate [mm / yr]'))
        string = "%s\n%s" % (string, fielddisplay(self, 'watercolumn_max', 'maximum water column height [m], recommended default: 2 m'))
        return string
        # }}}

    def marshall(self, prefix, md, fid):  # {{{
        yts = md.constants.yts

        WriteData(fid, prefix, 'name', 'md.hydrology.model', 'data', 4, 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'class', 'hydrology', 'fieldname', 'drainage_rate', 'format', 'DoubleMat', 'mattype', 1, 'scale', 1. / (1000. * yts))  #from mm / yr to m / s
        WriteData(fid, prefix, 'class', 'hydrology', 'object', self, 'fieldname', 'watercolumn_max', 'format', 'DoubleMat', 'mattype', 1)
        WriteData(fid, prefix, 'data', {'Watercolumn'}, 'name', 'md.hydrology.requested_outputs', 'format', 'StringArray')
    # }}}
