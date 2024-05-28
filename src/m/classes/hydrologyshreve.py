import numpy as np

from fielddisplay import fielddisplay
from checkfield import checkfield
from WriteData import WriteData


class hydrologyshreve(object):
    """HYDROLOGYSHREVE class definition

    Usage:
        hydrologyshreve = hydrologyshreve()
    """

    def __init__(self, *args):  # {{{
        self.spcwatercolumn = np.nan
        self.stabilization = 0
        self.requested_outputs = []

        nargs = len(args)
        if nargs == 0:
            self.setdefaultparameters()
        elif nargs == 1:
            self.setdefaultparameters(args)
        else:
            raise RuntimeError('constructor not supported')
    # }}}

    def __repr__(self):  # {{{
        s = '   hydrologyshreve solution parameters:\n'
        s += '{}\n'.format(fielddisplay(self, 'spcwatercolumn', 'water thickness constraints (NaN means no constraint) [m]'))
        s += '{}\n'.format(fielddisplay(self, 'stabilization', 'artificial diffusivity (default: 1). can be more than 1 to increase diffusivity.'))
        s += '{}\n'.format(fielddisplay(self, 'requested_outputs', 'additional outputs requested'))
        return s
    # }}}

    def setdefaultparameters(self):  # {{{
        # Type of stabilization to use 0:nothing 1:artificial_diffusivity
        self.stabilization = 1
        self.requested_outputs = ['default']
    # }}}

    def checkconsistency(self, md, solution, analyses):  # {{{
        #Early return
        if 'HydrologyShreveAnalysis' not in analyses or (solution == 'TransientSolution' and not md.transient.ishydrology):
            return md

        md = checkfield(md, 'fieldname', 'hydrology.spcwatercolumn', 'Inf', 1, 'timeseries', 1)
        md = checkfield(md, 'fieldname', 'hydrology.stabilization', '>=', 0)
        return md
    # }}}

    def defaultoutputs(self, md):  # {{{
        return ['Watercolumn', 'HydrologyWaterVx', 'HydrologyWaterVy']
    # }}}

    def marshall(self, prefix, md, fid):  # {{{
        WriteData(fid, prefix, 'name', 'md.hydrology.model', 'data', 2, 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'spcwatercolumn', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', md.constants.yts)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'stabilization', 'format', 'Double')
        outputs = self.requested_outputs
        indices = [i for i, x in enumerate(outputs) if x == 'default']
        if len(indices):
            outputscopy = outputs[0:max(0, indices[0] - 1)] + self.defaultoutputs(md) + outputs[indices[0] + 1:]
            outputs = outputscopy
        WriteData(fid, prefix, 'data', outputs, 'name', 'md.hydrology.requested_outputs', 'format', 'StringArray')
    # }}}

    def extrude(self, md):  # {{{
        return self
    # }}}
