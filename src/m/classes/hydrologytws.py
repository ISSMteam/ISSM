import numpy as np

from checkfield import checkfield
from fielddisplay import fielddisplay
from structtoobj import structtoobj
from WriteData import WriteData

class hydrologytws(object):
    """HYDROLOGYTWS class definition

    Usage:
        hydrologytws = hydrologytws()
    """

    def __init__(self, *args):  # {{{
        self.spcwatercolumn = np.nan
        self.requested_outputs = []

        nargs = len(args)
        if nargs == 0:
            self.setdefaultparameters()
        elif nargs == 1:
            self = structtoobj(self, args[0])
        else:
            raise RuntimeError('constructor not supported')
    # }}}

    def __repr__(self):  # {{{
        s = '   hydrologytws solution parameters:\n'
        s += '{}\n'.format(fielddisplay(self, 'spcwatercolumn', 'water thickness constraints (NaN means no constraint) [m]'))
        s += '{}\n'.format(fielddisplay(self, 'requested_outputs', 'additional outputs requested'))
        return s
    # }}}

    def defaultoutputs(self, md):  # {{{
        return ['Watercolumn']
    # }}}

    def setdefaultparameters(self):  # {{{
        self.requested_outputs = ['default']
        return self
    # }}}

    def extrude(self, md):  # {{{
        return self
    # }}}

    def checkconsistency(self, md, solution, analyses):  # {{{
        # Early return
        if 'HydrologyTwsAnalysis' not in analyses:
            return md
        md = checkfield(md, 'fieldname', 'hydrology.spcwatercolumn', 'Inf', 1, 'timeseries', 1)
        return md
    # }}}

    def marshall(self, prefix, md, fid):  # {{{
        WriteData(fid, prefix, 'name', 'md.hydrology.model', 'data', 6, 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'spcwatercolumn', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', md.constants.yts)
        outputs = self.requested_outputs
        indices = [i for i, x in enumerate(outputs) if x == 'default']
        if len(indices):
            outputscopy = outputs[:indices[0]] + self.defaultoutputs(md) + outputs[indices[0] + 1:]
            outputs = outputscopy
        WriteData(fid, prefix, 'data', outputs, 'name', 'md.hydrology.requested_outputs', 'format', 'StringArray')
    # }}}
