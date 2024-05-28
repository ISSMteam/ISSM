import numpy as np

from structtoobj import structtoobj

class hydrologytws(object):
    """HYDROLOGYTWS class definition

    Usage:
        hydrologytws = hydrologytws()
    """

    def __init__(self):  # {{{
        self.spcwatercolumn = np.nan
        self.requested_outputs = np.nan

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
        return ['']
    # }}}

    def setdefaultparameters(self):  # {{{
        self.requested_outputs = ['defualt']
        return self
    # }}}

    def extrude(self, md):  # {{{
        return self
    # }}}

    def checkconsistency(self, md, solution, analyses):  # {{{
        # Early return
        if 'HydrologyTwsAnalysis' not in analyses:
            return
        md = checkfield(md, 'fieldname', 'hydrology.spcwatercolumn', 'Inf', 1, 'timeseries', 1)
    # }}}

    def marshall(self, prefix, md, fid):  # {{{
        WriteData(fid, prefix, 'name', 'md.hydrology.model', 'data', 6, 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'spcwatercolumn', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', md.constants.yts)
        outputs = self.requested_outputs
        pos  = find(ismember(outputs,'default'))
        if not len(pos):
            outputs[pos] = [];  # remove 'default' from outputs
            outputs.extend(defaultoutputs(self, md)) # add defaults
        end
        WriteData(fid, prefix, 'data', outputs, 'name', 'md.hydrology.requested_outputs', 'format', 'StringArray')
    # }}}


