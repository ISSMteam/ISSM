import numpy as np

from checkfield import checkfield
from fielddisplay import fielddisplay
from WriteData import WriteData


class offlinesolidearthsolution(object):
    """OFFLINESOLIDEARTHSOLUTION class definition

    Usage:
        offlinesolidearthsolution = offlinesolidearthsolution()
    """

    def __init__(self, *args):  # {{{
        self.displacementeast = None
        self.displacementnorth = None
        self.displacementup = None
        self.geoid = None

        if len(args) == 0:
            self.setdefaultparameters()
        else:
            raise RuntimeError('constructor not supported')
    # }}}

    def __repr__(self):  # {{{
        s = '         units for time series is (yr)\n       external: offlinesolidearth solution\n'
        s += '{}\n'.format(fielddisplay(self, 'displacementeast', 'solid-Earth Eastwards bedrock displacement series (m)'))
        s += '{}\n'.format(fielddisplay(self, 'displacementnorth', 'solid-Earth Northwards bedrock displacement time series (m)'))
        s += '{}\n'.format(fielddisplay(self, 'displacementup', 'solid-Earth bedrock uplift time series (m)'))
        s += '{}\n'.format(fielddisplay(self, 'geoid', 'solid-Earth geoid time series (m)'))

        return s
    # }}}

    def setdefaultparameters(self):  # {{{
        self.displacementeast = []
        self.displacementnorth = []
        self.displacementup = []
        self.geoid = []
    # }}}

    def checkconsistency(self, md, solution, analyses):  # {{{
        if ('SealevelchangeAnalysis' not in analyses) or ((solution=='TransientSolution') and (md.solidearth.settings.isgrd==1)): 
            print('offlinesolidearthsolution checkconsistency error message: trying to run GRD patterns while supplying an offline solution for those patterns!')
            return md
        md = checkfield(md, 'fieldname', 'solidearth.external.displacementeast', 'Inf', 1, 'timeseries', 1)
        md = checkfield(md, 'fieldname', 'solidearth.external.displacementnorth', 'Inf', 1, 'timeseries', 1)
        md = checkfield(md, 'fieldname', 'solidearth.external.displacementup', 'Inf', 1, 'timeseries', 1)
        md = checkfield(md, 'fieldname', 'solidearth.external.geoid', 'Inf', 1, 'timeseries', 1)
    # }}}

    def marshall(self, prefix, md, fid):  # {{{
        yts = md.constants.yts

        # Transform our time series into time series rates
        if len(np.shape(self.displacementeast)) == 1:
            print('Warning: offlinesolidearthsolution.py::marshall: only one time step provided, assuming the values are rates per year')
            displacementeast_rate = np.append(np.array(self.displacementeast).reshape(-1, 1), 0)
            displacementnorth_rate = np.append(np.array(self.displacementnorth).reshape(-1, 1), 0)
            displacementup_rate = np.append(np.array(self.displacementup).reshape(-1, 1), 0)
            geoid_rate = np.append(np.array(self.geoid).reshape(-1, 1), 0)
        else:
            time = self.displacementeast[-1, :]
            dt = np.diff(time, axis=0)
            displacementeast_rate = np.diff(self.displacementeast[0:-2, :], 1, 1) / dt
            displacementeast_rate = np.append(displacementeast_rate,time[:-1].reshape(1,-1),axis=0)
            displacementnorth_rate = np.diff(self.displacementnorth[0:-2, :], 1, 1) / dt
            displacementnorth_rate = np.append(displacementnorth_rate,time[:-1].reshape(1,-1),axis=0)
            displacementup_rate = np.diff(self.displacementup[0:-2, :], 1, 1) / dt
            displacementup_rate = np.append(displacementup_rate,time[:-1].reshape(1,-1),axis=0)
            geoid_rate = np.diff(self.geoid[0:-2, :], 1, 1) / dt
            geoid_rate = np.append(geoid_rate,time[:-1].reshape(1,-1),axis=0)

        WriteData(fid, prefix, 'name', 'md.solidearth.external.nature', 'data', 2, 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'displacementeast', 'data', displacementeast_rate, 'format', 'DoubleMat', 'name', 'md.solidearth.external.displacementeast', 'mattype', 1, 'scale', 1 / yts,'timeserieslength', md.mesh.numberofvertices + 1, 'yts', yts);
        WriteData(fid, prefix, 'object', self, 'fieldname', 'displacementup', 'data', displacementup_rate,'format', 'DoubleMat', 'name', 'md.solidearth.external.displacementup', 'mattype', 1, 'scale', 1 / yts,'timeserieslength', md.mesh.numberofvertices + 1, 'yts', yts);
        WriteData(fid, prefix, 'object', self, 'fieldname', 'displacementnorth', 'data', displacementnorth_rate,'format', 'DoubleMat', 'name', 'md.solidearth.external.displacementnorth', 'mattype', 1, 'scale', 1 / yts,'timeserieslength', md.mesh.numberofvertices + 1, 'yts', yts);
        WriteData(fid, prefix, 'object', self, 'fieldname', 'geoid', 'data', geoid_rate,'format', 'DoubleMat', 'name', 'md.solidearth.external.geoid', 'mattype', 1, 'scale', 1 / yts,'timeserieslength', md.mesh.numberofvertices + 1, 'yts', yts);
    # }}}

    def extrude(self, md):  # {{{
        return self
    # }}}
