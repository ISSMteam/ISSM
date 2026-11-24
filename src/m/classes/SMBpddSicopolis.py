import numpy as np

from checkfield import checkfield
from fielddisplay import fielddisplay
from helpers import *
from MatlabFuncs import *
from project3d import project3d
from WriteData import WriteData


class SMBpddSicopolis(object):
    """SMBPDDSICOPOLIS class definition

    Usage:
        SMBpddSicopolis = SMBpddSicopolis()

    temperature_anomaly / precipitation_anomaly can be either
        - an (nv x 1) array (same anomaly for all 12 months), or
        - a list/tuple of 12 arrays (one per calendar month), each (nv+1) x nt
          with the time row (in years) as the last row.
    """

    def __init__(self, *args):  # {{{
        self.precipitation = np.nan
        self.monthlytemperatures = np.nan
        self.temperature_anomaly = np.nan
        self.precipitation_anomaly = np.nan
        self.smb_corr = np.nan
        self.desfac = 0
        self.s0p = np.nan
        self.s0t = np.nan
        self.rlaps = 0
        self.isfirnwarming = 0
        self.pdd_fac_ice = 0
        self.pdd_fac_snow = 0
        self.steps_per_step = 1
        self.averaging = 0
        self.requested_outputs = []

        if len(args) == 0:
            self.setdefaultparameters()
        else:
            raise Exception('constructor not supported')
    # }}}

    @staticmethod
    def _ismonthly(field):  # {{{
        return isinstance(field, (list, tuple))
    # }}}

    @staticmethod
    def _isunset(field):  # {{{
        return (not isinstance(field, (list, tuple))) and np.all(np.isnan(field))
    # }}}

    def __repr__(self):  # {{{
        s = '   surface forcings parameters:\n'
        s += '   SICOPOLIS PDD scheme (Calov & Greve, 2005):\n'
        s += '{}\n'.format(fielddisplay(self, 'monthlytemperatures', 'monthly surface temperatures [K]'))
        s += '{}\n'.format(fielddisplay(self, 'precipitation', 'monthly surface precipitation [m/yr water eq]'))
        s += '{}\n'.format(fielddisplay(self, 'temperature_anomaly', 'anomaly to monthly reference temperature (additive [K]): nv x 1 (same for all months) or list of 12 (nv+1) x nt monthly time series (time in years, last row)'))
        s += '{}\n'.format(fielddisplay(self, 'precipitation_anomaly', 'anomaly to monthly precipitation (multiplicative, e.g. q = q0*exp(0.070458*DeltaT) after Huybrechts (2002)) [no unit]): nv x 1 or list of 12 (nv+1) x nt monthly time series'))
        s += '{}\n'.format(fielddisplay(self, 'smb_corr', 'correction of smb after PDD call [m/a]'))
        s += '{}\n'.format(fielddisplay(self, 's0p', 'should be set to elevation from precip source (between 0 and a few 1000s m, default is 0) [m]'))
        s += '{}\n'.format(fielddisplay(self, 's0t', 'should be set to elevation from temperature source (between 0 and a few 1000s m, default is 0) [m]'))
        s += '{}\n'.format(fielddisplay(self, 'rlaps', 'present day lapse rate (default is 7.4 degree/km)'))
        s += '{}\n'.format(fielddisplay(self, 'desfac', 'desertification elevation factor (default is -log(2.0)/1000)'))
        s += '{}\n'.format(fielddisplay(self, 'isfirnwarming', 'is firnwarming (Reeh 1991) activated (0 or 1, default is 1)'))
        s += '{}\n'.format(fielddisplay(self, 'steps_per_step', 'number of smb steps per time step'))
        s += '{}\n'.format(fielddisplay(self, 'averaging', 'averaging methods from short to long steps'))        
        s += '{}\n'.format(fielddisplay(self, 'pdd_fac_ice', 'Pdd factor for ice for all the domain [mm water equiv/day/degree C]'))
        s += '{}\n'.format(fielddisplay(self, 'pdd_fac_snow', 'Pdd factor for snow for all the domain [mm water equiv/day/degree C]'))
        s += '\t\t{}\n'.format('0: Arithmetic (default)')
        s += '\t\t{}\n'.format('1: Geometric')
        s += '\t\t{}\n'.format('2: Harmonic')
        s += '{}\n'.format(fielddisplay(self, 'requested_outputs', 'additional outputs requested (TemperaturePDD, SmbAccumulation, SmbMelt, SmbRunoff)'))
        return s
    # }}}

    def extrude(self, md):  # {{{
        self.precipitation = project3d(md, 'vector', self.precipitation, 'type', 'node')
        self.monthlytemperatures = project3d(md, 'vector', self.monthlytemperatures, 'type', 'node')
        for name in ['temperature_anomaly', 'precipitation_anomaly']:
            value = getattr(self, name)
            if self._ismonthly(value):
                setattr(self, name, [project3d(md, 'vector', v, 'type', 'node') for v in value])
            else:
                setattr(self, name, project3d(md, 'vector', value, 'type', 'node'))
        self.smb_corr = project3d(md, 'vector', self.smb_corr, 'type', 'node')
        self.s0p = project3d(md, 'vector', self.s0p, 'type', 'node')
        self.s0t = project3d(md, 'vector', self.s0t, 'type', 'node')
    # }}}

    def defaultoutputs(self, md):  # {{{
        listing = ['SmbMassBalance']
        return listing
    # }}}

    def initialize(self, md):  # {{{
        if np.isnan(self.s0p):
            self.s0p = np.zeros((md.mesh.numberofvertices, ))
            print('      no SMBpddSicopolis.s0p specified: values set as zero')

        if np.isnan(self.s0t):
            self.s0t = np.zeros((md.mesh.numberofvertices, ))
            print('      no SMBpddSicopolis.s0t specified: values set as zero')

        if self._isunset(self.temperature_anomaly):
            self.temperature_anomaly = np.zeros((md.mesh.numberofvertices, ))
            print('      no SMBpddSicopolis.temperature_anomaly specified: values set as zero')

        if self._isunset(self.precipitation_anomaly):
            self.precipitation_anomaly = np.ones((md.mesh.numberofvertices, ))
            print('      no SMBpddSicopolis.precipitation_anomaly specified: values set as ones')

        if np.isnan(self.smb_corr):
            self.smb_corr = np.zeros((md.mesh.numberofvertices, ))
            print('      no SMBpddSicopolis.smb_corr specified: values set as zero')
        return self
    # }}}

    def setdefaultparameters(self):  # {{{
        self.isfirnwarming = 1
        self.desfac = -np.log(2.0) / 1000
        self.rlaps = 7.4
        self.pdd_fac_ice = 7.28
        self.pdd_fac_snow = 2.73
        self.requested_outputs = ['default']
        return self
    # }}}

    def _checkmonthly(self, md, name, *extra):  # {{{
        """Check a monthly anomaly: list of 12 (nv+1) x nt time series"""
        value = getattr(self, name)
        if len(value) != 12:
            md = md.checkmessage('smb.{} should be a list of 12 monthly time series (got {})'.format(name, len(value)))
            return md
        for i in range(12):
            md = checkfield(md, 'field', np.asarray(value[i]), 'fieldname', 'smb.{}[{}]'.format(name, i),
                            'NaN', 1, 'Inf', 1, 'timeseries', 1, 'size', [md.mesh.numberofvertices + 1, np.nan], *extra)
        return md
    # }}}

    def checkconsistency(self, md, solution, analyses):  # {{{
        if solution == 'TransientSolution' and not md.transient.issmb:
            return md
        if 'MasstransportAnalysis' in analyses:
            md = checkfield(md, 'fieldname', 'smb.desfac', '<=', 1, 'numel', 1)
            md = checkfield(md, 'fieldname', 'smb.s0p', '>=', 0, 'NaN', 1, 'Inf', 1, 'size', [md.mesh.numberofvertices, 1])
            md = checkfield(md, 'fieldname', 'smb.s0t', '>=', 0, 'NaN', 1, 'Inf', 1, 'size', [md.mesh.numberofvertices, 1])
            md = checkfield(md, 'fieldname', 'smb.rlaps', '>=', 0, 'numel', 1)
            md = checkfield(md, 'fieldname', 'smb.monthlytemperatures', 'NaN', 1, 'Inf', 1, 'size', [md.mesh.numberofvertices, 12], '>' ,0., '<', 300.)
            md = checkfield(md, 'fieldname', 'smb.precipitation', 'NaN', 1, 'Inf', 1, 'size', [md.mesh.numberofvertices, 12],'>=' ,0.)
            md = checkfield(md, 'fieldname', 'smb.pdd_fac_ice', '>', 0, 'numel', 1)
            md = checkfield(md, 'fieldname', 'smb.pdd_fac_snow', '>', 0, 'numel', 1)
            if self._ismonthly(self.temperature_anomaly):
                md = self._checkmonthly(md, 'temperature_anomaly')
            else:
                md = checkfield(md, 'fieldname', 'smb.temperature_anomaly', 'NaN', 1, 'Inf', 1)
            if self._ismonthly(self.precipitation_anomaly):
                md = self._checkmonthly(md, 'precipitation_anomaly', '>=', 0.)
            else:
                md = checkfield(md, 'fieldname', 'smb.precipitation_anomaly', 'NaN', 1, 'Inf', 1)
        md = checkfield(md, 'fieldname', 'smb.steps_per_step', '>=', 1, 'numel', [1])
        md = checkfield(md, 'fieldname', 'smb.averaging', 'numel', [1], 'values', [0, 1, 2])
        md = checkfield(md, 'fieldname', 'smb.requested_outputs', 'stringrow', 1)
        return md
    # }}}

    def _marshall_monthly(self, fid, prefix, md, name, scale):  # {{{
        """Write a list of 12 (nv+1) x nt matrices as a MatArray.

        WriteData applies scale / timeserieslength / yts to each matrix of a MatArray.
        We hand it copies: WriteData only makes a shallow copy of a list, so it would
        otherwise rescale the arrays stored in the model in place."""
        records = []
        for v in getattr(self, name):
            m = np.array(v, dtype=float, copy=True)
            if m.ndim == 1:
                m = m.reshape(-1, 1)
            records.append(m)
        WriteData(fid, prefix, 'data', records, 'name', 'md.smb.' + name, 'format', 'MatArray', 'scale', scale, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', md.constants.yts)
    # }}}

    def marshall(self, prefix, md, fid):  # {{{
        yts = md.constants.yts
        WriteData(fid, prefix, 'name', 'md.smb.model', 'data', 10, 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'isfirnwarming', 'format', 'Boolean')
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'desfac', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 's0p', 'format', 'DoubleMat', 'mattype', 1)
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 's0t', 'format', 'DoubleMat', 'mattype', 1)
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'rlaps', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'pdd_fac_ice', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'pdd_fac_snow', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'monthlytemperatures', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', yts)
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'precipitation', 'format', 'DoubleMat', 'mattype', 1, 'scale', 1. / yts, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', yts)
        if self._ismonthly(self.temperature_anomaly):
            self._marshall_monthly(fid, prefix, md, 'temperature_anomaly', 1.)
        else:
            WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'temperature_anomaly', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', yts)
        if self._ismonthly(self.precipitation_anomaly):
            self._marshall_monthly(fid, prefix, md, 'precipitation_anomaly', 1. / yts)
        else:
            WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'precipitation_anomaly', 'format', 'DoubleMat', 'mattype', 1, 'scale', 1. / yts, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', yts)
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'smb_corr', 'format', 'DoubleMat', 'mattype', 1, 'scale', 1. / yts, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', yts)
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
