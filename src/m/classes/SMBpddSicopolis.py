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

    def __repr__(self):  # {{{
        s = '   surface forcings parameters:\n'
        s += '   SICOPOLIS PDD scheme (Calov & Greve, 2005):\n'
        s += '{}\n'.format(fielddisplay(self, 'monthlytemperatures', 'monthly surface temperatures [K]'))
        s += '{}\n'.format(fielddisplay(self, 'precipitation', 'monthly surface precipitation [m/yr water eq]'))
        s += '{}\n'.format(fielddisplay(self, 'temperature_anomaly', 'anomaly to monthly reference temperature (additive [K])'))
        s += '{}\n'.format(fielddisplay(self, 'precipitation_anomaly', 'anomaly to monthly precipitation (multiplicative, e.g. q = q0*exp(0.070458*DeltaT) after Huybrechts (2002)) [no unit])'))
        s += '{}\n'.format(fielddisplay(self, 'smb_corr', 'correction of smb after PDD call [m/a]'))
        s += '{}\n'.format(fielddisplay(self, 's0p', 'should be set to elevation from precip source (between 0 and a few 1000s m, default is 0) [m]'))
        s += '{}\n'.format(fielddisplay(self, 's0t', 'should be set to elevation from temperature source (between 0 and a few 1000s m, default is 0) [m]'))
        s += '{}\n'.format(fielddisplay(self, 'rlaps', 'present day lapse rate (default is 7.4 degree/km)'))
        s += '{}\n'.format(fielddisplay(self, 'desfac', 'desertification elevation factor (default is -log(2.0)/1000)'))
        s += '{}\n'.format(fielddisplay(self, 'isfirnwarming', 'is firnwarming (Reeh 1991) activated (0 or 1, default is 1)'))
        s += '{}\n'.format(fielddisplay(self, 'steps_per_step', 'number of smb steps per time step'))
        s += '{}\n'.format(fielddisplay(self, 'averaging', 'averaging methods from short to long steps'))        
        s += '{}\n'.format(fielddisplay(self, 'pdd_fac_ice', 'Pdd factor for ice for all the domain [mm ice equiv/day/degree C]'))
        s += '{}\n'.format(fielddisplay(self, 'pdd_fac_snow', 'Pdd factor for snow for all the domain [mm ice equiv/day/degree C]'))
        s += '\t\t{}\n'.format('0: Arithmetic (default)')
        s += '\t\t{}\n'.format('1: Geometric')
        s += '\t\t{}\n'.format('2: Harmonic')
        s += '{}\n'.format(fielddisplay(self, 'requested_outputs', 'additional outputs requested (TemperaturePDD, SmbAccumulation, SmbMelt)'))
        return s
    # }}}

    def extrude(self, md):  # {{{
        self.precipitation = project3d(md, 'vector', self.precipitation, 'type', 'node')
        self.monthlytemperatures = project3d(md, 'vector', self.monthlytemperatures, 'type', 'node')
        self.temperature_anomaly = project3d(md, 'vector', self.temperature_anomaly, 'type', 'node')
        self.precipitation_anomaly = project3d(md, 'vector', self.precipitation_anomaly, 'type', 'node')
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

        if np.isnan(self.temperature_anomaly):
            self.temperature_anomaly = np.zeros((md.mesh.numberofvertices, ))
            print('      no SMBpddSicopolis.temperature_anomaly specified: values set as zero')

        if np.isnan(self.precipitation_anomaly):
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

    def checkconsistency(self, md, solution, analyses):  # {{{
        if solution == 'TransientSolution' and not md.transient.issmb:
            return
        if 'MasstransportAnalysis' in analyses:
            md = checkfield(md, 'fieldname', 'smb.desfac', '<=', 1, 'numel', 1)
            md = checkfield(md, 'fieldname', 'smb.s0p', '>=', 0, 'NaN', 1, 'Inf', 1, 'size', [md.mesh.numberofvertices, 1])
            md = checkfield(md, 'fieldname', 'smb.s0t', '>=', 0, 'NaN', 1, 'Inf', 1, 'size', [md.mesh.numberofvertices, 1])
            md = checkfield(md, 'fieldname', 'smb.rlaps', '>=', 0, 'numel', 1)
            md = checkfield(md, 'fieldname', 'smb.monthlytemperatures', 'NaN', 1, 'Inf', 1, 'size', [md.mesh.numberofvertices, 12])
            md = checkfield(md, 'fieldname', 'smb.precipitation', 'NaN', 1, 'Inf', 1, 'size', [md.mesh.numberofvertices, 12])
            md = checkfield(md, 'fieldname', 'smb.pdd_fac_ice', '>', 0, 'numel', 1)
            md = checkfield(md, 'fieldname', 'smb.pdd_fac_snow', '>', 0, 'numel', 1)
        md = checkfield(md, 'fieldname', 'smb.steps_per_step', '>=', 1, 'numel', [1])
        md = checkfield(md, 'fieldname', 'smb.averaging', 'numel', [1], 'values', [0, 1, 2])
        md = checkfield(md, 'fieldname', 'smb.requested_outputs', 'stringrow', 1)
        return md
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
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'temperature_anomaly', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', yts)
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
