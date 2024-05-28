import numpy as np
from fielddisplay import fielddisplay
from checkfield import checkfield
from WriteData import WriteData
from project3d import project3d


class SMBd18opdd(object):
    """
    SMBd18opdd Class definition

       Usage:
          SMBd18opdd = SMBd18opdd()
    """
    def __init__(self, *args):  # {{{
        self.desfac = 0.
        self.s0p = float('NaN')
        self.s0t = float('NaN')
        self.rlaps = 0.
        self.rlapslgm = 0.
        self.dpermil = 0.
        self.f = 0.
        self.Tdiff = float('NaN')
        self.sealev = float('NaN')
        self.ismungsm = 0
        self.isd18opd = 0
        self.issetpddfac = 0
        self.istemperaturescaled = 0
        self.isprecipscaled = 0
        self.delta18o = float('NaN')
        self.delta18o_surface = float('NaN')
        self.temperatures_presentday = float('NaN')
        self.precipitations_presentday = float('NaN')
        self.temperatures_reconstructed = float('NaN')
        self.precipitations_reconstructed = float('NaN')
        self.pddfac_snow = float('NaN')
        self.pddfac_ice = float('NaN')
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
        s += '{}\n'.format(fielddisplay(self, 'isd18opd', 'is delta18o parametrisation from present day temperature and precipitation activated (0 or 1, default is 0)'))
        s += '{}\n'.format(fielddisplay(self, 'issetpddfac', 'is user passing in defined pdd factors (0 or 1, default is 0)'))
        s += '{}\n'.format(fielddisplay(self, 'desfac', 'desertification elevation factor (between 0 and 1, default is 0.5) [m]'))
        s += '{}\n'.format(ielddisplay(self, 's0p', 'should be set to elevation from precip source (between 0 and a few 1000s m, default is 0) [m]'))
        s += '{}\n'.format(fielddisplay(self, 's0t', 'should be set to elevation from temperature source (between 0 and a few 1000s m, default is 0) [m]'))
        s += '{}\n'.format(fielddisplay(self, 'rlaps', 'present day lapse rate [degree/km]'))

        if self.isd18opd:
            s += '{}\n'.format(fielddisplay(self, 'temperatures_presentday', 'monthly present day surface temperatures [K], required if delta18o/mungsm is activated'))
            s += '{}\n'.format(fielddisplay(self, 'precipitations_presentday', 'monthly surface precipitation [m/yr water eq], required if delta18o or mungsm is activated'))
            s += '{}\n'.format(fielddisplay(self, 'istemperaturescaled', 'if delta18o parametrisation from present day temperature and precipitation is activated, is temperature scaled to delta18o value? (0 or 1, default is 1)'))
            s += '{}\n'.format(fielddisplay(self, 'isprecipscaled', 'if delta18o parametrisation from present day temperature and precipitation is activated, is precipitation scaled to delta18o value? (0 or 1, default is 1)'))

            if self.istemperaturescaled == 0:
                s += '{}\n'.format(fielddisplay(self, 'temperatures_reconstructed', 'monthly historical surface temperatures [K], required if delta18o/mungsm/d18opd is activated and istemperaturescaled is not activated'))

            if self.isprecipscaled == 0:
                s += '{}\n'.format(fielddisplay(self, 'precipitations_reconstructed', 'monthly historical precipitation [m/yr water eq], required if delta18o/mungsm/d18opd is activated and isprecipscaled is not activated'))

            s += '{}\n'.format(fielddisplay(self, 'delta18o', 'delta18o [per mil], required if pdd is activated and delta18o activated'))
            s += '{}\n'.format(fielddisplay(self, 'dpermil', 'degree per mil, required if d18opd is activated'))
            s += '{}\n'.format(fielddisplay(self, 'f', 'precip/temperature scaling factor, required if d18opd is activated'))

        if self.issetpddfac == 1:
            s += '{}\n'.format(fielddisplay(self, 'pddfac_snow', 'Pdd factor for snow for all the domain [mm ice equiv/day/degree C]'))
            s += '{}\n'.format(fielddisplay(self, 'pddfac_ice', 'Pdd factor for ice for all the domain [mm ice equiv/day/degree C]'))

        s += '{}\n'.format(fielddisplay(self, 'steps_per_step', 'number of smb steps per time step'))
        s += '{}\n'.format(fielddisplay(self, 'averaging', 'averaging methods from short to long steps'))
        s += '\t\t{}\n'.format('0: Arithmetic (default)')
        s += '\t\t{}\n'.format('1: Geometric')
        s += '\t\t{}\n'.format('2: Harmonic')
        s += '{}\n'.format(fielddisplay(self, 'requested_outputs', 'additional outputs requested'))
        return s
    # }}}
    def extrude(self, md):  # {{{
        if self.isd18opd:
            self.temperatures_presentday = project3d(md, 'vector', self.temperatures_presentday, 'type', 'node')
        if self.isd18opd:
            self.precipitations_presentday = project3d(md, 'vector', self.precipitations_presentday, 'type', 'node')
        if self.istemperaturescaled == 0:
            self.temperatures_reconstructed = project3d(md, 'vector', self.temperatures_reconstructed, 'type', 'node')
        if self.isprecipscaled == 0:
            self.precipitations_reconstructed = project3d(md, 'vector', self.precipitations_reconstructed, 'type', 'node')
        self.s0p = project3d(md, 'vector', self.s0p, 'type', 'node')
        self.s0t = project3d(md, 'vector', self.s0t, 'type', 'node')
        return self
    # }}}
    def defaultoutputs(self, md):  # {{{
        return ['SmbMassBalance']
    # }}}
    def initialize(self, md):  # {{{
        if np.all(np.isnan(self.s0p)):
            self.s0p = np.zeros((md.mesh.numberofvertices))
            print("      no SMBd18opdd.s0p specified: values set as zero")

        if np.all(np.isnan(self.s0t)):
            self.s0t = np.zeros((md.mesh.numberofvertices))
            print("      no SMBd18opdd.s0t specified: values set as zero")
        return self
    # }}}
    def setdefaultparameters(self):  # {{{
        # pdd method not used in default mode
        self.ismungsm = 0
        self.isd18opd = 1
        self.istemperaturescaled = 1
        self.isprecipscaled = 1
        self.desfac = 0.5
        self.rlaps = 6.5
        self.rlapslgm = 6.5
        self.dpermil = 2.4
        self.f = 0.169
        self.issetpddfac = 0
        self.requested_outputs = ['default']
        return self
    # }}}
    def checkconsistency(self, md, solution, analyses):  # {{{
        if 'MasstransportAnalysis' in analyses:
            md = checkfield(md, 'fieldname', 'smb.desfac', '<=', 1, 'numel', [1])
            md = checkfield(md, 'fieldname', 'smb.s0p', '>=', 0, 'NaN', 1, 'Inf', 1, 'size', [md.mesh.numberofvertices])
            md = checkfield(md, 'fieldname', 'smb.s0t', '>=', 0, 'NaN', 1, 'Inf', 1, 'size', [md.mesh.numberofvertices])
            md = checkfield(md, 'fieldname', 'smb.rlaps', '>=', 0, 'numel', [1])
            md = checkfield(md, 'fieldname', 'smb.rlapslgm', '>=', 0, 'numel', [1])
            md = checkfield(md, 'fieldname', 'smb.steps_per_step', '>=', 1, 'numel', [1])
            md = checkfield(md, 'fieldname', 'smb.averaging', 'numel', [1], 'values', [0, 1, 2])

            if self.isd18opd:
                lent = float(np.size(self.temperatures_presentday, 1))
                lenp = float(np.size(self.precipitations_presentday, 1))
                multt = np.ceil(lent / 12.) * 12.
                multp = np.ceil(lenp / 12.) * 12.
                md = checkfield(md, 'fieldname', 'smb.temperatures_presentday', 'size', [md.mesh.numberofvertices, 12], 'NaN', 1, 'Inf', 1)
                md = checkfield(md, 'fieldname', 'smb.precipitations_presentday', 'size', [md.mesh.numberofvertices, 12], 'NaN', 1, 'Inf', 1)

                if self.istemperaturescaled == 0:
                    lent = float(np.size(self.temperatures_reconstructed, 1))
                    multt = np.ceil(lent / 12.) * 12.
                    md = checkfield(md, 'fieldname', 'smb.temperatures_reconstructed', 'size', [md.mesh.numberofvertices + 1, multt], 'NaN', 1, 'Inf', 1, 'timeseries', 1)

                if self.isprecipscaled == 0:
                    lenp = float(np.size(self.precipitations_reconstructed, 1))
                    multp = np.ceil(lent / 12.) * 12.
                    md = checkfield(md, 'fieldname', 'smb.precipitations_reconstructed', 'size', [md.mesh.numberofvertices + 1, multp], 'NaN', 1, 'Inf', 1, 'timeseries', 1)

                md = checkfield(md, 'fieldname', 'smb.delta18o', 'NaN', 1, 'Inf', 1, 'size', [2, np.nan], 'singletimeseries', 1)
                md = checkfield(md, 'fieldname', 'smb.dpermil', '>=', 0, 'numel', [1])
                md = checkfield(md, 'fieldname', 'smb.f', '>=', 0, 'numel', [1])

            if self.issetpddfac:
                md = checkfield(md, 'fieldname', 'smb.pddfac_snow', '>=', 0, 'NaN', 1, 'Inf', 1)
                md = checkfield(md, 'fieldname', 'smb.pddfac_ice', '>=', 0, 'NaN', 1, 'Inf', 1)

        md = checkfield(md, 'fieldname', 'masstransport.requested_outputs', 'stringrow', 1)
        return md
    # }}}
    def marshall(self, prefix, md, fid):  # {{{
        yts = md.constants.yts

        WriteData(fid, prefix, 'name', 'md.smb.model', 'data', 5, 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'ismungsm', 'format', 'Boolean')
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'isd18opd', 'format', 'Boolean')
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'issetpddfac', 'format', 'Boolean')
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'desfac', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 's0p', 'format', 'DoubleMat', 'mattype', 1)
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 's0t', 'format', 'DoubleMat', 'mattype', 1)
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'rlaps', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'rlapslgm', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'Tdiff', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', 2, 'yts', md.constants.yts)
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'sealev', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', 2, 'yts', md.constants.yts)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'steps_per_step', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'averaging', 'format', 'Integer')

        if self.isd18opd:
            WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'temperatures_presentday', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', md.constants.yts)
            WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'precipitations_presentday', 'format', 'DoubleMat', 'mattype', 1, 'scale', 1. / yts, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', md.constants.yts)
            WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'istemperaturescaled', 'format', 'Boolean')
            WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'isprecipscaled', 'format', 'Boolean')

            if self.istemperaturescaled == 0:
                WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'temperatures_reconstructed', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', md.constants.yts)

            if self.isprecipscaled == 0:
                WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'precipitations_reconstructed', 'format', 'DoubleMat', 'mattype', 1, 'scale', 1. / yts, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', md.constants.yts)

            WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'delta18o', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', 2, 'yts', md.constants.yts)
            WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'dpermil', 'format', 'Double')
            WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'f', 'format', 'Double')

        if self.issetpddfac:
            WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'pddfac_snow', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', md.constants.yts)
            WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'pddfac_ice', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', md.constants.yts)

        # Process requested outputs
        outputs = self.requested_outputs
        indices = [i for i, x in enumerate(outputs) if x == 'default']
        if len(indices) > 0:
            outputscopy = outputs[0:max(0, indices[0] - 1)] + self.defaultoutputs(md) + outputs[indices[0] + 1:]
            outputs = outputscopy
        WriteData(fid, prefix, 'data', outputs, 'name', 'md.smb.requested_outputs', 'format', 'StringArray')
    # }}}
