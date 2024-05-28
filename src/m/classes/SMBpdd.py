import numpy as np
from fielddisplay import fielddisplay
from checkfield import checkfield
from WriteData import WriteData
from project3d import project3d


class SMBpdd(object):
    """
    SMBpdd Class definition

       Usage:
          SMBpdd = SMBpdd()
    """

    def __init__(self):  # {{{
        self.precipitation = float('NaN')
        self.monthlytemperatures = float('NaN')
        self.desfac = 0.
        self.s0p = float('NaN')
        self.s0t = float('NaN')
        self.rlaps = 0.
        self.rlapslgm = 0.
        self.Pfac = float('NaN')
        self.Tdiff = float('NaN')
        self.sealev = float('NaN')
        self.isdelta18o = 0
        self.ismungsm = 0
        self.issetpddfac = 0
        self.delta18o = float('NaN')
        self.delta18o_surface = float('NaN')
        self.temperatures_presentday = float('NaN')
        self.temperatures_lgm = float('NaN')
        self.precipitations_presentday = float('NaN')
        self.precipitations_lgm = float('NaN')
        self.pddfac_snow = float('NaN')
        self.pddfac_ice = float('NaN')
        self.steps_per_step = 1
        self.averaging = 0
        self.requested_outputs = []

        # Set defaults
        self.setdefaultparameters()
    # }}}

    def __repr__(self):  # {{{
        string = "   surface forcings parameters:"

        string = "%s\n%s" % (string, fielddisplay(self, 'isdelta18o', 'is temperature and precipitation delta18o parametrisation activated (0 or 1, default is 0)'))
        string = "%s\n%s" % (string, fielddisplay(self, 'ismungsm', 'is temperature and precipitation mungsm parametrisation activated (0 or 1, default is 0)'))
        string = "%s\n%s" % (string, fielddisplay(self, 'issetpddfac', 'is user passing in defined pdd factors (0 or 1, default is 0)'))
        string = "%s\n%s" % (string, fielddisplay(self, 'desfac', 'desertification elevation factor (between 0 and 1, default is 0.5) [m]'))
        string = "%s\n%s" % (string, fielddisplay(self, 's0p', 'should be set to elevation from precip source (between 0 and a few 1000s m, default is 0) [m]'))
        string = "%s\n%s" % (string, fielddisplay(self, 's0t', 'should be set to elevation from temperature source (between 0 and a few 1000s m, default is 0) [m]'))
        string = "%s\n%s" % (string, fielddisplay(self, 'rlaps', 'present day lapse rate [degree/km]'))
        string = "%s\n%s" % (string, fielddisplay(self, 'rlapslgm', 'LGM lapse rate [degree/km]'))
        if not (self.isdelta18o and self.ismungsm):
            string = "%s\n%s" % (string, fielddisplay(self, 'monthlytemperatures', 'monthly surface temperatures [K], required if pdd is activated and delta18o not activated'))
            string = "%s\n%s" % (string, fielddisplay(self, 'precipitation', 'monthly surface precipitation [m/yr water eq], required if pdd is activated and delta18o or mungsm not activated'))
            if self.isdelta18o:
                string = "%s\n%s" % (string, fielddisplay(self, 'delta18o', 'delta18o [per mil], required if pdd is activated and delta18o activated'))
                string = "%s\n%s" % (string, fielddisplay(self, 'delta18o_surface', 'surface elevation of the delta18o site, required if pdd is activated and delta18o activated'))
                string = "%s\n%s" % (string, fielddisplay(self, 'temperatures_presentday', 'monthly present day surface temperatures [K], required if delta18o/mungsm is activated'))
                string = "%s\n%s" % (string, fielddisplay(self, 'temperatures_lgm', 'monthly LGM surface temperatures [K], required if delta18o or mungsm is activated'))
                string = "%s\n%s" % (string, fielddisplay(self, 'precipitations_presentday', 'monthly surface precipitation [m/yr water eq], required if delta18o or mungsm is activated'))
                string = "%s\n%s" % (string, fielddisplay(self, 'precipitations_lgm', 'monthly surface precipitation [m/yr water eq], required if delta18o or mungsm is activated'))
                string = "%s\n%s" % (string, fielddisplay(self, 'Tdiff', 'time interpolation parameter for temperature, 1D(year), required if mungsm is activated'))
                string = "%s\n%s" % (string, fielddisplay(self, 'sealev', 'sea level [m], 1D(year), required if mungsm is activated'))

            if self.ismungsm:
                string = "%s\n%s" % (string, fielddisplay(self, 'temperatures_presentday', 'monthly present day surface temperatures [K], required if delta18o/mungsm is activated'))
                string = "%s\n%s" % (string, fielddisplay(self, 'temperatures_lgm', 'monthly LGM surface temperatures [K], required if delta18o or mungsm is activated'))
                string = "%s\n%s" % (string, fielddisplay(self, 'precipitations_presentday', 'monthly surface precipitation [m/yr water eq], required if delta18o or mungsm is activated'))
                string = "%s\n%s" % (string, fielddisplay(self, 'precipitations_lgm', 'monthly surface precipitation [m/yr water eq], required if delta18o or mungsm is activated'))
                string = "%s\n%s" % (string, fielddisplay(self, 'Pfac', 'time interpolation parameter for precipitation, 1D(year), required if mungsm is activated'))
                string = "%s\n%s" % (string, fielddisplay(self, 'Tdiff', 'time interpolation parameter for temperature, 1D(year), required if mungsm is activated'))
                string = "%s\n%s" % (string, fielddisplay(self, 'sealev', 'sea level [m], 1D(year), required if mungsm is activated'))
                
        string = "%s\n%s" % (string, fielddisplay(self, 'steps_per_step', 'number of smb steps per time step'))
        string = "%s\n%s" % (string, fielddisplay(self, 'averaging', 'averaging methods from short to long steps'))
        string = "%s\n\t\t%s" % (string, '0: Arithmetic (default)')
        string = "%s\n\t\t%s" % (string, '1: Geometric')
        string = "%s\n\t\t%s" % (string, '2: Harmonic')

        string = "%s\n%s" % (string, fielddisplay(self, 'requested_outputs', 'additional outputs requested'))

        return string
    # }}}

    def extrude(self, md):  # {{{
        if not (self.isdelta18o and self.ismungsm):
            self.precipitation = project3d(md, 'vector', self.precipitation, 'type', 'node')
            self.monthlytemperatures = project3d(md, 'vector', self.monthlytemperatures, 'type', 'node')

        if self.isdelta18o:
            self.temperatures_lgm = project3d(md, 'vector', self.temperatures_lgm, 'type', 'node')
            self.temperatures_presentday = project3d(md, 'vector', self.temperatures_presentday, 'type', 'node')
            self.precipitations_presentday = project3d(md, 'vector', self.precipitations_presentday, 'type', 'node')
            self.precipitations_lgm = project3d(md, 'vector', self.precipitations_lgm, 'type', 'node')

        if self.ismungsm:
            self.temperatures_lgm = project3d(md, 'vector', self.temperatures_lgm, 'type', 'node')
            self.temperatures_presentday = project3d(md, 'vector', self.temperatures_presentday, 'type', 'node')
            self.precipitations_presentday = project3d(md, 'vector', self.precipitations_presentday, 'type', 'node')
            self.precipitations_lgm = project3d(md, 'vector', self.precipitations_lgm, 'type', 'node')

        if self.issetpddfac:
            self.pddfac_snow = project3d(md, 'vector', self.pddfac_snow, 'type', 'node')
        if self.issetpddfac:
            self.pddfac_ice = project3d(md, 'vector', self.pddfac_ice, 'type', 'node')
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
            print("      no SMBpdd.s0p specified: values set as zero")

        if np.all(np.isnan(self.s0t)):
            self.s0t = np.zeros((md.mesh.numberofvertices))
            print("      no SMBpdd.s0t specified: values set as zero")

        return self
    # }}}

    def setdefaultparameters(self):  # {{{
        #pdd method not used in default mode
        self.isdelta18o = 0
        self.ismungsm = 0
        self.desfac = 0.5
        self.rlaps = 6.5
        self.rlapslgm = 6.5
        self.issetpddfac = 0
        self.requested_outputs = ['default']

        return self
    # }}}

    def checkconsistency(self, md, solution, analyses):    # {{{

        if 'MasstransportAnalysis' in analyses:
            md = checkfield(md, 'fieldname', 'smb.desfac', '<=', 1, 'numel', [1])
            md = checkfield(md, 'fieldname', 'smb.s0p', '>=', 0, 'NaN', 1, 'Inf', 1, 'size', [md.mesh.numberofvertices])
            md = checkfield(md, 'fieldname', 'smb.s0t', '>=', 0, 'NaN', 1, 'Inf', 1, 'size', [md.mesh.numberofvertices])
            md = checkfield(md, 'fieldname', 'smb.rlaps', '>=', 0, 'numel', [1])
            md = checkfield(md, 'fieldname', 'smb.rlapslgm', '>=', 0, 'numel', [1])

            if (self.isdelta18o == 0 and self.ismungsm == 0):
                md = checkfield(md, 'fieldname', 'smb.monthlytemperatures', 'size', [md.mesh.numberofvertices, 12], 'NaN', 1, 'Inf', 1)
                md = checkfield(md, 'fieldname', 'smb.precipitation', 'size', [md.mesh.numberofvertices, 12], 'NaN', 1, 'Inf', 1)
            elif self.isdelta18o:
                md = checkfield(md, 'fieldname', 'smb.delta18o', 'NaN', 1, 'Inf', 1, 'size', [2, np.nan], 'singletimeseries', 1)
                md = checkfield(md, 'fieldname', 'smb.delta18o_surface', 'NaN', 1, 'Inf', 1, 'size', [2, np.nan], 'singletimeseries', 1)
                md = checkfield(md, 'fieldname', 'smb.temperatures_presentday', 'size', [md.mesh.numberofvertices, 12], 'NaN', 1, 'Inf', 1)
                md = checkfield(md, 'fieldname', 'smb.temperatures_lgm', 'size', [md.mesh.numberofvertices, 12], 'NaN', 1, 'Inf', 1)
                md = checkfield(md, 'fieldname', 'smb.precipitations_presentday', 'size', [md.mesh.numberofvertices, 12], 'NaN', 1, 'Inf', 1)
                md = checkfield(md, 'fieldname', 'smb.precipitations_lgm', 'size', [md.mesh.numberofvertices, 12], 'NaN', 1, 'Inf', 1)
                md = checkfield(md, 'fieldname', 'smb.Tdiff', 'NaN', 1, 'Inf', 1, 'size', [2, np.nan], 'singletimeseries', 1)
                md = checkfield(md, 'fieldname', 'smb.sealev', 'NaN', 1, 'Inf', 1, 'size', [2, np.nan], 'singletimeseries', 1)
            elif self.ismungsm:
                md = checkfield(md, 'fieldname', 'smb.temperatures_presentday', 'size', [md.mesh.numberofvertices, 12], 'NaN', 1, 'Inf', 1)
                md = checkfield(md, 'fieldname', 'smb.temperatures_lgm', 'size', [md.mesh.numberofvertices, 12], 'NaN', 1, 'Inf', 1)
                md = checkfield(md, 'fieldname', 'smb.precipitations_presentday', 'size', [md.mesh.numberofvertices, 12], 'NaN', 1, 'Inf', 1)
                md = checkfield(md, 'fieldname', 'smb.precipitations_lgm', 'size', [md.mesh.numberofvertices, 12], 'NaN', 1, 'Inf', 1)
                md = checkfield(md, 'fieldname', 'smb.Pfac', 'NaN', 1, 'Inf', 1, 'size', [2, np.nan], 'singletimeseries', 1)
                md = checkfield(md, 'fieldname', 'smb.Tdiff', 'NaN', 1, 'Inf', 1, 'size', [2, np.nan], 'singletimeseries', 1)
                md = checkfield(md, 'fieldname', 'smb.sealev', 'NaN', 1, 'Inf', 1, 'size', [2, np.nan], 'singletimeseries', 1)

            if self.issetpddfac:
                md = checkfield(md, 'fieldname', 'smb.pddfac_snow', '>=', 0, 'NaN', 1, 'Inf', 1)
                md = checkfield(md, 'fieldname', 'smb.pddfac_ice', '>=', 0, 'NaN', 1, 'Inf', 1)

        md = checkfield(md, 'fieldname', 'smb.steps_per_step', '>=', 1, 'numel', [1])
        md = checkfield(md, 'fieldname', 'smb.averaging', 'numel', [1], 'values', [0, 1, 2])
        md = checkfield(md, 'fieldname', 'masstransport.requested_outputs', 'stringrow', 1)
        return md
    # }}}

    def marshall(self, prefix, md, fid):    # {{{
        yts = md.constants.yts

        WriteData(fid, prefix, 'name', 'md.smb.model', 'data', 4, 'format', 'Integer')

        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'isdelta18o', 'format', 'Boolean')
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'ismungsm', 'format', 'Boolean')
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'issetpddfac', 'format', 'Boolean')
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'desfac', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 's0p', 'format', 'DoubleMat', 'mattype', 1)
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 's0t', 'format', 'DoubleMat', 'mattype', 1)
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'rlaps', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'rlapslgm', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'steps_per_step', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'averaging', 'format', 'Integer')

        if (self.isdelta18o == 0 and self.ismungsm == 0):
            WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'monthlytemperatures', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', yts)
            WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'precipitation', 'format', 'DoubleMat', 'mattype', 1, 'scale', 1. / yts, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', yts)
        elif self.isdelta18o:
            WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'temperatures_presentday', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', yts)
            WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'temperatures_lgm', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', yts)
            WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'precipitations_presentday', 'format', 'DoubleMat', 'mattype', 1, 'scale', 1. / yts, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', yts)
            WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'precipitations_lgm', 'format', 'DoubleMat', 'mattype', 1, 'scale', 1. / yts, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', yts)
            WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'delta18o_surface', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', 2, 'yts', md.constants.yts)
            WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'delta18o', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', 2, 'yts', md.constants.yts)
            WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'Tdiff', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', 2, 'yts', md.constants.yts)
            WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'sealev', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', 2, 'yts', md.constants.yts)
        elif self.ismungsm:
            WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'temperatures_presentday', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', yts)
            WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'temperatures_lgm', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', yts)
            WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'precipitations_presentday', 'format', 'DoubleMat', 'mattype', 1, 'scale', 1. / yts, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', yts)
            WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'precipitations_lgm', 'format', 'DoubleMat', 'mattype', 1, 'scale', 1. / yts, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', yts)
            WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'Pfac', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', 2, 'yts', md.constants.yts)
            WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'Tdiff', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', 2, 'yts', md.constants.yts)
            WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'sealev', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', 2, 'yts', md.constants.yts)

        if self.issetpddfac:
            WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'pddfac_snow', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', md.constants.yts)
            WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'pddfac_ice', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', md.constants.yts)

        #process requested outputs
        outputs = self.requested_outputs
        indices = [i for i, x in enumerate(outputs) if x == 'default']
        if len(indices) > 0:
            outputscopy = outputs[0:max(0, indices[0] - 1)] + self.defaultoutputs(md) + outputs[indices[0] + 1:]
            outputs = outputscopy
        WriteData(fid, prefix, 'data', outputs, 'name', 'md.smb.requested_outputs', 'format', 'StringArray')

    # }}}
