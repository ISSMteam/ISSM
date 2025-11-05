from math import log

import numpy as np

from checkfield import checkfield
from fielddisplay import fielddisplay
from project3d import project3d
from WriteData import WriteData


class SMBsemic(object):
    """SMBsemic class definition

    Usage:
        SMBsemic = SMBsemic()
    """

    def __init__(self, *args):  # {{{
        self.dailysnowfall = np.nan
        self.dailyrainfall = np.nan
        self.dailydsradiation = np.nan
        self.dailydlradiation = np.nan
        self.dailypressure = np.nan
        self.dailyairdensity = np.nan
        self.dailyairhumidity = np.nan
        self.dailytemperature = np.nan

        self.Tamp = np.nan
        self.mask = np.nan
        self.hice = np.nan
        self.hsnow = np.nan
        self.desfac = 0
        self.desfacElevation = 0
        self.rlaps = 0
        self.rdl = 0
        self.s0gcm = np.nan
        self.steps_per_step = 1
        self.averaging = 0
        self.requested_outputs = []

        self.hcrit = 0
        self.rcrit = 0

        # albedo
        self.albedo = 0 # required for first energy balance calculation of SEMIC
        self.albedo_snow = 0 # required for ISBA method
        self.albedo_scheme = 0
        self.alb_smax = np.nan
        self.alb_smin = np.nan
        self.albi = np.nan
        self.albl = np.nan

        # albedo parameters depending on albedo_scheme
        # for slater
        self.tmin = np.nan
        self.tmax = np.nan

        # for isba & denby method
        self.mcrit = np.nan

        # for isba
        self.tau_a = np.nan
        self.tau_f = np.nan
        self.wcrit = np.nan

        # for alex
        self.tmid = np.nan
        self.afac = np.nan

        # method
        self.ismethod = 0
        self.isdesertification = 0
        self.isLWDcorrect = 0

        if len(args) == 0:
            self.setdefaultparameters()
        else:
            raise Exception('constructor not supported')
    # }}}

    def __repr__(self):  # {{{
        s = '   surface forcings parameters:\n'
        s += '   Interface for coupling GCM data to the energy balance model SEMIC (Krapp et al (2017) https://doi.org/10.5194/tc-11-1519-2017).\n'
        s += '   The implemented coupling uses daily mean GCM input to calculate yearly mean smb, accumulation, ablation, and surface temperature.\n'
        s += '   smb and temperatures are updated every year\n'
        s += '\n   SEMIC parameters:\n'
        s += '{}\n'.format(fielddisplay(self, 'dailysnowfall', 'daily surface dailysnowfall [m/s]'))
        s += '{}\n'.format(fielddisplay(self, 'dailyrainfall', 'daily surface dailyrainfall [m/s]'))
        s += '{}\n'.format(fielddisplay(self, 'dailydsradiation', 'daily downwelling shortwave radiation [W/m2]'))
        s += '{}\n'.format(fielddisplay(self, 'dailydlradiation', 'daily downwelling longwave radiation [W/m2]'))
        s += '{}\n'.format(fielddisplay(self, 'dailywindspeed', 'daily surface wind speed [m/s]'))
        s += '{}\n'.format(fielddisplay(self, 'dailypressure', 'daily surface pressure [Pa]'))
        s += '{}\n'.format(fielddisplay(self, 'dailyairdensity', 'daily air density [kg/m3]'))
        s += '{}\n'.format(fielddisplay(self, 'dailyairhumidity', 'daily air specific humidity [kg/kg]'))
        s += '{}\n'.format(fielddisplay(self, 'rlaps', 'present day lapse rate (default is 7.4 [degree/km]; Erokhina et al. 2017)'))
        s += '{}\n'.format(fielddisplay(self, 'desfac', 'desertification elevation factor (default is -log(2.0)/1000 [1/m]; Vizcaino et al. 2010)'))
        s += '{}\n'.format(fielddisplay(self, 'rdl', 'longwave downward radiation decrease (default is 29 [W/m^2/km]; Marty et al. 2002)'))
        s += '{}\n'.format(fielddisplay(self, 's0gcm', 'GCM reference elevation; (default is 0) [m]'))
        s += '{}\n'.format(fielddisplay(self, 'ismethod','method for calculating SMB with SEMIC. Default version of SEMIC is really slow. 0: steady, 1: transient (default: 0)'))
        if self.ismethod: # transient mode
            s += '{}\n'.format(fielddisplay(self,'desfacElevation','desertification elevation (default is 2000 m; Vizcaino et al. 2010)'))
            s += '{}\n'.format(fielddisplay(self,'Tamp','amplitude of diurnal cycle [K]'))
            s += '{}\n'.format(fielddisplay(self,'albedo','initial albedo [no unit]'))
            s += '{}\n'.format(fielddisplay(self,'albedo_snow','initial albedo for snow [no unit]'))
            s += '{}\n'.format(fielddisplay(self,'hice','initial thickness of ice [unit: m]'))
            s += '{}\n'.format(fielddisplay(self,'hsnow','initial thickness of snow [unit: m]'))
            s += '{}\n'.format(fielddisplay(self,'mask','masking for albedo. 0: ocean, 1: land, 2: ice (default: 2)'))
            s += '{}\n'.format(fielddisplay(self,'hcrit','critical snow height for albedo [unit: m]'))
            s += '{}\n'.format(fielddisplay(self,'rcrit','critical refreezing height for albedo [no unit]'))

            s += '\nSEMIC albedo parameters.\n'
            s += '{}\n'.format(fielddisplay(self,'albedo_scheme','albedo scheme for SEMIC. 0: none, 1: slater, 2: denby, 3: isba, 4: alex (default is 0)'))
            s += '{}\n'.format(fielddisplay(self,'alb_smax','maximum snow albedo (default: 0.79)'))
            s += '{}\n'.format(fielddisplay(self,'alb_smin','minimum snow albedo (default: 0.6)'))
            s += '{}\n'.format(fielddisplay(self,'albi','background albedo for bare ice (default: 0.41)'))
            s += '{}\n'.format(fielddisplay(self,'albl','background albedo for bare land (default: 0.07)'))
            
            s += '{}\n'.format(fielddisplay(self,'isdesertification','enable or disable desertification of Vizcaino et al. (2010). 0: off, 1: on (default: 1)'))
            s += '{}\n'.format(fielddisplay(self,'isLWDcorrect','enable or disable downward longwave correction of Marty et al. (2002). 0: off, 1: on (default: 1)'))
        # albedo_scheme - 0: none, 1: slater, 2: isba, 3: denby, 4: alex.
        if self.albedo_scheme == 0:
            s += '\n\tSEMIC snow albedo parameter of None.\n'
            s += '\t   albedo of snow is updated from albedo snow max (alb_smax).\n'
            s += '\t   alb_snow = abl_smax \n '
        elif self.albedo_scheme == 1:
            s += '\n\tSEMIC snow albedo parameters of Slater et al, (1998).\n'
            s += '\t   alb = alb_smax - (alb_smax - alb_smin)*tm^(3.0)\n'
            s += '\t   tm  = 1 (tsurf > 273.15 K)\n'
            s += '\t         tm = f*(tsurf-tmin) (tmin <= tsurf < 273.15)\n'
            s += '\t         0 (tsurf < tmin)\n'
            s += '\t   f = 1/(273.15-tmin)\n'
            s += '{}\n'.format(fielddisplay(self, 'tmin', 'minimum temperature for which albedo decline become effective. (default: 263.15 K)[unit: K])'))
            s += '{}\n'.format(fielddisplay(self, 'tmax', 'maxmium temperature for which albedo decline become effective. This value should be fixed. (default: 273.15 K)[unit: K])'))
        elif self.albedo_scheme == 2:
            s += '\n\tSEMIC snow albedo parameters of Denby et al. (2002 Tellus)\n'
            s += '{}\n'.format(fielddisplay(self,'mcrit','critical melt rate (defaut: 6e-8) [unit: m/sec]'))
        elif self.albedo_scheme == 3:
            s += '\n\tSEMIC snow albedo parameters of ISB (Douville et al., 1995).\n'
            s += '{}\n'.format(fielddisplay(self, 'mcrit', 'critical melt rate (default: 6e-8) [unit: m/sec]'))
            s += '{}\n'.format(fielddisplay(self, 'wcrit', 'critical liquid water content (default: 15) [unit: kg/m2]'))
            s += '{}\n'.format(fielddisplay(self, 'tau_a', 'dry albedo decline [unit: 1/day]'))
            s += '{}\n'.format(fielddisplay(self, 'tau_f', 'wet albedo decline [unit: 1/day]'))
            s += '\n\tReference'
            s += '\tDouville, H., Royer, J.-F., and Mahfouf, J.-F.: A new snow parameterization for the Météo-France climate model. Part I: validation in stand-alone experiments, Climate Dynamics, 12, 21–35, https://doi.org/10.1007/s003820050092, 1995.'
        elif self.albedo_scheme == 4:
            s += '\n\tSEMIC snow albedo parameters of Alex.?\n'
            s += '{}\n'.format(fielddisplay(self,'afac','[unit: ?]'))
            s += '{}\n'.format(fielddisplay(self,'tmid','[unit: ?]'))
        else:
            raise Exception('ERROR: {} is not supported albedo scheme.'.format(self.albedo_scheme))

        s += '{}\n'.format(fielddisplay(self, 'steps_per_step', 'number of smb steps per time step'))
        s += '{}\n'.format(fielddisplay(self, 'averaging', 'averaging methods from short to long steps'))
        s += '\t\t{}\n'.format('0: Arithmetic (default)')
        s += '\t\t{}\n'.format('1: Geometric')
        s += '\t\t{}\n'.format('2: Harmonic')
        s += '{}\n'.format(fielddisplay(self, 'requested_outputs', 'additional outputs requested'))
        return s
    # }}}

    def extrude(self, md):  # {{{
        self.dailysnowfall = project3d(md, 'vector', self.dailysnowfall, 'type', 'node')
        self.dailyrainfall = project3d(md, 'vector', self.dailyrainfall, 'type', 'node')
        self.dailydsradiation = project3d(md, 'vector', self.dailydsradiation, 'type', 'node')
        self.dailydlradiation = project3d(md, 'vector', self.dailydlradiation, 'type', 'node')
        self.dailywindspeed = project3d(md, 'vector', self.dailywindspeed, 'type', 'node')
        self.dailypressure = project3d(md, 'vector', self.dailypressure, 'type', 'node')
        self.dailyairdensity = project3d(md, 'vector', self.dailyairdensity, 'type', 'node')
        self.dailyairhumidity = project3d(md, 'vector', self.dailyairhumidity, 'type', 'node')
        self.dailytemperature = project3d(md, 'vector', self.dailytemperature, 'type', 'node')
        self.s0gcm = project3d(md, 'vector', self.s0gcm, 'type', 'node')
        return self
    # }}}

    def defaultoutputs(self, md):  # {{{
        return ['SmbMassBalance']
    # }}}

    def outputlists(self, md):  # {{{
        if self.ismethod:
            list = ['default','SmbMassBalance', 'SmbMassBalanceSnow', 'SmbMassBalanceIce',
                  'SmbMelt', 'SmbRefreeze','SmbAccumulation',
                  'SmbHIce', 'SmbHSnow', 'SmbAlbedo', 'SmbAlbedoSnow', 'TemperatureSEMIC',
                  'SmbSemicQmr', 'TotalSmb', 'TotalSmbMelt', 'TotalSmbRefreeze', 
                  'SmbRunoff','SmbEvaporation']
        else:
            list = ['default','SmbMassBalance']
        return list
    # }}}

    def initialize(self, md):  # {{{
        if np.all(np.isnan(self.s0gcm)):
            self.s0gcm = np.zeros((md.mesh.numberofvertices))
            print('      no SMBsemic.s0gcm specified: values set as zero')

        self.Tamp = 3 * np.ones((md.mesh.numberofvertices,))
        #self.albedo = 0.8 * np.ones((md.mesh.numberofvertices,))
        #self.albedo_snow = 0.5 * np.ones((md.mesh.numberofvertices,))
        self.hice = np.zeros((md.mesh.numberofvertices,))
        self.hsnow = 5 * np.ones((md.mesh.numberofvertices,))

        return self
    # }}}

    def setdefaultparameters(self):  #{{{
        # albedo parameters
        self.albedo_scheme = 0
        self.alb_smax = 0.79
        self.alb_smin = 0.6
        self.albi = 0.41
        self.albl = 0.07

        # albedo parameters for?
        # for slater
        self.tmin = 263.15
        self.tmax = 273.15

        # for isba & denby
        self.mcrit = 6e-8

        # for isba
        self.tau_a = 0.008
        self.tau_f = 0.24
        self.wcrit = 15.0

        # for alex
        self.tmid = 273.35
        self.afac = -0.18

        self.hcrit = 0.028 # from Krapp et al. (2017)
        self.rcrit = 0.85 # from Krapp et al. (2017)

        self.desfac = -log(2.0) / 1000
        self.desfacElevation = 2000
        self.rlaps = 7.4
        self.rdl   = 29 # from Marty et al. (2002)

        self.ismethod = 0
        self.isdesertification = 1
        self.isLWDcorrect      = 1
        self.requested_outputs = ['default']
        return self
    # }}}

    def checkconsistency(self, md, solution, analyses):  # {{{
        if 'MasstransportAnalysis' in analyses:
            md = checkfield(md, 'fieldname', 'smb.desfac', '<=', 1, 'numel', 1)
            md = checkfield(md, 'fieldname', 'smb.s0gcm', '>=', 0, 'NaN', 1, 'Inf', 1, 'size', [md.mesh.numberofvertices, 1])
            md = checkfield(md, 'fieldname', 'smb.rlaps', '>=', 0, 'numel', 1);
            md = checkfield(md, 'fieldname', 'smb.rdl', '>=', 0, 'numel', 1);
            md = checkfield(md, 'fieldname', 'smb.dailysnowfall','timeseries',1,'NaN',1,'Inf',1,'>=',0);
            md = checkfield(md, 'fieldname', 'smb.dailyrainfall','timeseries',1,'NaN',1,'Inf',1,'>=',0);
            md = checkfield(md, 'fieldname', 'smb.dailydsradiation','timeseries',1,'NaN',1,'Inf',1,'>=',0);
            md = checkfield(md, 'fieldname', 'smb.dailydlradiation','timeseries',1,'NaN',1,'Inf',1,'>=',0);
            md = checkfield(md, 'fieldname', 'smb.dailywindspeed','timeseries',1,'NaN',1,'Inf',1,'>=',0);
            md = checkfield(md, 'fieldname', 'smb.dailypressure','timeseries',1,'NaN',1,'Inf',1,'>=',0);
            md = checkfield(md, 'fieldname', 'smb.dailyairdensity','timeseries',1,'NaN',1,'Inf',1,'>=',0);
            md = checkfield(md, 'fieldname', 'smb.dailyairhumidity','timeseries',1,'NaN',1,'Inf',1,'>=',0);
            md = checkfield(md, 'fieldname', 'smb.dailytemperature','timeseries',1,'NaN',1,'Inf',1,'>=',0);

            # TODO: transient model should be merged with SEMIC developed by Ruckamp et al. (2018)
            md = checkfield(md, 'fieldname', 'smb.ismethod', 'numel', 1, 'values', [0, 1])
            md = checkfield(md, 'fieldname', 'smb.isdesertification', 'Nan',1, 'Inf', 1, 'numel', 1, 'values', [0, 1])
            md = checkfield(md, 'fieldname', 'smb.isLWDcorrect', 'Nan',1, 'Inf',1, 'numel', 1, 'values', [0, 1])
            if self.ismethod: # transient mode
                md = checkfield(md, 'fieldname', 'smb.desfacElevation', '>=', 0, 'numel', 1)
                md = checkfield(md, 'fieldname', 'smb.albedo_scheme', 'NaN', 1, 'Inf', 1, 'numel', 1, 'values', [0, 1, 2, 3, 4])
                md = checkfield(md, 'fieldname', 'smb.alb_smax', '>=', 0, 'NaN', 1, 'Inf', 1, 'numel', 1)
                md = checkfield(md, 'fieldname', 'smb.mask', 'NaN', 1, 'Inf', 1, 'size', [md.mesh.numberofvertices, 1], 'values', [0, 1, 2])

                # initial values
                md = checkfield(md, 'fieldname', 'smb.albedo', 'NaN', 1, 'Inf', 1, 'size', [md.mesh.numberofvertices, 1])
                md = checkfield(md, 'fieldname', 'smb.albedo_snow', 'NaN', 1, 'Inf', 1, 'size', [md.mesh.numberofvertices, 1])
                md = checkfield(md, 'fieldname', 'smb.alb_smax', '>=', 0, '<=', 1, 'NaN', 1, 'Inf', 1, 'numel', 1)
                md = checkfield(md, 'fieldname', 'smb.alb_smin', '>=', 0, '<=', 1, 'NaN', 1, 'Inf', 1, 'numel', 1)
                md = checkfield(md, 'fieldname', 'smb.albi', '>=', 0, '<=', 1, 'NaN', 1, 'Inf', 1, 'numel', 1)
                md = checkfield(md, 'fieldname', 'smb.albl', '>=', 0, '<=', 1, 'NaN', 1, 'Inf', 1, 'numel', 1)
                md = checkfield(md, 'fieldname', 'smb.hice', 'NaN', 1, 'Inf', 1, 'size', [md.mesh.numberofvertices, 1])
                md = checkfield(md, 'fieldname', 'smb.hsnow', 'NaN', 1, 'Inf', 1, 'size', [md.mesh.numberofvertices, 1])
        md = checkfield(md, 'fieldname', 'smb.steps_per_step', '>=', 1, 'numel', [1])
        md = checkfield(md, 'fieldname', 'smb.averaging', 'numel', [1], 'values', [0, 1, 2])
        md = checkfield(md, 'fieldname', 'smb.requested_outputs', 'stringrow', 1)
        return md
    # }}}

    def marshall(self, prefix, md, fid):  # {{{
        yts = md.constants.yts
        WriteData(fid, prefix, 'name', 'md.smb.model', 'data', 12, 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'ismethod', 'format', 'Integer', 'values', [0, 1])
        WriteData(fid, prefix, 'object', self, 'class' ,'smb', 'fieldname', 'desfac', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'desfacElevation', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'class' ,'smb', 'fieldname', 's0gcm', 'format', 'DoubleMat', 'mattype', 1)
        WriteData(fid, prefix, 'object', self, 'class' ,'smb', 'fieldname', 'rlaps', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'class' ,'smb', 'fieldname', 'rdl', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'class' ,'smb', 'fieldname', 'dailysnowfall', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', yts)
        WriteData(fid, prefix, 'object', self, 'class' ,'smb', 'fieldname', 'dailyrainfall', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', yts)
        WriteData(fid, prefix, 'object', self, 'class' ,'smb', 'fieldname', 'dailydsradiation', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', yts)
        WriteData(fid, prefix, 'object', self, 'class' ,'smb', 'fieldname', 'dailydlradiation', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', yts)
        WriteData(fid, prefix, 'object', self, 'class' ,'smb', 'fieldname', 'dailywindspeed', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', yts)
        WriteData(fid, prefix, 'object', self, 'class' ,'smb', 'fieldname', 'dailypressure', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', yts)
        WriteData(fid, prefix, 'object', self, 'class' ,'smb', 'fieldname', 'dailyairdensity', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', yts)
        WriteData(fid, prefix, 'object', self, 'class' ,'smb', 'fieldname', 'dailyairhumidity', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', yts)
        WriteData(fid, prefix, 'object', self, 'class' ,'smb', 'fieldname', 'dailytemperature', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', yts)
        # TODO: transient mode should be merged with SEMIC developed by Ruckamp et al. (2018).
        if self.ismethod: # transient mode
            WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'Tamp', 'format', 'DoubleMat', 'mattype', 1)
            WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'mask', 'format', 'DoubleMat', 'mattype', 1)
            WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'hice', 'format', 'DoubleMat', 'mattype', 1)
            WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'hsnow', 'format', 'DoubleMat', 'mattype', 1)
            WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'qmr', 'format', 'DoubleMat', 'mattype', 1)

            WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'hcrit', 'format', 'Double')
            WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'rcrit', 'format', 'Double')

            # albedo
            WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'albedo', 'format', 'DoubleMat', 'mattype', 1)
            WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'albedo_snow', 'format', 'DoubleMat', 'mattype', 1)
            WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'albedo_scheme', 'format', 'Integer')
            WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'albi', 'format', 'Double')
            WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'albl', 'format', 'Double')
            WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'alb_smin', 'format', 'Double')
            WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'alb_smax', 'format', 'Double')

            # albedo parameters for ?
            # for slater
            WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'tmin', 'format', 'Double')
            WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'tmax', 'format', 'Double')
            # for isba & denby
            WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'mcrit', 'format', 'Double')
            # for isba
            WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'wcrit', 'format', 'Double')
            WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'tau_a', 'format', 'Double')
            WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'tau_f', 'format', 'Double')
            # for alex
            WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'tmid', 'format', 'Double')
            WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'afac', 'format', 'Double')
        #specific parameterization
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname','isdesertification', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname','isLWDcorrect', 'format', 'Integer')

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
