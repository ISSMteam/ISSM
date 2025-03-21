import numpy as np

from checkfield import checkfield
from fielddisplay import fielddisplay
from project3d import project3d
from WriteData import WriteData


class SMBgemb(object):
    """SMBGEMB class definition

    Usage:
        SMB = SMBgemb(md.mesh)
    """

    def __init__(self, *args):  # {{{
        """Each one of these properties is a transient forcing to the GEMB 
        model, loaded from meteorological data derived from an automatic 
        weather stations (AWS). Each property is therefore a matrix, of size 
        (numberofvertices x number of time steps.
        """

        #solution choices
        self.isgraingrowth       = 0
        self.isalbedo            = 0
        self.isshortwave         = 0
        self.isthermal           = 0
        self.isaccumulation      = 0
        self.ismelt              = 0
        self.isdensification     = 0
        self.isturbulentflux     = 0
        self.isconstrainsurfaceT = 0
        self.isdeltaLWup         = 0
        self.ismappedforcing     = 0
        self.iscompressedforcing = 0

        # Inputs
        self.Ta                     = np.nan    # 2 m air temperature, in Kelvin
        self.V                      = np.nan    # wind speed (m/s-1)
        self.dswrf                  = np.nan    # downward shortwave radiation flux [W/m^2]
        self.dlwrf                  = np.nan    # downward longwave radiation flux [W/m^2]
        self.P                      = np.nan    # precipitation [mm w.e. / m^2]
        self.eAir                   = np.nan    # screen level vapor pressure [Pa]
        self.pAir                   = np.nan    # surface pressure [Pa]

        self.Tmean                  = np.nan    # mean annual temperature [K]
        self.Vmean                  = np.nan    # mean annual wind velocity [m s-1]
        self.C                      = np.nan    # mean annual snow accumulation [kg m-2 yr-1]
        self.Tz                     = np.nan    # height above ground at which temperature (T) was sampled [m]
        self.Vz                     = np.nan    # height above ground at which wind (V) was sampled [m]

        # Optional inputs
        self.aValue                 = np.nan    # Albedo forcing at every element. Used only if aIdx == 0, or density exceeds adThresh.
        self.teValue                = np.nan    # Outward longwave radiation thermal emissivity forcing at every element (default in code is 1), Used only if eIdx== 0, or effective grain radius exceeds teThresh
        self.dulwrfValue            = np.nan    #Delta with which to perturb the long wave radiation upwards. Use if isdeltaLWup is true.
        self.mappedforcingpoint     = np.nan    #Mapping of which forcing point will map to each mesh element (integer). Of size number of elements. Use if ismappedforcing is true.
        self.mappedforcingelevation = np.nan    #The elevation of each mapped forcing location (m above sea level). Of size number of forcing points. Use if ismappedforcing is true.
        self.lapseTaValue           = np.nan    #Temperature lapse rate if forcing has different grid and should be remapped. Use if ismappedforcing is true. (Default value is -0.006 K m-1.)
        self.lapsedlwrfValue        = np.nan    #Longwave down lapse rate if forcing has different grid and should be remapped. Use if ismappedforcing is true. (Default value is -0.032 W m-2 m-1.)

        # Initialization of snow properties
        self.Dzini                  = np.nan    # cell depth (m)
        self.Dini                   = np.nan    # snow density (kg m-3)
        self.Reini                  = np.nan    # effective grain size (mm)
        self.Gdnini                 = np.nan    # grain dendricity (0-1)
        self.Gspini                 = np.nan    # grain sphericity (0-1)
        self.ECini                  = np.nan    # evaporation/condensation (kg m-2)
        self.Wini                   = np.nan    # Water content (kg m-2)
        self.Aini                   = np.nan    # albedo (0-1)
        self.Adiffini               = np.nan    # albedo, diffusive radiation (0-1)
        self.Tini                   = np.nan    # snow temperature (K)
        self.Sizeini                = np.nan    # Number of layers

        # Settings
        self.aIdx                   = np.nan    # method for calculating albedo and subsurface absorption (default is 1)
        # 0: direct input from aValue parameter, no use of adThresh
        # 1: effective grain radius [Gardner & Sharp, 2009]
        # 2: effective grain radius [Brun et al., 1992; LeFebre et al., 2003], with swIdx=1, SW penetration follows grain size in 3 spectral bands (Brun et al., 1992)
        # 3: density and cloud amount [Greuell & Konzelmann, 1994]
        # 4: exponential time decay & wetness [Bougamont & Bamber, 2005]
        
        self.eIdx                   = np.nan    #method for calculating emissivity (default is 1)
        # 0: direct input from teValue parameter, no use of teThresh
        # 1: default value of 1, in areas with grain radius below teThresh
        # 2: default value of 1, in areas with grain radius below teThresh and areas of dry snow (not bare ice or wet) at the surface

        self.tcIdx                   = np.nan    #method for calculating thermal conductivity (default is 1)
        # 1: after Sturm et al, 1997
        # 2: after Calonne et al., 2011 

        self.swIdx                  = np.nan    # apply all SW to top grid cell (0) or allow SW to penetrate surface (1) (default 0, if swIdx=1 and aIdx=2, function of effective radius (Brun et al., 1992) or else dependent on snow density (taken from Bassford, 2002)) 

        self.denIdx                 = np.nan    # densification model to use (default is 2):
        # 1 = emperical model of Herron and Langway (1980)
        # 2 = semi-emperical model of Anthern et al. (2010)
        # 3 = DO NOT USE: physical model from Appendix B of Anthern et al. (2010)
        # 4 = DO NOT USE: emperical model of Li and Zwally (2004)
        # 5 = DO NOT USE: modified emperical model (4) by Helsen et al. (2008)
        # 6 = Antarctica semi-emperical model of Ligtenberg et al. (2011)
        # 7 = Greenland semi-emperical model of Kuipers Munneke et al. (2015)

        self.dsnowIdx               = np.nan    # model for fresh snow accumulation density (default is 1):
        # 0 = Original GEMB value, 150 kg/m^3
        # 1 = Antarctica value of fresh snow density, 350 kg/m^3
        # 2 = Greenland value of fresh snow density, 315 kg/m^3, Fausto et al. (2018)
        # 3 = Antarctica model of Kaspers et al. (2004)
        # 4 = Greenland model of Kuipers Munneke et al. (2015)

        self.zTop                   = np.nan    # depth over which grid length is constant at the top of the snopack (default 10) [m]
        self.dzTop                  = np.nan    # initial top vertical grid spacing (default .05) [m]
        self.dzMin                  = np.nan    # initial min vertical allowable grid spacing (default dzMin/2) [m]

        self.zY                     = np.nan    # stretch grid cells bellow top_z by a [top_dz * y ^ (cells bellow top_z)]
        self.zMax                   = np.nan    # initial max model depth (default is min(thickness, 250)) [m]
        self.zMin                   = np.nan    # initial min model depth (default is min(thickness, 130)) [m]
        self.outputFreq             = np.nan    # output frequency in days (default is monthly, 30)

        # Specific albedo parameters
        # Method 1
        self.dswdiffrf              = np.nan    # downward diffusive shortwave radiation flux [W/m^2]
        self.szaValue               = np.nan    # Solar Zenith Angle [degree]
        self.cotValue               = np.nan    # Cloud Optical Thickness
        self.ccsnowValue            = np.nan    # concentration of light absorbing carbon for snow [ppm1]
        self.cciceValue             = np.nan    # concentration of light absorbing carbon for ice [ppm1]
        # Method 1 and 2
        self.aSnow                  = np.nan    # new snow albedo (0.64 - 0.89)
        self.aIce                   = np.nan    # range 0.27-0.58 for old snow
        #Method 3: Radiation Correction Factors -> only used for met station data and Greuell & Konzelmann, 1994 albedo
        self.cldFrac                = np.nan    # average cloud amount
        #Method 4: additonal tuning parameters albedo as a funtion of age and water content (Bougamont et al., 2005)
        self.t0wet                  = np.nan    # time scale for wet snow (15-21.9)
        self.t0dry                  = np.nan    # warm snow timescale (30)
        self.K                      = np.nan    # time scale temperature coef. (7)
        self.adThresh               = np.nan    # Apply aIdx method to all areas with densities below this value,
        # or else apply direct input value from aValue, allowing albedo to be altered.
        # Default value is rho water (1023 kg m-3).
        teThresh                    = np.nan    #Apply eIdx method to all areas with grain radii above this value (mm),
        #or else apply direct input value from teValue, allowing emissivity to be altered.
        #Default value is a effective grain radius of 10 mm.

        # Densities
        self.InitDensityScaling     = np.nan    # initial scaling factor multiplying the density of ice, which describes the density of the snowpack.

        # Thermal
        self.ThermoDeltaTScaling    = np.nan    # scaling factor to multiply the thermal diffusion timestep (delta t)

        self.steps_per_step         = 1
        self.averaging              = 0
        self.requested_outputs      = []

        #Several fields are missing from the standard GEMB model, which are 
        #captured intrinsically by ISSM.
        #dateN: that's the last row of the above fields.
        #dt:    included in dateN. Not an input.
        #elev:  this is taken from the ISSM surface itself.

        nargin = len(args)
        if nargin == 1:
            mesh = args[0]
            self.setdefaultparameters(mesh)
        else:
            raise Exception('constructor not supported: need mesh to set defaults')
        # }}}

    def __repr__(self):  # {{{
        #string = "   surface forcings parameters:"
        #string = "#s\n#s"%(string, fielddisplay(self, 'mass_balance', 'surface mass balance [m/yr ice eq]'))
        #string = "#s\n#s"%(string, fielddisplay(self, 'requested_outputs', 'additional outputs requested'))
        s = '   surface forcings for SMB GEMB model :\n'
        s += '{}\n'.format(fielddisplay(self, 'isgraingrowth', 'run grain growth module (default true)'))
        s += '{}\n'.format(fielddisplay(self, 'isalbedo', 'run albedo module (default true)'))
        s += '{}\n'.format(fielddisplay(self, 'isshortwave', 'run short wave module (default true)'))
        s += '{}\n'.format(fielddisplay(self, 'isthermal', 'run thermal module (default true)'))
        s += '{}\n'.format(fielddisplay(self, 'isaccumulation', 'run accumulation module (default true)'))
        s += '{}\n'.format(fielddisplay(self, 'ismelt', 'run melting  module (default true)'))
        s += '{}\n'.format(fielddisplay(self, 'isdensification', 'run densification module (default true)'))
        s += '{}\n'.format(fielddisplay(self, 'isturbulentflux', 'run turbulant heat fluxes module (default true)'))
        s += '{}\n'.format(fielddisplay(self, 'isconstrainsurfaceT', 'constrain surface temperatures to air temperature, turn off EC and surface flux contribution to surface temperature change (default false)'))
        s += '{}\n'.format(fielddisplay(self, 'isdeltaLWup', 'set to true to invoke a bias in the long wave upward spatially, specified by dulwrfValue (default false)'))
        s += '{}\n'.format(fielddisplay(self,'ismappedforcing','set to true if forcing grid does not match model mesh, mapping specified by mappedforcingpoint (default false)'))
        s += '{}\n'.format(fielddisplay(self,'iscompressedforcing','set to true to compress the input matrices when writing to binary (default false)'))
        s += '{}\n'.format(fielddisplay(self, 'Ta', '2 m air temperature, in Kelvin'))
        s += '{}\n'.format(fielddisplay(self, 'V', 'wind speed (m s-1)'))
        s += '{}\n'.format(fielddisplay(self, 'dswrf', 'downward shortwave radiation flux [W/m^2]'))
        s += '{}\n'.format(fielddisplay(self, 'dswdiffrf', 'downward diffusive portion of shortwave radiation flux (default to 0) [W/m^2]'))
        s += '{}\n'.format(fielddisplay(self, 'dlwrf', 'downward longwave radiation flux [W/m^2]'))
        s += '{}\n'.format(fielddisplay(self, 'P', 'precipitation [mm w.e. / m^2]'))
        s += '{}\n'.format(fielddisplay(self, 'eAir', 'screen level vapor pressure [Pa]'))
        s += '{}\n'.format(fielddisplay(self, 'pAir', 'surface pressure [Pa]'))
        s += '{}\n'.format(fielddisplay(self, 'Tmean', 'mean annual temperature [K]'))
        s += '{}\n'.format(fielddisplay(self, 'C', 'mean annual snow accumulation [kg m-2 yr-1]'))
        s += '{}\n'.format(fielddisplay(self, 'Vmean', 'mean annual temperature [m s-1] (default 10 m/s)'))
        s += '{}\n'.format(fielddisplay(self, 'Tz', 'height above ground at which temperature (T) was sampled [m]'))
        s += '{}\n'.format(fielddisplay(self, 'Vz', 'height above ground at which wind (V) eas sampled [m]'))
        s += '{}\n'.format(fielddisplay(self, 'zTop', 'depth over which grid length is constant at the top of the snopack (default 10) [m]'))
        s += '{}\n'.format(fielddisplay(self, 'dzTop', 'initial top vertical grid spacing (default .05) [m] '))
        s += '{}\n'.format(fielddisplay(self, 'dzMin', 'initial min vertical allowable grid spacing (default dzMin/2) [m] '))
        s += '{}\n'.format(fielddisplay(self, 'zMax', 'initial max model depth (default is min(thickness, 500)) [m]'))
        s += '{}\n'.format(fielddisplay(self, 'zMin', 'initial min model depth (default is min(thickness, 30)) [m]'))
        s += '{}\n'.format(fielddisplay(self, 'zY', 'stretch grid cells bellow top_z by a [top_dz * y ^ (cells bellow top_z)]'))
        s += '{}\n'.format(fielddisplay(self, 'InitDensityScaling', ['initial scaling factor multiplying the density of ice', 'which describes the density of the snowpack.']))
        s += '{}\n'.format(fielddisplay(self, 'ThermoDeltaTScaling', 'scaling factor to multiply the thermal diffusion timestep (delta t)'))
        s += '{}\n'.format(fielddisplay(self, 'outputFreq', 'output frequency in days (default is monthly, 30)'))
        s += '{}\n'.format(fielddisplay(self, 'adThresh', 'Apply aIdx method to all areas with densities below this value, or else apply direct input value from aValue, allowing albedo to be altered.'))
        s += '{}\n'.format(fielddisplay(self, 'aIdx', ['method for calculating albedo and subsurface absorption (default is 1)',
            '0: direct input from aValue parameter',
            '1: effective grain radius [Gardner & Sharp, 2009]',
            '2: effective grain radius [Brun et al., 1992; LeFebre et al., 2003], with swIdx=1, SW penetration follows grain size in 3 spectral bands (Brun et al., 1992)',
            '3: density and cloud amount [Greuell & Konzelmann, 1994]',
            '4: exponential time decay & wetness [Bougamont & Bamber, 2005]']))

        s += '{}\n'.format(fielddisplay(self, 'dulwrfValue', 'Specified bias to be applied to the outward long wave radiation at every element (W/m-2, +upward)'))
        s += '{}\n'.format(fielddisplay(self, 'teValue', 'Outward longwave radiation thermal emissivity forcing at every element (default in code is 1)'))
        s += '{}\n'.format(fielddisplay(self, 'teThresh', ['Apply eIdx method to all areas with effective grain radius above this value (mm),', 'or else apply direct input value from teValue, allowing emissivity to be altered.']))
        s += '{}\n'.format(fielddisplay(self, 'eIdx', ['method for calculating emissivity (default is 1)',
            '0: direct input from teValue parameter, no use of teThresh',
            '1: default value of 1, in areas with grain radius below teThresh',
            '2: default value of 1, in areas with grain radius below teThresh and areas of dry snow (not bare ice or wet) at the surface']))
        s += '{}\n'.format(fielddisplay(self, 'tcIdx', ['method for calculating thermal conductivity (default is 1)',
            '1: after Sturm et al, 1997',
            '2: after Calonne et al., 2011']))

        s += '{}\n'.format(fielddisplay(self,'mappedforcingpoint','Mapping of which forcing point will map to each mesh element for ismappedforcing option (integer). Size number of elements.'))
        s += '{}\n'.format(fielddisplay(self,'mappedforcingelevation','The elevation of each mapped forcing location (m above sea level) for ismappedforcing option. Size number of forcing points.'))
        s += '{}\n'.format(fielddisplay(self,'lapseTaValue','Temperature lapse rate if forcing has different grid and should be remapped for ismappedforcing option. (Default value is -0.006 K m-1.)'))
        s += '{}\n'.format(fielddisplay(self,'lapsedlwrfValue','Longwave down lapse rate if forcing has different grid and should be remapped for ismappedforcing option. (Default value is -0.032 W m-2 m-1.)'))

        # Snow properties init
        s += '{}\n'.format(fielddisplay(self, 'Dzini', 'Initial cell depth when restart [m]'))
        s += '{}\n'.format(fielddisplay(self, 'Dini', 'Initial snow density when restart [kg m-3]'))
        s += '{}\n'.format(fielddisplay(self, 'Reini', 'Initial grain size when restart [mm]'))
        s += '{}\n'.format(fielddisplay(self, 'Gdnini', 'Initial grain dricity when restart [-]'))
        s += '{}\n'.format(fielddisplay(self, 'Gspini', 'Initial grain sphericity when restart [-]'))
        s += '{}\n'.format(fielddisplay(self, 'ECini', 'Initial evaporation/condensation when restart [kg m-2]'))
        s += '{}\n'.format(fielddisplay(self, 'Wini', 'Initial snow water content when restart [kg m-2]'))
        s += '{}\n'.format(fielddisplay(self, 'Aini', 'Initial albedo when restart [-]'))
        s += '{}\n'.format(fielddisplay(self, 'Adiffini', 'Initial diffusive radiation albedo when restart (default to 1) [-]'))
        s += '{}\n'.format(fielddisplay(self, 'Tini', 'Initial snow temperature when restart [K]'))
        s += '{}\n'.format(fielddisplay(self, 'Sizeini', 'Initial number of layers when restart [-]'))

        # Additional albedo parameters
        s += '{}\n'.format(fielddisplay(self, 'aValue', 'Albedo forcing at every element'))
        if isinstance(self.aIdx, (list, type(np.array([1, 2])))) and (self.aIdx == [1, 2] or (1 in self.aIdx and 2 in self.aIdx)):
            s += '{}\n'.format(fielddisplay(self, 'aSnow', 'new snow albedo (0.64 - 0.89)'))
            s += '{}\n'.format(fielddisplay(self, 'aIce', 'albedo of ice (0.27-0.58)'))
            if self.aIdx == 1:
                s += '{}\n'.format(fielddisplay(self,'szaValue','Solar Zenith Angle [degree]'))
                s += '{}\n'.format(fielddisplay(self,'cotValue','Cloud Optical Thickness'))
                s += '{}\n'.format(fielddisplay(self,'ccsnowValue','concentration of light absorbing carbon for snow [ppm1]'))
                s += '{}\n'.format(fielddisplay(self,'cciceValue','concentration of light absorbing carbon for ice [ppm1]'))
        elif self.aIdx == 3:
            s += '{}\n'.format(fielddisplay(self, 'cldFrac', 'average cloud amount'))
        elif self.aIdx == 4:
            s += '{}\n'.format(fielddisplay(self, 't0wet', 'time scale for wet snow (15-21.9) [d]'))
            s += '{}\n'.format(fielddisplay(self, 't0dry', 'warm snow timescale (30) [d]'))
            s += '{}\n'.format(fielddisplay(self, 'K', 'time scale temperature coef. (7) [d]'))

        s += '{}\n'.format(fielddisplay(self, 'swIdx', 'apply all SW to top grid cell (0) or allow SW to penetrate surface (1) [default 0, if swIdx=1 and aIdx=2 function of effective radius (Brun et al., 1992) or else dependent on snow density (taken from Bassford, 2002)]'))
        s += '{}\n'.format(fielddisplay(self, 'denIdx', ['densification model to use (default is 2):',
            '1 = emperical model of Herron and Langway (1980)',
            '2 = semi-emperical model of Anthern et al. (2010)',
            '3 = DO NOT USE: physical model from Appix B of Anthern et al. (2010)',
            '4 = DO NOT USE: emperical model of Li and Zwally (2004)',
            '5 = DO NOT USE: modified emperical model (4) by Helsen et al. (2008)',
            '6 = Antarctica semi-emperical model of Ligtenberg et al. (2011)',
            '7 = Greenland semi-emperical model of Kuipers Munneke et al. (2015)']))
        s += '{}\n'.format(fielddisplay(self, 'dsnowIdx', ['model for fresh snow accumulation density (default is 1):',
            '0 = Original GEMB value, 150 kg/m^3',
            '1 = Antarctica value of fresh snow density, 350 kg/m^3',
            '2 = Greenland value of fresh snow density, 315 kg/m^3, Fausto et al. (2018)',
            '3 = Antarctica model of Kaspers et al. (2004), Make sure to set Vmean accurately',
            '4 = Greenland model of Kuipers Munneke et al. (2015)']))

        s += '{}\n'.format(fielddisplay(self, 'steps_per_step', 'number of smb steps per time step'))
        s += '{}\n'.format(fielddisplay(self, 'averaging', 'averaging methods from short to long steps'))
        s += '\t\t{}\n'.format('0: Arithmetic (default)')
        s += '\t\t{}\n'.format('1: Geometric')
        s += '\t\t{}\n'.format('2: Harmonic')
        s += '{}\n'.format(fielddisplay(self, 'requested_outputs', 'additional outputs requested'))
        return s
    # }}}

    def extrude(self, md):  # {{{
        if np.shape(self.Ta)[0] == md.mesh.numberofelements or np.shape(self.Ta)[0] == md.mesh.numberofelements + 1 :
            self.Ta = project3d(md, 'vector', self.Ta, 'type', 'element')
            self.V = project3d(md, 'vector', self.V, 'type', 'element')
            self.dswrf = project3d(md, 'vector', self.dswrf, 'type', 'element')
            self.dlwrf = project3d(md, 'vector', self.dlwrf, 'type', 'element')
            self.P = project3d(md, 'vector', self.P, 'type', 'element')
            self.eAir = project3d(md, 'vector', self.eAir, 'type', 'element')
            self.pAir = project3d(md, 'vector', self.pAir, 'type', 'element')

        if not np.isnan(self.Dzini):
            self.self.Dzini=project3d(md,'vector',self.self.Dzini,'type','element');
        if not np.isnan(self.Dini):
            self.self.Dini=project3d(md,'vector',self.Dini,'type','element');
        if not np.isnan(self.Reini):
            self.self.Reini=project3d(md,'vector',self.Reini,'type','element');
        if not np.isnan(self.Gdnini):
            self.Gdnini=project3d(md,'vector',self.Gdnini,'type','element');
        if not np.isnan(self.Gspini):
            self.Gspini=project3d(md,'vector',self.Gspini,'type','element');
        if not np.isnan(self.ECini):
            self.ECini=project3d(md,'vector',self.ECini,'type','element');
        if not np.isnan(self.Wini):
            self.Wini=project3d(md,'vector',self.Wini,'type','element');
        if not np.isnan(self.Aini):
            self.Aini=project3d(md,'vector',self.Aini,'type','element');
        if not np.isnan(self.Adiffini):
            self.Adiffini=project3d(md,'vector',self.Adiffini,'type','element');
        if not np.isnan(self.Tini):
            self.Tini=project3d(md,'vector',self.Tini,'type','element');

        if not np.isnan(self.dswdiffrf):
            self.dswdiffrf=project3d(md,'vector',self.dswdiffrf,'type','element');
        if not np.isnan(self.szaValue):
            self.szaValue=project3d(md,'vector',self.szaValue,'type','element');
        if not np.isnan(self.cotValue):
            self.cotValue=project3d(md,'vector',self.cotValue,'type','element');
        if not np.isnan(self.ccsnowValue):
            self.ccsnowValue=project3d(md,'vector',self.ccsnowValue,'type','element');
        if not np.isnan(self.cciceValue):
            self.cciceValue=project3d(md,'vector',self.cciceValue,'type','element');

        if not np.isnan(self.aValue):
            self.aValue = project3d(md, 'vector', self.aValue, 'type', 'element')
        if not np.isnan(self.teValue):
            self.teValue = project3d(md, 'vector', self.teValue, 'type', 'element')
        if not np.isnan(self.mappedforcingpoint):
            self.mappedforcingpoint = project3d(md, 'vector', self.mappedforcingpoint, 'type', 'element')

        return self
    # }}}

    def defaultoutputs(self, md):  # {{{
        return ['SmbMassBalance','SmbAccumulatedMassBalance']
    # }}}

    def setdefaultparameters(self, mesh):  # {{{
        self.isgraingrowth = 1
        self.isalbedo = 1
        self.isshortwave = 1
        self.isthermal = 1
        self.isaccumulation = 1
        self.ismelt = 1
        self.isdensification = 1
        self.isturbulentflux = 1
        self.isconstrainsurfaceT = 0
        self.isdeltaLWup = 0
        self.ismappedforcing = 0
        self.iscompressedforcing = 0

        self.aIdx = 1
        self.eIdx = 1
        self.tcIdx = 1
        self.swIdx = 0 
        self.denIdx = 2
        self.dsnowIdx = 1
        self.zTop = 10 * np.ones((mesh.numberofelements,))
        self.dzTop = 0.05 * np.ones((mesh.numberofelements,))
        self.dzMin = self.dzTop / 2
        self.InitDensityScaling = 1.0
        self.ThermoDeltaTScaling = 1 / 11.0

        self.Vmean = 10 * np.ones((mesh.numberofelements,))

        self.zMax = 250 * np.ones((mesh.numberofelements,))
        self.zMin = 130 * np.ones((mesh.numberofelements,))
        self.zY = 1.025 * np.ones((mesh.numberofelements,))
        self.outputFreq = 30

        # Additional albedo parameters
        self.aSnow = 0.85
        self.aIce = 0.48
        self.cldFrac = 0.1
        self.t0wet = 15
        self.t0dry = 30
        self.K = 7
        self.adThresh = 1023
        self.teThresh = 10

        self.teValue = np.ones((mesh.numberofelements,))
        self.aValue = self.aSnow * np.ones(mesh.numberofelements,)
        self.dulwrfValue = np.zeros((mesh.numberofelements,))
        self.lapseTaValue = -0.006
        self.lapsedlwrfValue = -0.032

        self.dswdiffrf = 0.0 * np.ones(mesh.numberofelements,)
        self.szaValue = 0.0 * np.ones(mesh.numberofelements,)
        self.cotValue = 0.0 * np.ones(mesh.numberofelements,)
        self.ccsnowValue = 0.0 * np.ones(mesh.numberofelements,)
        self.cciceValue = 0.0 * np.ones(mesh.numberofelements,)

        self.Dzini = 0.05 * np.ones((mesh.numberofelements, 2))
        self.Dini = 910.0 * np.ones((mesh.numberofelements, 2))
        self.Reini = 2.5 * np.ones((mesh.numberofelements, 2))
        self.Gdnini = 0.0 * np.ones((mesh.numberofelements, 2))
        self.Gspini = 0.0 * np.ones((mesh.numberofelements, 2))
        self.ECini = 0.0 * np.ones((mesh.numberofelements,))
        self.Wini = 0.0 * np.ones((mesh.numberofelements, 2))
        self.Aini = 0.0 * np.ones((mesh.numberofelements, 2))
        self.Adiffini = np.ones((mesh.numberofelements, 2))
        self.Tini = 273.15 * np.ones((mesh.numberofelements, 2))
        #       /!\ Default value of Tini must be equal to Tmean but don't know 
        #           Tmean yet (computed when atmospheric forcings are 
        #           interpolated on mesh).
        #           If initialization without restart, this value will be 
        #           overwritten when snow parameters are retrieved in 
        #           Element.cpp
        self.Sizeini = 2 * np.ones((mesh.numberofelements,))
    # }}}

    def checkconsistency(self, md, solution, analyses):    # {{{

        md = checkfield(md, 'fieldname', 'smb.isgraingrowth', 'values', [0, 1])
        md = checkfield(md, 'fieldname', 'smb.isalbedo', 'values', [0, 1])
        md = checkfield(md, 'fieldname', 'smb.isshortwave', 'values', [0, 1])
        md = checkfield(md, 'fieldname', 'smb.isthermal', 'values', [0, 1])
        md = checkfield(md, 'fieldname', 'smb.isaccumulation', 'values', [0, 1])
        md = checkfield(md, 'fieldname', 'smb.ismelt', 'values', [0, 1])
        md = checkfield(md, 'fieldname', 'smb.isdensification', 'values', [0, 1])
        md = checkfield(md, 'fieldname', 'smb.isturbulentflux', 'values', [0, 1])
        md = checkfield(md, 'fieldname', 'smb.isdeltaLWup', 'values',[0, 1])
        md = checkfield(md, 'fieldname', 'smb.isconstrainsurfaceT', 'values', [0, 1])
        md = checkfield(md, 'fieldname', 'smb.ismappedforcing', 'values',[0, 1])
        md = checkfield(md, 'fieldname', 'smb.iscompressedforcing', 'values',[0, 1])

        sizeta=np.shape(self.Ta)
        md = checkfield(md, 'fieldname', 'smb.Ta', 'mappedtimeseries', 1, 'NaN', 1, 'Inf', 1, '>', 273-100, '<', 273+100) #-100/100 celsius min/max value
        md = checkfield(md, 'fieldname', 'smb.V', 'mappedtimeseries', 1, 'NaN', 1, 'Inf', 1, '>=', 0, '<', 45, 'size', sizeta) #max 500 km/h
        md = checkfield(md, 'fieldname', 'smb.dswrf', 'mappedtimeseries', 1, 'NaN', 1, 'Inf', 1, '>=', 0, '<=', 1400, 'size', sizeta)
        md = checkfield(md, 'fieldname', 'smb.dswdiffrf', 'mappedtimeseries', 1, 'NaN', 1, 'Inf', 1, '>=', 0, '<=', 1400)
        md = checkfield(md, 'fieldname', 'smb.dlwrf', 'mappedtimeseries', 1, 'NaN', 1, 'Inf', 1, '>=', 0, 'size', sizeta)
        md = checkfield(md, 'fieldname', 'smb.P', 'mappedtimeseries', 1, 'NaN', 1, 'Inf', 1, '>=', 0, '<=', 200, 'size', sizeta)
        md = checkfield(md, 'fieldname', 'smb.eAir', 'mappedtimeseries',1,'NaN',1,'Inf', 1, 'size', sizeta)

        md = checkfield(md, 'fieldname', 'smb.Tmean', 'size', [sizeta[0]-1], 'NaN', 1, 'Inf', 1, '>', 273-100, '<', 273+100) #-100/100 celsius min/max value
        md = checkfield(md, 'fieldname', 'smb.C', 'size', [sizeta[0]-1], 'NaN', 1, 'Inf', 1, '>', 0)
        md = checkfield(md, 'fieldname', 'smb.Vmean', 'size', [sizeta[0]-1], 'NaN', 1, 'Inf', 1, '>=', 0)
        md = checkfield(md, 'fieldname', 'smb.Tz', 'size', [sizeta[0]-1], 'NaN', 1, 'Inf', 1, '>=', 0, '<=', 5000)
        md = checkfield(md, 'fieldname', 'smb.Vz', 'size', [sizeta[0]-1], 'NaN', 1, 'Inf', 1, '>=', 0, '<=', 5000)

        md = checkfield(md, 'fieldname', 'smb.teValue', 'timeseries', 1, 'NaN', 1, 'Inf', 1, '>=', 0, '<=', 1)
        md = checkfield(md, 'fieldname', 'smb.dulwrfValue', 'timeseries', 1, 'NaN', 1, 'Inf', 1)

        if self.ismappedforcing:
            md = checkfield(md, 'fieldname', 'smb.mappedforcingpoint', 'size',[md.mesh.numberofelements], 'NaN', 1, 'Inf', 1, '>', 0, '<=' ,sizeta[0]-1)
            md = checkfield(md, 'fieldname', 'smb.mappedforcingelevation', 'size', [sizeta[0]-1], 'NaN', 1, 'Inf', 1)
            md = checkfield(md, 'fieldname', 'smb.lapseTaValue', 'NaN', 1, 'Inf', 1)
            md = checkfield(md, 'fieldname', 'smb.lapsedlwrfValue', 'NaN', 1, 'Inf', 1)

        md = checkfield(md, 'fieldname', 'smb.aIdx', 'NaN', 1, 'Inf', 1, 'values', [0, 1, 2, 3, 4])
        md = checkfield(md, 'fieldname', 'smb.eIdx', 'NaN', 1, 'Inf', 1, 'values', [0, 1, 2])
        md = checkfield(md, 'fieldname', 'smb.tcIdx', 'NaN', 1, 'Inf', 1, 'values', [1, 2])
        md = checkfield(md, 'fieldname', 'smb.swIdx', 'NaN', 1, 'Inf', 1, 'values', [0, 1])
        md = checkfield(md, 'fieldname', 'smb.denIdx', 'NaN', 1, 'Inf', 1, 'values', [1, 2, 3, 4, 5, 6, 7])
        md = checkfield(md, 'fieldname', 'smb.dsnowIdx', 'NaN', 1, 'Inf', 1, 'values', [0, 1, 2, 3, 4])

        md = checkfield(md, 'fieldname', 'smb.zTop', 'NaN', 1, 'Inf', 1, '> = ', 0)
        md = checkfield(md, 'fieldname', 'smb.dzTop', 'NaN', 1, 'Inf', 1, '>', 0)
        md = checkfield(md, 'fieldname', 'smb.dzMin', 'NaN', 1, 'Inf', 1, '>', 0)
        md = checkfield(md, 'fieldname', 'smb.zY', 'NaN', 1, 'Inf', 1, '> = ', 1)
        md = checkfield(md, 'fieldname', 'smb.outputFreq', 'NaN', 1, 'Inf', 1, '>', 0, '<', 10 * 365)  #10 years max
        md = checkfield(md, 'fieldname', 'smb.InitDensityScaling', 'NaN', 1, 'Inf', 1, '> = ', 0, '< = ', 1)
        md = checkfield(md, 'fieldname', 'smb.ThermoDeltaTScaling', 'NaN', 1, 'Inf', 1, '> = ', 0, '< = ', 1)
        md = checkfield(md, 'fieldname', 'smb.adThresh', 'NaN', 1, 'Inf', 1, '>=', 0)
        md = checkfield(md, 'fieldname', 'smb.teThresh', 'NaN', 1, 'Inf',1,'>=',0)
        
        md = checkfield(md, 'fieldname', 'smb.aValue', 'timeseries', 1, 'NaN', 1, 'Inf', 1, '>=', 0, '<=', 1)
        if isinstance(self.aIdx, (list, type(np.array([1, 2])))) and (self.aIdx == [1, 2] or (1 in self.aIdx and 2 in self.aIdx)):
            md = checkfield(md, 'fieldname', 'smb.aSnow', 'NaN', 1, 'Inf', 1, '> = ', .64, '< = ', .89)
            md = checkfield(md, 'fieldname', 'smb.aIce', 'NaN', 1, 'Inf', 1, '> = ', .27, '< = ', .58)
            if self.aIdx == 1:
                md = checkfield(md,'fieldname','smb.szaValue','timeseries',1,'NaN',1,'Inf',1,'>=',0,'<=',90)
                md = checkfield(md,'fieldname','smb.cotValue','timeseries',1,'NaN',1,'Inf',1,'>=',0)
                md = checkfield(md,'fieldname','smb.ccsnowValue','timeseries',1,'NaN',1,'Inf',1,'>=',0)
                md = checkfield(md,'fieldname','smb.cciceValue','timeseries',1,'NaN',1,'Inf',1,'>=',0)
        elif self.aIdx == 3:
            md = checkfield(md, 'fieldname', 'smb.cldFrac', 'NaN', 1, 'Inf', 1, '> = ', 0, '< = ', 1)
        elif self.aIdx == 4:
            md = checkfield(md, 'fieldname', 'smb.t0wet', 'NaN', 1, 'Inf', 1, '> = ', 15, '< = ', 21.9)
            md = checkfield(md, 'fieldname', 'smb.t0dry', 'NaN', 1, 'Inf', 1, '> = ', 30, '< = ', 30)
            md = checkfield(md, 'fieldname', 'smb.K', 'NaN', 1, 'Inf', 1, '> = ', 7, '< = ', 7)

        # Check zTop is < local thickness
        he = np.sum(md.geometry.thickness[md.mesh.elements - 1], axis=1) / np.size(md.mesh.elements, 1)
        if np.any(he < self.zTop):
            raise IOError('SMBgemb consistency check error: zTop should be smaller than local ice thickness')
        md = checkfield(md, 'fieldname', 'smb.steps_per_step', '>=', 1, 'numel', [1])
        md = checkfield(md, 'fieldname', 'smb.averaging', 'numel', [1], 'values', [0, 1, 2])
        md = checkfield(md, 'fieldname', 'smb.requested_outputs', 'stringrow', 1)
        return md
    # }}}

    def marshall(self, prefix, md, fid):    # {{{
        yts = md.constants.yts

        WriteData(fid, prefix, 'name', 'md.smb.model', 'data', 8, 'format', 'Integer')

        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'isgraingrowth', 'format', 'Boolean')
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'isalbedo', 'format', 'Boolean')
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'isshortwave', 'format', 'Boolean')
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'isthermal', 'format', 'Boolean')
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'isaccumulation', 'format', 'Boolean')
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'ismelt', 'format', 'Boolean')
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'isdensification', 'format', 'Boolean')
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'isturbulentflux', 'format', 'Boolean')
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'isconstrainsurfaceT', 'format', 'Boolean')
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'isdeltaLWup', 'format', 'Boolean')
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'ismappedforcing', 'format', 'Boolean')

        if self.iscompressedforcing:
            writetype='CompressedMat'
        else:
            writetype='DoubleMat'

        WriteData(fid,prefix, 'object', self, 'class', 'smb', 'fieldname', 'Ta', 'format', writetype, 'mattype', 2, 'timeserieslength', np.shape(self.Ta)[0], 'yts', md.constants.yts)
        WriteData(fid,prefix, 'object', self, 'class', 'smb', 'fieldname', 'V', 'format', writetype, 'mattype', 2, 'timeserieslength', np.shape(self.Ta)[0], 'yts', md.constants.yts)
        WriteData(fid,prefix, 'object', self, 'class', 'smb', 'fieldname', 'dswrf', 'format', writetype, 'mattype', 2, 'timeserieslength', np.shape(self.Ta)[0], 'yts', md.constants.yts)
        WriteData(fid,prefix, 'object', self, 'class', 'smb', 'fieldname', 'dswdiffrf', 'format', writetype, 'mattype', 2, 'timeserieslength', np.shape(self.Ta)[0], 'yts', md.constants.yts)
        WriteData(fid,prefix, 'object', self, 'class', 'smb', 'fieldname', 'dlwrf', 'format', writetype, 'mattype', 2, 'timeserieslength', np.shape(self.Ta)[0], 'yts', md.constants.yts)
        WriteData(fid,prefix, 'object', self, 'class', 'smb', 'fieldname', 'P', 'format', writetype, 'mattype', 2, 'timeserieslength', np.shape(self.Ta)[0], 'yts', md.constants.yts)
        WriteData(fid,prefix, 'object', self, 'class', 'smb', 'fieldname', 'eAir', 'format', writetype, 'mattype', 2, 'timeserieslength', np.shape(self.Ta)[0], 'yts', md.constants.yts)
        WriteData(fid,prefix, 'object', self, 'class', 'smb', 'fieldname', 'pAir', 'format', writetype, 'mattype', 2, 'timeserieslength', np.shape(self.Ta)[0], 'yts', md.constants.yts)

        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'Tmean', 'format', 'DoubleMat', 'mattype', 2)
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'C', 'format', 'DoubleMat', 'mattype', 2)
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'Vmean', 'format', 'DoubleMat', 'mattype', 2)
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'Tz', 'format', 'DoubleMat', 'mattype', 2)
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'Vz', 'format', 'DoubleMat', 'mattype', 2)
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'zTop', 'format', 'DoubleMat', 'mattype', 2)
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'dzTop', 'format', 'DoubleMat', 'mattype', 2)
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'dzMin', 'format', 'DoubleMat', 'mattype', 2)
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'zY', 'format', 'DoubleMat', 'mattype', 2)
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'zMax', 'format', 'DoubleMat', 'mattype', 2)
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'zMin', 'format', 'DoubleMat', 'mattype', 2)

        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'aIdx', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'eIdx', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'tcIdx', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'swIdx', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'denIdx', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'dsnowIdx', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'InitDensityScaling', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'ThermoDeltaTScaling', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'outputFreq', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'aSnow', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'aIce', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'cldFrac', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 't0wet', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 't0dry', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'K', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'adThresh', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'teThresh', 'format', 'Double')

        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'aValue', 'format', 'DoubleMat', 'mattype', 2, 'timeserieslength', md.mesh.numberofelements + 1, 'yts', yts)
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'teValue', 'format', 'DoubleMat', 'mattype', 2, 'timeserieslength', md.mesh.numberofelements + 1, 'yts', yts)
        WriteData(fid,prefix,'object',self,'class','smb','fieldname','dulwrfValue','format','DoubleMat','mattype',2,'timeserieslength',md.mesh.numberofelements+1,'yts',md.constants.yts)
        WriteData(fid,prefix,'object',self,'class','smb','fieldname','szaValue','format','DoubleMat','mattype',2,'timeserieslength',md.mesh.numberofelements+1,'yts',md.constants.yts)
        WriteData(fid,prefix,'object',self,'class','smb','fieldname','cotValue','format','DoubleMat','mattype',2,'timeserieslength',md.mesh.numberofelements+1,'yts',md.constants.yts)
        WriteData(fid,prefix,'object',self,'class','smb','fieldname','ccsnowValue','format','DoubleMat','mattype',2,'timeserieslength',md.mesh.numberofelements+1,'yts',md.constants.yts)
        WriteData(fid,prefix,'object',self,'class','smb','fieldname','cciceValue','format','DoubleMat','mattype',2,'timeserieslength',md.mesh.numberofelements+1,'yts',md.constants.yts)

        #snow properties init
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'Dzini', 'format', 'DoubleMat', 'mattype', 3)
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'Dini', 'format', 'DoubleMat', 'mattype', 3)
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'Reini', 'format', 'DoubleMat', 'mattype', 3)
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'Gdnini', 'format', 'DoubleMat', 'mattype', 3)
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'Gspini', 'format', 'DoubleMat', 'mattype', 3)
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'ECini', 'format', 'DoubleMat', 'mattype', 2)
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'Wini', 'format', 'DoubleMat', 'mattype', 3)
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'Aini', 'format', 'DoubleMat', 'mattype', 3)
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'Adiffini', 'format', 'DoubleMat', 'mattype', 3)
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'Tini', 'format', 'DoubleMat', 'mattype', 3)
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'Sizeini', 'format', 'IntMat', 'mattype', 2)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'steps_per_step', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'averaging', 'format', 'Integer')

        if self.ismappedforcing:
            WriteData(fid,prefix,'object',self,'class','smb','fieldname','mappedforcingpoint','format','IntMat','mattype',2)
            WriteData(fid,prefix,'object',self,'class','smb','fieldname','mappedforcingelevation','format','DoubleMat','mattype',3)
            WriteData(fid,prefix,'object',self,'class','smb','fieldname','lapseTaValue','format','Double')
            WriteData(fid,prefix,'object',self,'class','smb','fieldname','lapsedlwrfValue','format','Double')

        # Figure out dt from forcings
        if (np.any(self.P[-1] - self.Ta[-1] != 0) | np.any(self.V[-1] - self.Ta[-1] != 0) | np.any(self.dswrf[-1] - self.Ta[-1] != 0) | np.any(self.dlwrf[-1] - self.Ta[-1] != 0) | np.any(self.eAir[-1] - self.Ta[-1] != 0) | np.any(self.pAir[-1] - self.Ta[-1] != 0)):
            raise IOError('All GEMB forcings (Ta, P, V, dswrf, dlwrf, eAir, pAir) must have the same time steps in the final row!')

        if ((np.ndim(self.teValue)>1) & np.any(self.teValue[-1] - self.Ta[-1] != 0)):
            raise IOError('If GEMB forcing teValue is transient, it must have the same time steps as input Ta in the final row!')
        if ((np.ndim(self.dswdiffrf)>1) & np.any(self.dswdiffrf[-1] - self.Ta[-1] != 0)):
            raise IOError('If GEMB forcing dswdiffrf is transient, it must have the same time steps as input Ta in the final row!')
        if ((np.ndim(self.aValue)>1) & np.any(self.aValue[-1] - self.Ta[-1] != 0)):
            raise IOError('If GEMB forcing aValue is transient, it must have the same time steps as input Ta in the final row!')
        if ((np.ndim(self.dulwrfValue)>1) & np.any(self.dulwrfValue[-1] - self.Ta[-1] != 0)):
            raise IOError('If GEMB forcing dulwrfValue is transient, it must have the same time steps as input Ta in the final row!')
        if ((np.ndim(self.szaValue)>1) & np.any(self.szaValue[-1] - self.Ta[-1] != 0)):
            raise IOError('If GEMB forcing szaValue is transient, it must have the same time steps as input Ta in the final row!')
        if ((np.ndim(self.cotValue)>1) & np.any(self.cotValue[-1] - self.Ta[-1] != 0)):
            raise IOError('If GEMB forcing cotValue is transient, it must have the same time steps as input Ta in the final row!')
        if ((np.ndim(self.ccsnowValue)>1) & np.any(self.ccsnowValue[-1] - self.Ta[-1] != 0)):
            raise IOError('If GEMB forcing ccsnowValue is transient, it must have the same time steps as input Ta in the final row!')
        if ((np.ndim(self.cciceValue)>1) & np.any(self.cciceValue[-1] - self.Ta[-1] != 0)):
            raise IOError('If GEMB forcing cciceValue is transient, it must have the same time steps as input Ta in the final row!')

        time = self.Ta[-1]  # Assume all forcings are on the same time step
        dtime = np.diff(time, n=1, axis=0)
        dt = min(dtime)

        WriteData(fid, prefix, 'data', dt, 'name', 'md.smb.dt', 'format', 'Double', 'scale', yts)

        # Check if smb_dt goes evenly into transient core time step
        if (md.timestepping.time_step % dt >= 1e-10):
            raise IOError('smb_dt/dt = {}. The number of SMB time steps in one transient core time step has to be an an integer'.format(md.timestepping.time_step / dt))
        # Make sure that adaptive time step is off
        if md.timestepping.__class__.__name__ == 'timesteppingadaptive':
            raise IOError('GEMB cannot be run with adaptive timestepping.  Check class type of md.timestepping')

        # Process requested outputs
        outputs = self.requested_outputs
        indices = [i for i, x in enumerate(outputs) if x == 'default']
        if len(indices) > 0:
            outputscopy = outputs[0:max(0, indices[0] - 1)] + self.defaultoutputs(md) + outputs[indices[0] + 1:]
            outputs = outputscopy

        WriteData(fid, prefix, 'data', outputs, 'name', 'md.smb.requested_outputs', 'format', 'StringArray')
    # }}}
