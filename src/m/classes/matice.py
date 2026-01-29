import numpy as np

from fielddisplay import fielddisplay
from project3d import project3d
from checkfield import checkfield
from WriteData import WriteData


class matice(object):
    """MATICE class definition

    Usage:
            matice = matice()
    """

    def __init__(self):  # {{{
        self.rho_ice = 0
        self.rho_water = 0
        self.rho_freshwater = 0
        self.mu_water = 0
        self.heatcapacity = 0
        self.latentheat = 0
        self.thermalconductivity = 0
        self.temperateiceconductivity = 0
        self.effectiveconductivity_averaging = 0
        self.meltingpoint = 0
        self.beta = 0
        self.mixed_layer_capacity = 0
        self.thermal_exchange_velocity = 0
        self.rheology_B = np.nan
        self.rheology_n = np.nan
        self.rheology_law = ''

        # SLC
        self.earth_density = 0

        self.setdefaultparameters()
    # }}}

    def __repr__(self):  # {{{
        s = '   Materials:\n'
        s += '{}\n'.format(fielddisplay(self, 'rho_ice', 'ice density [kg/m^3]'))
        s += '{}\n'.format(fielddisplay(self, 'rho_water', 'water density [kg/m^3]'))
        s += '{}\n'.format(fielddisplay(self, 'rho_freshwater', 'fresh water density [kg/m^3]'))
        s += '{}\n'.format(fielddisplay(self, 'mu_water', 'water viscosity [Ns/m^2]'))
        s += '{}\n'.format(fielddisplay(self, 'heatcapacity', 'heat capacity [J/kg/K]'))
        s += '{}\n'.format(fielddisplay(self, 'thermalconductivity', 'ice thermal conductivity [W/m/K]'))
        s += '{}\n'.format(fielddisplay(self, 'temperateiceconductivity', 'temperate ice thermal conductivity [W/m/K]'))
        s += '{}\n'.format(fielddisplay(self, 'effectiveconductivity_averaging', 'computation of effectiveconductivity: (0) arithmetic mean, (1) harmonic mean, (2) geometric mean (default)'))
        s += '{}\n'.format(fielddisplay(self, 'meltingpoint', 'melting point of ice at 1atm in K'))
        s += '{}\n'.format(fielddisplay(self, 'latentheat', 'latent heat of fusion [J/m^3]'))
        s += '{}\n'.format(fielddisplay(self, 'beta', 'rate of change of melting point with pressure [K/Pa]'))
        s += '{}\n'.format(fielddisplay(self, 'mixed_layer_capacity', 'mixed layer capacity [W/kg/K]'))
        s += '{}\n'.format(fielddisplay(self, 'thermal_exchange_velocity', 'thermal exchange velocity [m/s]'))
        s += '{}\n'.format(fielddisplay(self, 'rheology_B', 'flow law parameter [Pa s^(1/n)]'))
        s += '{}\n'.format(fielddisplay(self, 'rheology_n', 'Glen\'s flow law exponent'))
        s += '{}\n'.format(fielddisplay(self, 'rheology_law', 'law for the temperature dependance of the rheology: \'None\', \'BuddJacka\', \'Cuffey\', \'CuffeyTemperate\', \'Paterson\', \'Arrhenius\', \'LliboutryDuval\', \'NyeCO2\', or \'NyeH2O\''))
        s += '{}\n'.format(fielddisplay(self, 'earth_density', 'Mantle density [kg/m^-3]'))
        return s
    # }}}

    def extrude(self, md):  # {{{
        self.rheology_B = project3d(md, 'vector', self.rheology_B, 'type', 'node')
        self.rheology_n = project3d(md, 'vector', self.rheology_n, 'type', 'element')
        return self
    # }}}

    def setdefaultparameters(self):  # {{{
        # Ice density (kg/m^3)
        self.rho_ice = 917.0
        # Ocean water density (kg/m^3)
        self.rho_water = 1023.0
        # Fresh water density (kg/m^3)
        self.rho_freshwater = 1000.0
        # Water viscosity (N.s/m^2)
        self.mu_water = 0.001787
        # Ice heat capacity cp (J/kg/K)
        self.heatcapacity = 2093.0
        # Ice latent heat of fusion L (J/kg)
        self.latentheat = 3.34 * pow(10, 5)
        # Ice thermal conductivity (W/m/K)
        self.thermalconductivity = 2.4
        # Temperate ice thermal conductivity (W/m/K)
        self.temperateiceconductivity = 0.24
        # Computation of effective conductivity
        self.effectiveconductivity_averaging = 1
        # The melting point of ice at 1 atmosphere of pressure in K
        self.meltingpoint = 273.15
        # Rate of change of melting point with pressure (K/Pa)
        self.beta = 9.8 * pow(10, -8)
        # Mixed layer (ice-water interface) heat capacity (J/kg/K)
        self.mixed_layer_capacity = 3974.0
        # Thermal exchange velocity (ice-water interface) (m/s)
        self.thermal_exchange_velocity = 1.00 * pow(10, -4)
        # Rheology law: what is the temperature dependence of B with T
        # available: none, paterson and arrhenius
        self.rheology_law = 'Paterson'

        # Rheology for ice
        self.rheology_B = 2.1 * 1e8
        self.rheology_n = 3

        # SLR
        self.earth_density = 5512  # average density of the Earth, (kg/m^3)
    # }}}

    def checkconsistency(self, md, solution, analyses):  # {{{
        if solution == 'TransientSolution' and md.transient.isslc:
            md = checkfield(md, 'fieldname', 'materials.earth_density', '>', 0, 'numel', [1])
        else:
            md = checkfield(md, 'fieldname', 'materials.rho_ice', '>', 0)
            md = checkfield(md, 'fieldname', 'materials.rho_water', '>', 0)
            md = checkfield(md, 'fieldname', 'materials.rho_freshwater', '>', 0)
            md = checkfield(md, 'fieldname', 'materials.mu_water', '>', 0)
            md = checkfield(md, 'fieldname', 'materials.rheology_B', '>', 0, 'size', 'universal', 'NaN', 1, 'Inf', 1)
            md = checkfield(md, 'fieldname', 'materials.rheology_n', '>', 0, 'size', 'universal', 'NaN', 1, 'Inf', 1)
            md = checkfield(md, 'fieldname', 'materials.rheology_law', 'values', ['None', 'BuddJacka', 'Cuffey', 'CuffeyTemperate', 'Paterson', 'Arrhenius', 'LliboutryDuval', 'NyeCO2', 'NyeH2O'])
            md = checkfield(md, 'fieldname', 'materials.effectiveconductivity_averaging', 'numel', [1], 'values', [0, 1, 2])

        return md
    # }}}

    def marshall(self, prefix, md, fid):  # {{{
        WriteData(fid, prefix, 'name', 'md.materials.type', 'data', 3, 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'class', 'materials', 'fieldname', 'rho_ice', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'class', 'materials', 'fieldname', 'rho_water', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'class', 'materials', 'fieldname', 'rho_freshwater', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'class', 'materials', 'fieldname', 'mu_water', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'class', 'materials', 'fieldname', 'heatcapacity', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'class', 'materials', 'fieldname', 'latentheat', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'class', 'materials', 'fieldname', 'thermalconductivity', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'class', 'materials', 'fieldname', 'temperateiceconductivity', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'class', 'materials', 'fieldname', 'effectiveconductivity_averaging', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'class', 'materials', 'fieldname', 'meltingpoint', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'class', 'materials', 'fieldname', 'beta', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'class', 'materials', 'fieldname', 'mixed_layer_capacity', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'class', 'materials', 'fieldname', 'thermal_exchange_velocity', 'format', 'Double')
        # NOTE: We first have to check if we have a NumPy array here
        if np.size(self.rheology_B)==1 or (((np.shape(self.rheology_B)[0] == md.mesh.numberofvertices) or (np.shape(self.rheology_B)[0] == md.mesh.numberofvertices + 1)) or ((len(np.shape(self.rheology_B)) == 2) and (np.shape(self.rheology_B)[0] == md.mesh.numberofelements) and (np.shape(self.rheology_B)[1] > 1))):
            mattype = 1
            tsl = md.mesh.numberofvertices
        else:
            mattype = 2
            tsl = md.mesh.numberofelements
        WriteData(fid, prefix, 'object', self, 'class', 'materials', 'fieldname', 'rheology_B', 'format', 'DoubleMat', 'mattype', mattype, 'timeserieslength', tsl + 1, 'yts', md.constants.yts)
        WriteData(fid, prefix, 'object', self, 'class', 'materials', 'fieldname', 'rheology_n', 'format', 'DoubleMat', 'mattype', 2)
        WriteData(fid, prefix, 'data', self.rheology_law, 'name', 'md.materials.rheology_law', 'format', 'String')
        WriteData(fid, prefix, 'object', self, 'class', 'materials', 'fieldname', 'earth_density', 'format', 'Double')
    # }}}
