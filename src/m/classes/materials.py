import numpy as np
from checkfield import checkfield
from fielddisplay import fielddisplay
from project3d import project3d
from WriteData import WriteData


class materials(object):
    """MATERIALS class definition

    Usage:
        materials = materials()
    """

    def __init__(self, *args):  # {{{
        self.nature = []
        if len(args) == 0:
            self.nature = ['ice']
        else:
            self.nature = args

        for i in range(len(self.nature)):
            if not(self.nature[i] == 'litho' or self.nature[i] == 'ice' or self.nature[i] == 'hydro'):
                raise RuntimeError("materials constructor error message: nature of the material not supported yet! ('ice' or 'litho' or 'hydro')")

        # Start filling in the dynamic fields (not truly dynamic under Python)
        for i in range(len(self.nature)):
            nat = self.nature[i]
            if nat == 'ice':
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
                self.rheology_B = 0
                self.rheology_n = 0
                self.rheology_law = 0
            elif nat == 'litho':
                self.numlayers = 0
                self.radius = 0
                self.viscosity = 0
                self.lame_lambda = 0
                self.lame_mu = 0
                self.burgers_viscosity = 0
                self.burgers_mu = 0
                self.ebm_alpha = 0
                self.ebm_delta = 0
                self.ebm_taul = 0
                self.ebm_tauh = 0
                self.rheologymodel = 0
                self.density = 0
                self.issolid = 0
            elif nat == 'hydro':
                self.rho_ice = 0
                self.rho_water = 0
                self.rho_freshwater = 0
            else:
                raise RuntimeError("materials constructor error message: nature of the material not supported yet! ('ice' or 'litho' or 'hydro')")
        self.earth_density = 0

        self.setdefaultparameters()
    # }}}

    def __repr__(self):  # {{{
        s = '   Materials:\n'
        for i in range(len(self.nature)):
            nat = self.nature[i]
            if nat == 'ice':
                s += '\n      Ice:\n'
                s += '{}\n'.format(fielddisplay(self, 'rho_ice', 'ice density [kg/m^3]'))
                s += '{}\n'.format(fielddisplay(self, 'rho_water', 'ocean water density [kg/m^3]'))
                s += '{}\n'.format(fielddisplay(self, 'rho_freshwater', 'fresh water density [kg/m^3]'))
                s += '{}\n'.format(fielddisplay(self, 'mu_water', 'water viscosity [N s/m^2]'))
                s += '{}\n'.format(fielddisplay(self, 'heatcapacity', 'heat capacity [J/kg/K]'))
                s += '{}\n'.format(fielddisplay(self, 'thermalconductivity', 'ice thermal conductivity [W/m/K]'))
                s += '{}\n'.format(fielddisplay(self, 'temperateiceconductivity', 'temperate ice thermal conductivity [W/m/K]'))
                s += '{}\n'.format(fielddisplay(self, 'meltingpoint', 'melting point of ice at 1atm in K'))
                s += '{}\n'.format(fielddisplay(self, 'latentheat', 'latent heat of fusion [J/m^3]'))
                s += '{}\n'.format(fielddisplay(self, 'beta', 'rate of change of melting point with pressure [K/Pa]'))
                s += '{}\n'.format(fielddisplay(self, 'mixed_layer_capacity', 'mixed layer capacity [W/kg/K]'))
                s += '{}\n'.format(fielddisplay(self, 'thermal_exchange_velocity', 'thermal exchange velocity [m/s]'))
                s += '{}\n'.format(fielddisplay(self, 'rheology_B', 'flow law parameter [Pa s^(1/n)]'))
                s += '{}\n'.format(fielddisplay(self, 'rheology_n', 'Glen\'s flow law exponent'))
                s += '{}\n'.format(fielddisplay(self, 'rheology_law', 'law for the temperature dependance of the rheology: \'None\', \'BuddJacka\', \'Cuffey\', \'CuffeyTemperate\', \'Paterson\', \'Arrhenius\', \'LliboutryDuval\', \'NyeCO2\', or \'NyeH2O\''))
            elif nat == 'litho':
                s += '\n      Litho:\n'
                s += '{}\n'.format(fielddisplay(self, 'numlayers', 'number of layers (default: 2)'))
                s += '{}\n'.format(fielddisplay(self, 'radius', 'array describing the radius for each interface (numlayers + 1) [m]'))
                s += '{}\n'.format(fielddisplay(self, 'viscosity', 'array describing each layer\'s viscosity (numlayers) [Pa.s]'))
                s += '{}\n'.format(fielddisplay(self, 'lame_lambda', 'array describing the lame lambda parameter (numlayers) [Pa]'))
                s += '{}\n'.format(fielddisplay(self, 'lame_mu', 'array describing the shear modulus for each layers (numlayers) [Pa]'))
                s += '{}\n'.format(fielddisplay(self, 'burgers_viscosity', 'array describing each layer\'s transient viscosity, only for Burgers rheologies  (numlayers) [Pa.s]'))
                s += '{}\n'.format(fielddisplay(self, 'burgers_mu', 'array describing each layer\'s transient shear modulus, only for Burgers rheologies  (numlayers) [Pa]'))

                s += '{}\n'.format(fielddisplay(self, 'ebm_alpha', 'array describing each layer\'s exponent parameter controlling the shape of shear modulus curve between taul and tauh, only for EBM rheology (numlayers)'))
                s += '{}\n'.format(fielddisplay(self, 'ebm_delta', 'array describing each layer\'s amplitude of the transient relaxation (ratio between elastic rigity to pre-maxwell relaxation rigity), only for EBM rheology (numlayers)'))
                s += '{}\n'.format(fielddisplay(self, 'ebm_taul', 'array describing each layer\'s starting period for transient relaxation, only for EBM rheology  (numlayers) [s]'))
                s += '{}\n'.format(fielddisplay(self, 'ebm_tauh', 'array describing each layer''s array describing each layer\'s end period for transient relaxation, only for Burgers rheology (numlayers) [s]'))

                s += '{}\n'.format(fielddisplay(self, 'rheologymodel', 'array describing whether we adopt a Maxwell (0), Burgers (1) or EBM (2) rheology (default: 0)'))
                s += '{}\n'.format(fielddisplay(self, 'density', 'array describing each layer\'s density (numlayers) [kg/m^3]'))
                s += '{}\n'.format(fielddisplay(self, 'issolid', 'array describing whether the layer is solid or liquid (default: 1) (numlayers)'))
            elif nat == 'hydro':
                s += '\n      Hydro:\n'
                s += '{}\n'.format(fielddisplay(self, 'rho_ice', 'ice density [kg/m^3]'))
                s += '{}\n'.format(fielddisplay(self, 'rho_water', 'ocean water density [kg/m^3]'))
                s += '{}\n'.format(fielddisplay(self, 'earth_density', 'mantle density [kg/m^3]'))
                s += '{}\n'.format(fielddisplay(self, 'rho_freshwater', 'fresh water density [kg/m^3]'))

            else:
                raise RuntimeError('materials constructor error message: nature of the material not supported yet! (\'ice\' or \'litho\' or \'hydro\')')
        return s
    # }}}

    def setdefaultparameters(self):  # {{{
        for i in range(len(self.nature)):
            nat = self.nature[i]
            if nat == 'ice':
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
                self.latentheat = 3.34 * 1e5

                # Ice thermal conductivity (W/m/K)
                self.thermalconductivity = 2.4

                # Wet ice thermal conductivity (W/m/K)
                self.temperateiceconductivity = 0.24

                # Computation of effective conductivity
                self.effectiveconductivity_averaging = 1

                # The melting point of ice at 1 atmosphere of pressure in K
                self.meltingpoint = 273.15

                # Rate of change of melting point with pressure (K/Pa)
                self.beta = 9.8 * 1e-8

                # Mixed layer (ice-water interface) heat capacity (J/kg/K)
                self.mixed_layer_capacity = 3974.0

                # Thermal exchange velocity (ice-water interface) (m/s)
                self.thermal_exchange_velocity = 1.00 * 1e-4

                # Rheology law: what is the temperature dependence of B with T
                # available: none, paterson and arrhenius
                self.rheology_law = 'Paterson'

                # Rheology fields default
                self.rheology_B = 1 * 1e8
                self.rheology_n = 3
            elif nat == 'litho':
                # We default to a configuration that enables running GIA
                # solutions using giacaron and/or giaivins
                self.numlayers = 2

                # Center of the earth (approximation, must not be 0), then the
                # lab (lithosphere/asthenosphere boundary) then the surface
                # (with 1d3 to avoid numerical singularities)
                self.radius = [1e3, 6278 * 1e3, 6378 * 1e3]

                self.viscosity = [1e21, 1e40] # Mantle and lithosphere viscosity (respectively) [Pa.s]
                self.lame_mu = [1.45 * 1e11, 6.7 * 1e10] # (Pa) # Lithosphere and mantle shear modulus (respectively) [Pa]
                self.lame_lambda = self.lame_mu # (Pa) # Mantle and lithosphere lamba parameter (respectively) [Pa]
                self.burgers_viscosity = [np.nan, np.nan]
                self.burgers_mu = [np.nan, np.nan]

                self.ebm_alpha = [np.nan, np.nan]
                self.ebm_delta = [np.nan, np.nan]
                self.ebm_taul = [np.nan, np.nan]
                self.ebm_tauh = [np.nan, np.nan]
                self.rheologymodel = [0, 0]
                self.density = [5.51 * 1e3, 5.50 * 1e3] # (Pa) # Mantle and lithosphere density [kg/m^3]
                self.issolid = [1, 1] # Is layer solid or liquid?
            elif nat == 'hydro':
                # Ice density (kg/m^3)
                self.rho_ice = 917.0

                # Ocean water density (kg/m^3)
                self.rho_water = 1023.0

                # Fresh water density (kg/m^3)
                self.rho_freshwater = 1000.0
            else:
                raise RuntimeError("materials setdefaultparameters error message: nature of the material not supported yet! ('ice' or 'litho' or 'hydro')")

            # Average density of the Earth (kg/m^3)
            self.earth_density = 5512
    # }}}

    def checkconsistency(self, md, solution, analyses):  # {{{
        for i in range(len(self.nature)):
            nat = self.nature[i]
            if nat == 'ice':
                md = checkfield(md, 'fieldname', 'materials.rho_ice', '>', 0)
                md = checkfield(md, 'fieldname', 'materials.rho_water', '>', 0)
                md = checkfield(md, 'fieldname', 'materials.rho_freshwater', '>', 0)
                md = checkfield(md, 'fieldname', 'materials.mu_water', '>', 0)
                md = checkfield(md, 'fieldname', 'materials.rheology_B', '>', 0, 'timeseries', 1, 'NaN', 1, 'Inf', 1)
                md = checkfield(md, 'fieldname', 'materials.rheology_n', '>', 0, 'size', [md.mesh.numberofelements])
                md = checkfield(md, 'fieldname', 'materials.rheology_law', 'values', ['None', 'BuddJacka', 'Cuffey', 'CuffeyTemperate', 'Paterson', 'Arrhenius', 'LliboutryDuval', 'NyeCO2', 'NyeH2O'])
            elif nat == 'litho':
                if 'LoveAnalysis' not in analyses:
                    return md

                md = checkfield(md, 'fieldname', 'materials.numlayers', 'NaN', 1, 'Inf', 1, '>', 0, 'numel', [1])
                md = checkfield(md, 'fieldname', 'materials.radius', 'NaN', 1, 'Inf', 1, 'size', [md.materials.numlayers + 1, 1], '>', 0)
                md = checkfield(md, 'fieldname', 'materials.lame_mu', 'NaN', 1, 'Inf', 1, 'size', [md.materials.numlayers, 1], '>=', 0)
                md = checkfield(md, 'fieldname', 'materials.lame_lambda', 'NaN', 1, 'Inf', 1, 'size', [md.materials.numlayers, 1], '>=', 0)
                md = checkfield(md, 'fieldname', 'materials.issolid', 'NaN', 1, 'Inf', 1, 'size', [md.materials.numlayers, 1], '>=', 0, '<', 2)
                md = checkfield(md, 'fieldname', 'materials.density', 'NaN', 1, 'Inf', 1, 'size', [md.materials.numlayers, 1], '>', 0)
                md = checkfield(md, 'fieldname', 'materials.viscosity', 'NaN', 1, 'Inf', 1, 'size', [md.materials.numlayers, 1], '>=', 0)
                md = checkfield(md, 'fieldname', 'materials.rheologymodel', 'NaN', 1, 'Inf', 1, 'size', [md.materials.numlayers, 1], '>=', 0, '<=', 2)

                if np.any(self.rheologymodel == 1):
                    md = checkfield(md, 'fieldname', 'materials.burgers_viscosity', 'Inf', 1, 'size', [md.materials.numlayers, 1], '>=', 0)
                    md = checkfield(md, 'fieldname', 'materials.burgers_mu', 'Inf', 1, 'size', [md.materials.numlayers, 1], '>=', 0)

                if np.any(self.rheologymodel == 2):
                    md = checkfield(md, 'fieldname', 'materials.ebm_alpha', 'Inf', 1, 'size', [md.materials.numlayers, 1], '>=', 0)
                    md = checkfield(md, 'fieldname', 'materials.ebm_delta', 'Inf', 1, 'size', [md.materials.numlayers, 1], '>=', 0)
                    md = checkfield(md, 'fieldname', 'materials.ebm_taul', 'Inf', 1, 'size', [md.materials.numlayers, 1], '>=', 0)
                    md = checkfield(md, 'fieldname', 'materials.ebm_tauh', 'Inf', 1, 'size', [md.materials.numlayers, 1], '>=', 0)

                for i in range(md.materials.numlayers):
                    if md.materials.rheologymodel[i] == 1 and (np.isnan(md.materials.burgers_viscosity[i] or np.isnan(md.materials.burgers_mu[i]))):
                        raise RuntimeError('materials checkconsistency error message: Litho burgers_viscosity or burgers_mu has NaN values, inconsistent with rheologymodel choice')

                    if md.materials.rheologymodel[i] == 2 and (np.isnan(md.materials.ebm_alpha[i]) or np.isnan(md.materials.ebm_delta[i]) or np.isnan(md.materials.ebm_taul[i]) or np.isnan(md.materials.ebm_tauh[i])):
                        raise RuntimeError('materials checkconsistency error message: Litho ebm_alpha, ebm_delta, ebm_taul or ebm_tauh has NaN values, inconsistent with rheologymodel choice')
                if md.materials.issolid[0] == 0 or md.materials.lame_mu[0] == 0:
                    raise RuntimeError('First layer must be solid (issolid[0] > 0 AND lame_mu[0] > 0). Add a weak inner core if necessary.')
                ind = np.where(md.materials.issolid == 0)[0]

            elif nat == 'hydro':
                md = checkfield(md, 'fieldname', 'materials.rho_ice', '>', 0)
                md = checkfield(md, 'fieldname', 'materials.rho_water', '>', 0)
                md = checkfield(md, 'fieldname', 'materials.earth_density', '>', 0, 'numel', [1])
                md = checkfield(md, 'fieldname', 'materials.rho_freshwater', '>', 0)
            else:
                raise RuntimeError('materials checkconsistency error message: nature of the material not supported yet! (\'ice\' or \'litho\' or \'hydro\')')

        return md
    # }}}

    def marshall(self, prefix, md, fid):  # {{{
        #1: MatdamageiceEnum 2: MatestarEnum 3: MaticeEnum 4: MatenhancediceEnum 5: MaterialsEnum
        WriteData(fid, prefix, 'name', 'md.materials.nature', 'data', naturetointeger(self.nature), 'format', 'IntMat', 'mattype', 3)
        WriteData(fid, prefix, 'name', 'md.materials.type', 'data', 5, 'format', 'Integer') #DANGER: this can evolve if you have classes
        for i in range(len(self.nature)):
            nat = self.nature[i]
            if nat == 'ice':
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
                WriteData(fid, prefix, 'object', self, 'class', 'materials', 'fieldname', 'rheology_B', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', md.constants.yts)
                WriteData(fid, prefix, 'object', self, 'class', 'materials', 'fieldname', 'rheology_n', 'format', 'DoubleMat', 'mattype', 2)
                WriteData(fid, prefix, 'data', self.rheology_law, 'name', 'md.materials.rheology_law', 'format', 'String')
            elif nat == 'litho':
                WriteData(fid, prefix, 'object', self, 'class', 'materials', 'fieldname', 'numlayers', 'format', 'Integer')
                WriteData(fid, prefix, 'object', self, 'class', 'materials', 'fieldname', 'radius', 'format', 'DoubleMat', 'mattype', 3)
                WriteData(fid, prefix, 'object', self, 'class', 'materials', 'fieldname', 'lame_mu', 'format', 'DoubleMat', 'mattype', 3)
                WriteData(fid, prefix, 'object', self, 'class', 'materials', 'fieldname', 'lame_lambda', 'format', 'DoubleMat', 'mattype', 3)
                WriteData(fid, prefix, 'object', self, 'class', 'materials', 'fieldname', 'issolid', 'format', 'DoubleMat', 'mattype', 3)
                WriteData(fid, prefix, 'object', self, 'class', 'materials', 'fieldname', 'density', 'format', 'DoubleMat', 'mattype', 3)
                WriteData(fid, prefix, 'object', self, 'class', 'materials', 'fieldname', 'viscosity', 'format', 'DoubleMat', 'mattype', 3)
                WriteData(fid, prefix, 'object', self, 'class', 'materials', 'fieldname', 'rheologymodel', 'format', 'DoubleMat', 'mattype', 3)
                WriteData(fid, prefix, 'object', self, 'class', 'materials', 'fieldname', 'burgers_viscosity', 'format', 'DoubleMat', 'mattype', 3)
                WriteData(fid, prefix, 'object', self, 'class', 'materials', 'fieldname', 'burgers_mu', 'format', 'DoubleMat', 'mattype', 3)
                WriteData(fid, prefix, 'object', self, 'class', 'materials', 'fieldname', 'ebm_alpha', 'format', 'DoubleMat', 'mattype', 3)
                WriteData(fid, prefix, 'object', self, 'class', 'materials', 'fieldname', 'ebm_delta', 'format', 'DoubleMat', 'mattype', 3)
                WriteData(fid, prefix, 'object', self, 'class', 'materials', 'fieldname', 'ebm_taul', 'format', 'DoubleMat', 'mattype', 3)
                WriteData(fid, prefix, 'object', self, 'class', 'materials', 'fieldname', 'ebm_tauh', 'format', 'DoubleMat', 'mattype', 3)
                # Compute earth density compatible with our layer density distribution
                earth_density = 0
                for i in range(self.numlayers):
                    earth_density = earth_density + (pow(self.radius[i + 1], 3) - pow(self.radius[i], 3)) * self.density[i]
                earth_density = earth_density / pow(self.radius[self.numlayers], 3)
                self.earth_density = earth_density
            elif nat == 'hydro':
                WriteData(fid, prefix, 'object', self, 'class', 'materials', 'fieldname', 'rho_ice', 'format', 'Double')
                WriteData(fid, prefix, 'object', self, 'class', 'materials', 'fieldname', 'rho_water', 'format', 'Double')
                WriteData(fid, prefix, 'object', self, 'class', 'materials', 'fieldname', 'rho_freshwater', 'format', 'Double')
            else:
                raise RuntimeError('materials constructor error message: nature of the material not supported yet! (\'ice\' or \'litho\' or \'hydro\')')
        WriteData(fid, prefix, 'data', self.earth_density, 'name', 'md.materials.earth_density', 'format', 'Double')
    # }}}

    def extrude(self, md):  # {{{
        for i in range(len(self.nature)):
            nat = self.nature[i]
            if nat == 'ice':
                self.rheology_B = project3d(md, 'vector', self.rheology_B, 'type', 'node')
                self.rheology_n = project3d(md, 'vector', self.rheology_n, 'type', 'element')
            return self
    # }}}
# }}}

def naturetointeger(strnat):  # {{{
    intnat = np.zeros(len(strnat))

    for i in range(len(intnat)):
        str_nat = strnat[i]
        if str_nat == 'damageice':
            intnat[i] = 1
        elif str_nat == 'estar':
            intnat[i] = 2
        elif str_nat == 'ice':
            intnat[i] = 3
        elif str_nat == 'enhancedice':
            intnat[i] = 4
        # elif str_nat == 'materials':
        #     intnat[i] = 5 #this case will never happen, kept to ensure equivalent of codes between IoCodeToMeterialsEnum and IoCodeToNatureEnum
        elif str_nat == 'litho':
            intnat[i] = 6
        elif str_nat == 'hydro':
            intnat[i] = 7
        else:
            raise RuntimeError('materials constructor error message: nature of the material not supported yet! (\'ice\' or \'litho\' or \'hydro\')')

    return intnat
# }}}
