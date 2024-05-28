import numpy as np

from checkfield import checkfield
from fielddisplay import fielddisplay
from project3d import project3d
from WriteData import WriteData


class initialization(object):
    """initialization class definition

    Usage:
        initialization = initialization()
    """

    def __init__(self):  #{{{
        self.vx                  = np.nan
        self.vy                  = np.nan
        self.vz                  = np.nan
        self.vel                 = np.nan
        self.pressure            = np.nan
        self.temperature         = np.nan
        self.enthalpy            = np.nan
        self.waterfraction       = np.nan
        self.sediment_head       = np.nan
        self.epl_head            = np.nan
        self.epl_thickness       = np.nan
        self.watercolumn         = np.nan
        self.hydraulic_potential = np.nan
        self.channelarea         = np.nan
        self.sealevel            = np.nan
        self.bottompressure      = np.nan
        self.dsl                 = np.nan
        self.str                 = np.nan
        self.sample              = np.nan
        self.debris              = np.nan
        self.age                 = np.nan

        self.setdefaultparameters()
    # }}}

    def __repr__(self):  #{{{
        s = '   initial field values:\n'
        s += '{}\n'.format(fielddisplay(self, 'vx', 'x component of velocity [m/yr]'))
        s += '{}\n'.format(fielddisplay(self, 'vy', 'y component of velocity [m/yr]'))
        s += '{}\n'.format(fielddisplay(self, 'vz', 'z component of velocity [m/yr]'))
        s += '{}\n'.format(fielddisplay(self, 'vel', 'velocity norm [m/yr]'))
        s += '{}\n'.format(fielddisplay(self, 'pressure', 'pressure [Pa]'))
        s += '{}\n'.format(fielddisplay(self, 'temperature', 'temperature [K]'))
        s += '{}\n'.format(fielddisplay(self, 'enthalpy', 'enthalpy [J]'))
        s += '{}\n'.format(fielddisplay(self, 'waterfraction', 'fraction of water in the ice'))
        s += '{}\n'.format(fielddisplay(self, 'watercolumn', 'thickness of subglacial water [m]'))
        s += '{}\n'.format(fielddisplay(self, 'sediment_head', 'sediment water head of subglacial system [m]'))
        s += '{}\n'.format(fielddisplay(self, 'epl_head', 'epl water head of subglacial system [m]'))
        s += '{}\n'.format(fielddisplay(self, 'epl_thickness', 'thickness of the epl [m]'))
        s += '{}\n'.format(fielddisplay(self, 'hydraulic_potential', 'Hydraulic potential (for GlaDS) [Pa]'))
        s += '{}\n'.format(fielddisplay(self, 'channelarea', 'subglaciale water channel area (for GlaDS) [m2]'))
        s += '{}\n'.format(fielddisplay(self, 'sample', 'Realization of a Gaussian random field'))
        s += '{}\n'.format(fielddisplay(self, 'debris', 'Surface debris layer [m]'))
        s += '{}\n'.format(fielddisplay(self, 'age', 'Initial age [yr]'))
        return s
    # }}}

    def setdefaultparameters(self):  #{{{
        return
    # }}}

    def checkconsistency(self, md, solution, analyses):  #{{{
        if 'StressbalanceAnalysis' in analyses and not solution == 'TransientSolution' and not md.transient.isstressbalance:
            if not np.any(np.logical_or(np.isnan(md.initialization.vx), np.isnan(md.initialization.vy))):
                if np.size(md.initialization.vx) > 1 or np.size(md.initialization.vy) > 1:
                    md = checkfield(md, 'fieldname', 'initialization.vx', 'NaN', 1, 'Inf', 1, 'size', [md.mesh.numberofvertices])
                    md = checkfield(md, 'fieldname', 'initialization.vy', 'NaN', 1, 'Inf', 1, 'size', [md.mesh.numberofvertices])
        if 'MasstransportAnalysis' in analyses and not solution == 'TransientSolution' and not md.transient.ismasstransport:
            md = checkfield(md, 'fieldname', 'initialization.vx', 'NaN', 1, 'Inf', 1, 'size', [md.mesh.numberofvertices])
            md = checkfield(md, 'fieldname', 'initialization.vy', 'NaN', 1, 'Inf', 1, 'size', [md.mesh.numberofvertices])
        if 'OceantransportAnalysis' in analyses:
            if solution == 'TransientSolution' and md.transient.isslc and md.transient.isoceantransport:
                md = checkfield(md, 'fieldname', 'initialization.bottompressure', 'NaN', 1, 'Inf', 1, 'size', [md.mesh.numberofvertices])
                md = checkfield(md, 'fieldname', 'initialization.dsl', 'NaN', 1, 'Inf', 1, 'size', [md.mesh.numberofvertices])
                md = checkfield(md, 'fieldname', 'initialization.str', 'NaN', 1, 'Inf', 1, 'size', [1])
        if 'BalancethicknessAnalysis' in analyses and solution == 'BalancethicknessSolution':
            md = checkfield(md, 'fieldname', 'initialization.vx', 'NaN', 1, 'Inf', 1, 'size', [md.mesh.numberofvertices])
            md = checkfield(md, 'fieldname', 'initialization.vy', 'NaN', 1, 'Inf', 1, 'size', [md.mesh.numberofvertices])
            # Triangle with zero velocity
            if np.any(np.logical_and(np.sum(np.abs(md.initialization.vx[md.mesh.elements - 1]), axis=1).reshape(-1, 1) == 0, np.sum(np.abs(md.initialization.vy[md.mesh.elements - 1]), axis=1).reshape(-1, 1) == 0, np.min(md.mask.ice_levelset[md.mesh.elements - 1], axis=1).reshape(-1, 1) < 0)):
                md.checkmessage('at least one triangle has all its vertices with a zero velocity')
        if 'ThermalAnalysis' in analyses and not solution == 'TransientSolution' and not md.transient.isthermal:
            md = checkfield(md, 'fieldname', 'initialization.vx', 'NaN', 1, 'Inf', 1, 'size', [md.mesh.numberofvertices])
            md = checkfield(md, 'fieldname', 'initialization.vy', 'NaN', 1, 'Inf', 1, 'size', [md.mesh.numberofvertices])
            md = checkfield(md, 'fieldname', 'initialization.temperature', 'NaN', 1, 'Inf', 1, 'size', 'universal')
            if md.mesh.dimension() == 3:
                md = checkfield(md, 'fieldname', 'initialization.vz', 'NaN', 1, 'Inf', 1, 'size', [md.mesh.numberofvertices])
            md = checkfield(md, 'fieldname', 'initialization.pressure', 'NaN', 1, 'Inf', 1, 'size', [md.mesh.numberofvertices])
        if 'EnthalpyAnalysis' in analyses and md.thermal.isenthalpy:
            md = checkfield(md, 'fieldname', 'initialization.waterfraction', '>=', 0, 'size', [md.mesh.numberofvertices])
            md = checkfield(md, 'fieldname', 'initialization.watercolumn', '>=', 0, 'size', [md.mesh.numberofvertices])
            pos = np.nonzero(md.initialization.waterfraction > 0.)[0]
            if(pos.size):
                md = checkfield(md, 'fieldname', 'delta Tpmp', 'field', np.absolute(md.initialization.temperature[pos] - (md.materials.meltingpoint - md.materials.beta * md.initialization.pressure[pos])), '<', 1e-11, 'message', 'set temperature to pressure melting point at locations with waterfraction > 0')
        if 'HydrologyShreveAnalysis' in analyses:
            if type(md.hydrology).__name__ == 'hydrologyshreve':
                if (solution == 'TransientSolution' and md.transient.ishydrology) or solution == 'HydrologySolution':
                    md = checkfield(md, 'fieldname', 'initialization.watercolumn', 'NaN', 1, 'Inf', 1, 'size', [md.mesh.numberofvertices])
        if 'HydrologyTwsAnalysis' in analyses:
            if type(md.hydrology).__name__ == 'hydrologytws':
                md = checkfield(md, 'fieldname', 'initialization.watercolumn', 'NaN', 1, 'Inf', 1, 'size', [md.mesh.numberofvertices])
        if 'SealevelchangeAnalysis' in analyses:
            if solution == 'TransientSolution' and md.transient.isslc:
                md = checkfield(md, 'fieldname', 'initialization.sealevel', 'NaN', 1, 'Inf', 1, 'size', [md.mesh.numberofvertices])
        if 'HydrologyGlaDSAnalysis' in analyses:
            if type(md.hydrology).__name__ == 'hydrologyglads':
                md = checkfield(md, 'fieldname', 'initialization.watercolumn', 'NaN', 1, 'Inf', 1, 'size', [md.mesh.numberofvertices])
                md = checkfield(md, 'fieldname', 'initialization.hydraulic_potential', 'NaN', 1, 'Inf', 1, 'size', [md.mesh.numberofvertices])
                md = checkfield(md, 'fieldname', 'initialization.channelarea', 'NaN', 1, 'Inf', 1, '>=', 0, 'size', [md.mesh.numberofelements])
        if 'HydrologyDCInefficientAnalysis' in analyses:
            if type(md.hydrology).__name__ == 'hydrologydc':
                md = checkfield(md, 'fieldname', 'initialization.sediment_head', 'NaN', 1, 'Inf', 1, 'size', [md.mesh.numberofvertices])
        if 'HydrologyDCEfficientAnalysis' in analyses:
            if type(md.hydrology).__name__ == 'hydrologydc':
                if md.hydrology.isefficientlayer:
                    md = checkfield(md, 'fieldname', 'initialization.epl_head', 'NaN', 1, 'Inf', 1, 'size', [md.mesh.numberofvertices])
                    md = checkfield(md, 'fieldname', 'initialization.epl_thickness', 'NaN', 1, 'Inf', 1, 'size', [md.mesh.numberofvertices])
        if 'SamplingAnalysis' in analyses and not solution == 'TransientSolution' and not md.transient.issampling:
            if np.any(np.isnan(md.initialization.sample)):
                md = checkfield(md, 'fieldname', 'initialization.sample', 'NaN', 1,'Inf', 1, 'size', [md.mesh.numberofvertices])
        if 'DebrisAnalysis' in analyses:
            if not np.isnan(md.initialization.debris):
                if (solution == 'TransientSolution' and md.transient.ishydrology) or solution == 'HydrologySolution':
                    md = checkfield(md, 'fieldname', 'initialization.debris', 'NaN', 1,'Inf', 1, 'size', [md.mesh.numberofvertices])
        if 'AgeAnalysis' in analyses:
            if not np.isnan(md.initialization.age):
                if (solution == 'TransientSolution' and md.transient.ishydrology) or solution == 'HydrologySolution':
                    md = checkfield(md, 'fieldname', 'initialization.age', 'NaN', 1,'Inf', 1, 'size', [md.mesh.numberofvertices])
        return md
    # }}}

    def marshall(self, prefix, md, fid):  # {{{
        yts = md.constants.yts

        WriteData(fid, prefix, 'object', self, 'fieldname', 'vx', 'format', 'DoubleMat', 'mattype', 1, 'scale', 1 / yts, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', yts)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'vy', 'format', 'DoubleMat', 'mattype', 1, 'scale', 1 / yts, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', yts)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'vz', 'format', 'DoubleMat', 'mattype', 1, 'scale', 1 / yts)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'pressure', 'format', 'DoubleMat', 'mattype', 1)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'sealevel', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', yts)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'bottompressure', 'format', 'DoubleMat', 'mattype', 1)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'str', 'format', 'DoubleMat', 'mattype', 1)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'dsl', 'format', 'DoubleMat', 'mattype', 1)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'temperature', 'format', 'DoubleMat', 'mattype', 1)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'waterfraction', 'format', 'DoubleMat', 'mattype', 1)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'sediment_head', 'format', 'DoubleMat', 'mattype', 1)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'epl_head', 'format', 'DoubleMat', 'mattype', 1)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'epl_thickness', 'format', 'DoubleMat', 'mattype', 1)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'watercolumn', 'format', 'DoubleMat', 'mattype', 1)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'channelarea', 'format', 'DoubleMat', 'mattype', 1)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'hydraulic_potential', 'format', 'DoubleMat', 'mattype', 1)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'sample', 'format', 'DoubleMat', 'mattype', 1)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'debris', 'format', 'DoubleMat', 'mattype', 1)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'age', 'format', 'DoubleMat', 'mattype', 1, 'scale', yts)

        if md.thermal.isenthalpy:
            if (np.size(self.enthalpy) <= 1):
                # Reconstruct enthalpy
                tpmp = md.materials.meltingpoint - md.materials.beta * md.initialization.pressure
                pos = np.where(md.initialization.waterfraction > 0)[0]
                self.enthalpy = md.materials.heatcapacity * (md.initialization.temperature - md.constants.referencetemperature)
                self.enthalpy[pos] = md.materials.heatcapacity * (tpmp[pos].reshape(-1,) - md.constants.referencetemperature) + md.materials.latentheat * md.initialization.waterfraction[pos].reshape(-1,)

            WriteData(fid, prefix, 'data', self.enthalpy, 'format', 'DoubleMat', 'mattype', 1, 'name', 'md.initialization.enthalpy')
    # }}}

    def extrude(self, md):  # {{{
        self.vx = project3d(md, 'vector', self.vx, 'type', 'node')
        self.vy = project3d(md, 'vector', self.vy, 'type', 'node')
        self.vz = project3d(md, 'vector', self.vz, 'type', 'node')
        self.vel = project3d(md, 'vector', self.vel, 'type', 'node')
        self.temperature = project3d(md, 'vector', self.temperature, 'type', 'node')
        self.enthalpy = project3d(md, 'vector', self.enthalpy, 'type', 'node')
        self.waterfraction = project3d(md, 'vector', self.waterfraction, 'type', 'node')
        self.watercolumn = project3d(md, 'vector', self.watercolumn, 'type', 'node')
        self.sediment_head = project3d(md, 'vector', self.sediment_head, 'type', 'node', 'layer', 1)
        self.epl_head = project3d(md, 'vector', self.epl_head, 'type', 'node', 'layer', 1)
        self.epl_thickness = project3d(md, 'vector', self.epl_thickness, 'type', 'node', 'layer', 1)
        self.sealevel = project3d(md, 'vector', self.sealevel, 'type', 'node', 'layer', 1)
        self.bottompressure = project3d(md, 'vector', self.bottompressure, 'type', 'node', 'layer', 1)
        self.dsl = project3d(md, 'vector', self.dsl, 'type', 'node', 'layer', 1)
        self.str = project3d(md, 'vector', self.str, 'type', 'node', 'layer', 1)
        self.debris = project3d(md, 'vector', self.debris, 'type', 'node', 'layer', 1)
        self.age = project3d(md, 'vector', self.age, 'type', 'node', 'layer', 1)

        # Lithostatic pressure by default
        if np.ndim(md.geometry.surface) == 2:
            print('Reshaping md.geometry.surface for your convenience but you should fix it in your model set up')
            self.pressure = md.constants.g * md.materials.rho_ice * (md.geometry.surface.reshape(-1, 1) - md.mesh.z)
        else:
            self.pressure = md.constants.g * md.materials.rho_ice * (md.geometry.surface - md.mesh.z)

        return self
    # }}}
