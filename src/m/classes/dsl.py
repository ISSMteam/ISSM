import numpy as np

from checkfield import checkfield
from fielddisplay import fielddisplay
from project3d import project3d
from WriteData import WriteData


class dsl(object):
    """DSL - class definition

    Usage:
        dsl = dsl() # dynamic sea level class, based on CMIP5 outputs
    """

    def __init__(self, *args):  #{{{
        self.global_average_thermosteric_sea_level = np.nan  # Corresponds to zostoga field in CMIP5 archives. Specified as a temporally variable quantity (in m).
        self.sea_surface_height_above_geoid = np.nan  # Corresponds to zos field in CMIP5 archives. Spatial average is 0. Specified as a spatio-temporally variable quantity (in m).
        self.sea_water_pressure_at_sea_floor = np.nan  #  Corresponds to bpo field in CMIP5 archives. Specified as a spatio-temporally variable quantity (in m equivalent, not in Pa!).

        if len(args) == 0:
            self.setdefaultparameters()
        else:
            raise Exception('constructor not supported')
    # }}}

    def __repr__(self):  #{{{
        s = '   dsl parameters:\n'
        s += '{}\n'.format(fielddisplay(self, 'global_average_thermosteric_sea_level', 'Corresponds to zostoga field in CMIP5 archives. Specified as a temporally variable quantity (in m).'))
        s += '{}\n'.format(fielddisplay(self, 'sea_surface_height_above_geoid', 'Corresponds to zos field in CMIP5 archives. Spatial average is 0. Specified as a spatio-temporally variable quantity (in m).'))
        s += '{}\n'.format(fielddisplay(self, 'sea_water_pressure_at_sea_floor', 'Corresponds to bpo field in CMIP5 archives. Specified as a spatio-temporally variable quantity (in m equivalent, not in Pa!).'))
        return s
    # }}}

    def setdefaultparameters(self):  #{{{
        self.global_average_thermosteric_sea_level = np.nan
        self.sea_surface_height_above_geoid = np.nan
        self.sea_water_pressure_at_sea_floor = np.nan
    # }}}

    def checkconsistency(self, md, solution, analyses):  #{{{
        # Early return
        if ('SealevelchangeAnalysis' not in analyses) or (solution == 'TransientSolution' and not md.transient.isslc) or (not md.transient.isoceantransport):
            return md
        md = checkfield(md, 'fieldname', 'dsl.global_average_thermosteric_sea_level', 'NaN', 1, 'Inf', 1)
        md = checkfield(md, 'fieldname', 'dsl.sea_surface_height_above_geoid', 'NaN', 1, 'Inf', 1, 'timeseries', 1)
        md = checkfield(md, 'fieldname', 'dsl.sea_water_pressure_at_sea_floor', 'NaN', 1, 'Inf', 1, 'timeseries', 1)

        if md.solidearth.settings.compute_bp_grd:
            md = checkfield(md, 'fieldname', 'dsl.sea_water_pressure_at_sea_floor', 'empty', 1)

        return md
    # }}}

    def marshall(self, prefix, md, fid):   #{{{
        yts = md.constants.yts
        WriteData(fid, prefix, 'name', 'md.dsl.model', 'data', 1, 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'global_average_thermosteric_sea_level', 'format', 'DoubleMat', 'mattype', 2, 'timeseries', 1, 'timeserieslength', 2, 'yts', yts) # mattype 2, because we are sending a GMSL value identical everywhere on each element.
        WriteData(fid, prefix, 'object', self, 'fieldname', 'sea_surface_height_above_geoid', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', yts)  # mattype 1 because we specify DSL at vertex locations.
        WriteData(fid, prefix, 'object', self, 'fieldname', 'sea_water_pressure_at_sea_floor', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', yts)  # mattype 1 because we specify bottom pressure at vertex locations.
    # }}}

    def extrude(self, md):  #{{{
        self.sea_surface_height_above_geoid = project3d(md, 'vector', self.sea_surface_height_above_geoid, 'type', 'node', 'layer', 1)
        self.sea_water_pressure_at_sea_floor = project3d(md, 'vector', self.sea_water_pressure_at_sea_floor, 'type', 'node', 'layer', 1)
        return self
    # }}}

    def initialize(self, md):  #{{{
        if np.all(np.isnan(self.global_average_thermosteric_sea_level)):
            self.global_average_thermosteric_sea_level = np.array([0, 0]).reshape(-1, 1)
            print('      no dsl.global_average_thermosteric_sea_level specified: transient values set to zero')

        if np.all(np.isnan(self.sea_surface_height_above_geoid)):
            self.sea_surface_height_above_geoid = np.append(np.zeros((md.mesh.numberofvertices, 1)), 0).reshape(-1, 1)
            print('      no dsl.sea_surface_height_above_geoid specified: transient values set to zero')

        if np.all(np.isnan(self.sea_water_pressure_at_sea_floor)):
            self.sea_water_pressure_at_sea_floor = np.append(np.zeros((md.mesh.numberofvertices, 1)), 0).reshape(-1, 1)
            print('      no dsl.sea_water_pressure_at_sea_floor specified: transient values set to zero')
    # }}}
