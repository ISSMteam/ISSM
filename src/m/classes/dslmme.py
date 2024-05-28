from checkfield import checkfield
from fielddisplay import fielddisplay
from project3d import project3d
from WriteData import WriteData


class dslmme(object):
    """DSLMME class definition

    Usage:
        dsl = dslmme() #dynamic sea level class based on a multi-model ensemble of CMIP5 outputs
    """

    def __init__(self, *args):  #{{{
        self.modelid = 0  # Index into the multi-model ensemble, determine which field will be used
        self.global_average_thermosteric_sea_level = [] # Corresponds to zostoga field in CMIP5 archives. Specified as a temporally variable quantity (in m) for each ensemble.
        self.sea_surface_height_above_geoid = [] # Corresponds to zos field in CMIP5 archives. Spatial average is 0. Specified as a spatio-temporally variable quantity (in m) for each ensemble.
        self.sea_water_pressure_at_sea_floor = [] #Corresponds to bpo field in CMIP5 archives. Specified as a spatio-temporally variable quantity (in m equivalent, not in Pa!) for each ensemble.

        if len(args) == 0:
            self.setdefaultparameters()
        else:
            raise Exception('constructor not supported')
    # }}}

    def __repr__(self):  # {{{
        s = '   dsl mme parameters:\n'
        s += '{}\n'.format(fielddisplay(self, 'modelid', 'index into the multi-model ensemble, determines which field will be used.'))
        s += '{}\n'.format(fielddisplay(self, 'global_average_thermosteric_sea_level', 'Corresponds to zostoga field in CMIP5 archives. Specified as a temporally variable quantity (in m) for each ensemble.'))
        s += '{}\n'.format(fielddisplay(self, 'sea_surface_height_above_geoid', 'Corresponds to zos field in CMIP5 archives. Spatial average is 0. Specified as a spatio-temporally variable quantity (in m) for each ensemble.'))
        s += '{}\n'.format(fielddisplay(self, 'sea_water_pressure_at_sea_floor', 'Corresponds to bpo field in CMIP5 archives. Specified as a spatio-temporally variable quantity (in m equivalent, not in Pa!) for each ensemble.'))
        return s
    # }}}

    def setdefaultparameters(self):  #{{{
        return self
    # }}}

    def checkconsistency(self, md, solution, analyses):  # {{{
        # Early return
        if ('SealevelchangeAnalysis' not in analyses) or (solution == 'TransientSolution' and not md.transient.isslc) or (not md.transient.isoceantransport):
            return md

        for i in range(len(self.global_average_thermosteric_sea_level)):
            md = checkfield(md, 'field', self.global_average_thermosteric_sea_level[i], 'NaN', 1, 'Inf', 1)
            md = checkfield(md, 'field', self.sea_surface_height_above_geoid[i], 'NaN', 1, 'Inf', 1, 'timeseries', 1)
            md = checkfield(md, 'field', self.sea_water_pressure_at_sea_floor[i], 'NaN', 1, 'Inf', 1, 'timeseries', 1)
        md = checkfield(md, 'field', self.modelid, 'NaN', 1, 'Inf', 1, '>=', 1, '<=', len(self.global_average_thermosteric_sea_level))

        if self.solidearth.settings.compute_bp_grd:
            md = checkfield(md, 'field', self.sea_water_pressure_at_sea_floor, 'empty', 1)

        return md
    # }}}

    def marshall(self, prefix, md, fid):  #{{{
        WriteData(fid, prefix, 'name', 'md.dsl.model', 'data', 2, 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'modelid', 'format', 'Double')
        WriteData(fid, prefix, 'name', 'md.dsl.nummodels', 'data', len(self.global_average_thermosteric_sea_level), 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'global_average_thermosteric_sea_level', 'format', 'MatArray', 'timeseries', 1, 'timeserieslength', 2)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'sea_water_pressure_at_sea_floor', 'format', 'MatArray', 'timeserieslength', md.mesh.numberofvertices + 1)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'sea_surface_height_above_geoid', 'format', 'MatArray', 'timeserieslength', md.mesh.numberofvertices + 1)
    # }}}

    def extrude(self, md):  #{{{
        for i in range(len(self.global_average_thermosteric_sea_level)):
            self.sea_surface_height_above_geoid[i] = project3d(md, 'vector', self.self.sea_surface_height_above_geoid[i], 'type', 'node', 'layer', 1)
            self.sea_water_pressure_at_sea_floor[i] = project3d(md, 'vector', self.sea_water_pressure_at_sea_floor[i], 'type', 'node', 'layer', 1)

        return self
    # }}}
