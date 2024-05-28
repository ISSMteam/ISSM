import numpy as np

from checkfield import checkfield
from fielddisplay import fielddisplay
from project3d import project3d
from structtoobj import structtoobj
from WriteData import WriteData


class spatiallinearbasalforcings(object):
    """SPATIAL LINEAR BASAL FORCINGS class definition

    Usage:
        spatiallinearbasalforcings = spatiallinearbasalforcings()
    """

    def __init__(self, *args):  # {{{
        nargs = len(args)
        if nargs == 0:
            self.groundedice_melting_rate   = np.nan
            self.deepwater_melting_rate     = np.nan
            self.deepwater_elevation        = np.nan
            self.upperwater_melting_rate    = np.nan
            self.upperwater_elevation       = np.nan
            self.geothermalflux             = np.nan
            self.perturbation_melting_rate  = np.nan

            self.setdefaultparameters()
        elif nargs == 1:
            lb = args[0]
            if lb.__module__ == 'linearbasalforcings':
                nvertices = len(lb.groundedice_melting_rate)
                self.groundedice_melting_rate   = lb.groundedice_melting_rate
                self.geothermalflux             = lb.geothermalflux
                self.deepwater_elevation        = lb.deepwater_elevation * np.ones((nvertices, ))
                self.deepwater_melting_rate     = lb.deepwater_melting_rate * np.ones((nvertices, ))
                self.upperwater_melting_rate    = lb.upperwater_melting_rate * np.ones((nvertices, ))
                self.upperwater_elevation       = lb.upperwater_elevation * np.ones((nvertices, ))
                self.perturbation_melting_rate  = lb.perturbation_melting_rate * np.ones((nvertices, ))
            else:
                # TODO: This has not been tested
                self = structtoobj(self, lb)
        else:
            raise Exception('constructor not supported')
    # }}}

    def __repr__(self):  # {{{
        s = '   spatial linear basal forcings parameters:\n'
        s += '{}\n'.format(fielddisplay(self, 'groundedice_melting_rate', 'basal melting rate (positive if melting) [m/yr]'))
        s += '{}\n'.format(fielddisplay(self, 'deepwater_melting_rate', 'basal melting rate (positive if melting applied for floating ice whith base < deepwater_elevation) [m/yr]'))
        s += '{}\n'.format(fielddisplay(self, 'deepwater_elevation', 'elevation of ocean deepwater [m]'))
        s += '{}\n'.format(fielddisplay(self, 'upperwater_melting_rate', 'basal melting rate (positive if melting applied for floating ice whith base >= upperwater_elevation) [m/yr]'))
        s += '{}\n'.format(fielddisplay(self, 'upperwater_elevation', 'elevation of ocean upperwater [m]'))
        s += '{}\n'.format(fielddisplay(self, 'perturbation_melting_rate', 'basal melting rate perturbation added to computed melting rate (positive if melting) [m/yr]'))
        s += '{}\n'.format(fielddisplay(self, 'geothermalflux', 'geothermal heat flux [W/m^2]'))
        return s
    # }}}

    def extrude(self, md):  # {{{
        self.groundedice_melting_rate = project3d(md, 'vector', self.groundedice_melting_rate, 'type', 'node', 'layer', 1) 
        self.deepwater_melting_rate = project3d(md, 'vector', self.deepwater_melting_rate, 'type', 'node', 'layer', 1) 
        self.deepwater_elevation = project3d(md, 'vector', self.deepwater_elevation, 'type', 'node', 'layer', 1)
        self.upperwater_melting_rate = project3d(md, 'vector', self.upperwater_melting_rate, 'type', 'node', 'layer', 1) 
        self.upperwater_elevation = project3d(md, 'vector', self.upperwater_elevation, 'type', 'node', 'layer', 1) 
        self.geothermalflux = project3d(md, 'vector', self.geothermalflux, 'type', 'node', 'layer', 1) # Bedrock only gets geothermal flux
        self.perturbation_melting_rate = project3d(md, 'vector', self.upperwater_melting_rate, 'type', 'node', 'layer', 1) 
        return self
    # }}}

    def initialize(self, md):  # {{{
        if np.all(np.isnan(self.groundedice_melting_rate)):
            self.groundedice_melting_rate = np.zeros((md.mesh.numberofvertices))
            print('      no basalforcings.groundedice_melting_rate specified: values set as zero')
        return self
    # }}}

    def setdefaultparameters(self):  # {{{
        return self
    # }}}

    def checkconsistency(self, md, solution, analyses):  # {{{
        if not np.all(np.isnan(self.perturbation_melting_rate)):
            md = checkfield(md, 'fieldname', 'basalforcings.perturbation_melting_rate', 'NaN', 1, 'Inf', 1, 'timeseries', 1)
        if 'MasstransportAnalysis' in analyses and not solution == 'TransientSolution' and not md.transient.ismasstransport:
            md = checkfield(md, 'fieldname', 'basalforcings.groundedice_melting_rate', 'NaN', 1, 'Inf', 1, 'timeseries', 1)
            md = checkfield(md, 'fieldname', 'basalforcings.deepwater_melting_rate', 'NaN', 1, 'Inf', 1, 'timeseries', 1, '>=', 0)
            md = checkfield(md, 'fieldname', 'basalforcings.deepwater_elevation', 'NaN', 1, 'Inf', 1, 'timeseries', 1)
            md = checkfield(md, 'fieldname', 'basalforcings.upperwater_melting_rate', 'NaN', 1, 'Inf', 1, 'timeseries', 1, '>=', 0)
            md = checkfield(md, 'fieldname', 'basalforcings.upperwater_elevation', 'NaN', 1, 'Inf', 1, 'timeseries', 1, '<', 0)
        if 'BalancethicknessAnalysis' in analyses:
            raise Exception('not implemented yet!')
            md = checkfield(md, 'fieldname', 'basalforcings.groundedice_melting_rate', 'NaN', 1, 'Inf', 1, 'size', [md.mesh.numberofvertices])
            md = checkfield(md, 'fieldname', 'basalforcings.deepwater_melting_rate', '>=', 0, 'numel', 1)
            md = checkfield(md, 'fieldname', 'basalforcings.deepwater_elevation', '<', 'basalforcings.upperwater_elevation', 'numel', 1)
            md = checkfield(md, 'fieldname', 'basalforcings.upperwater_melting_rate', '>=', 0, 'numel', 1)
            md = checkfield(md, 'fieldname', 'basalforcings.upperwater_elevation', '<=', 0, 'numel', 1)
        if 'ThermalAnalysis' in analyses and not solution == 'TransientSolution' and not md.transient.isthermal:
            raise Exception('not implemented yet!')
            md = checkfield(md, 'fieldname', 'basalforcings.groundedice_melting_rate', 'NaN', 1, 'Inf', 1, 'timeseries', 1)
            md = checkfield(md, 'fieldname', 'basalforcings.deepwater_melting_rate', '>=', 0, 'numel', 1)
            md = checkfield(md, 'fieldname', 'basalforcings.deepwater_elevation', '<', 'basalforcings.upperwater_elevation', 'numel', 1)
            md = checkfield(md, 'fieldname', 'basalforcings.upperwater_elevation', '<', 0, 'numel', 1)
            md = checkfield(md, 'fieldname', 'basalforcings.deepwater_melting_rate', '>=', 0, 'numel', 1)
            md = checkfield(md, 'fieldname', 'basalforcings.geothermalflux', 'NaN', 1, 'Inf', 1, 'timeseries', 1, '>=', 0)
        return md
    # }}}

    def marshall(self, prefix, md, fid):  # {{{
        yts = md.constants.yts

        WriteData(fid, prefix, 'name', 'md.basalforcings.model', 'data', 6, 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'groundedice_melting_rate', 'format', 'DoubleMat', 'name', 'md.basalforcings.groundedice_melting_rate', 'mattype', 1, 'scale', 1. / yts, 'timeserieslength', md.mesh.numberofvertices + 1,'yts', md.constants.yts)
        WriteData(fid, prefix,'object', self, 'fieldname', 'geothermalflux', 'name', 'md.basalforcings.geothermalflux', 'format', 'DoubleMat', 'mattype', 1,'timeserieslength', md.mesh.numberofvertices + 1, 'yts', md.constants.yts)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'deepwater_melting_rate', 'format', 'DoubleMat', 'name', 'md.basalforcings.deepwater_melting_rate', 'scale', 1. / yts, 'mattype', 1)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'deepwater_elevation', 'format', 'DoubleMat', 'name', 'md.basalforcings.deepwater_elevation', 'mattype', 1)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'upperwater_melting_rate', 'format', 'DoubleMat', 'name', 'md.basalforcings.upperwater_melting_rate', 'scale', 1. / yts, 'mattype', 1)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'upperwater_elevation', 'format', 'DoubleMat', 'name', 'md.basalforcings.upperwater_elevation', 'mattype', 1)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'perturbation_melting_rate', 'format', 'DoubleMat', 'name', 'md.basalforcings.perturbation_melting_rate', 'scale', 1. / yts, 'mattype', 1)
    # }}}
