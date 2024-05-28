import numpy as np

from checkfield import checkfield
from fielddisplay import fielddisplay
from project3d import project3d
from WriteData import WriteData


class linearbasalforcings(object):
    """LINEAR BASAL FORCINGS class definition

    Usage:
        basalforcings = linearbasalforcings()
    """

    def __init__(self, *args):  # {{{
        nargs = len(args)
        if nargs == 0:
            print('empty init')
            self.deepwater_melting_rate = 0.
            self.deepwater_elevation = 0.
            self.upperwater_melting_rate = 0.
            self.upperwater_elevation = 0.
            self.groundedice_melting_rate = float('NaN')
            self.perturbation_melting_rate = float('NaN')
            self.geothermalflux = float('NaN')

            # set defaults
            self.setdefaultparameters()
        elif nargs == 1 and args[0].__module__ == 'basalforcings':
            print('converting basalforings to linearbasalforcings')
            inv = args[0]
            self.groundedice_melting_rate = inv.groundedice_melting_rate
            self.perturbation_melting_rate = float('NaN')
            self.geothermalflux = inv.geothermalflux
            self.deepwater_melting_rate = 0.
            self.deepwater_elevation = 0.
            self.upperwater_melting_rate = 0.
            self.upperwater_elevation = 0.

            # Set defaults
            self.setdefaultparameters()
        else:
            raise Exception('constructor not supported')
    # }}}

    def __repr__(self):  # {{{
        s = '   linear basal forcings parameters:\n'
        s += '{}\n'.format(fielddisplay(self, "deepwater_melting_rate", "basal melting rate (positive if melting applied for floating ice whith base < deepwater_elevation) [m/yr]"))
        s += '{}\n'.format(fielddisplay(self, "deepwater_elevation", "elevation of ocean deepwater [m]"))
        s += '{}\n'.format(fielddisplay(self, "upperwater_melting_rate", "upper melting rate (positive if melting applied for floating ice whith base >= upperwater_elevation) [m/yr]"))
        s += '{}\n'.format(fielddisplay(self, "upperwater_elevation", "elevation of ocean upper water [m]"))
        s += '{}\n'.format(fielddisplay(self, "groundedice_melting_rate", "basal melting rate (positive if melting) [m/yr]"))
        s += '{}\n'.format(fielddisplay(self, "perturbation_melting_rate", "perturbation applied to computed melting rate (positive if melting) [m/yr]"))
        s += '{}\n'.format(fielddisplay(self, "geothermalflux", "geothermal heat flux [W/m^2]"))
        return s
    # }}}

    def extrude(self, md):  # {{{
        self.perturbation_melting_rate = project3d(md, 'vector', self.perturbation_melting_rate, 'type', 'node', 'layer', 1)
        self.groundedice_melting_rate = project3d(md, 'vector', self.groundedice_melting_rate, 'type', 'node', 'layer', 1)
        self.geothermalflux = project3d(md, 'vector', self.geothermalflux, 'type', 'node', 'layer', 1) # Bedrock only gets geothermal flux
        return self
    # }}}

    def initialize(self, md):  # {{{
        if np.all(np.isnan(self.groundedice_melting_rate)):
            self.groundedice_melting_rate = np.zeros((md.mesh.numberofvertices))
            print("      no basalforcings.groundedice_melting_rate specified: values set as zero")
        return self
    # }}}

    def setdefaultparameters(self):  # {{{
        self.deepwater_melting_rate = 50.0
        self.deepwater_elevation = -800.0
        self.upperwater_melting_rate = 0.0
        self.upperwater_elevation = -400.0
        return self
    # }}}

    def checkconsistency(self, md, solution, analyses):  # {{{
        if not np.all(np.isnan(self.perturbation_melting_rate)):
            md = checkfield(md, 'fieldname', 'basalforcings.perturbation_melting_rate', 'NaN', 1, 'Inf', 1, 'timeseries', 1)
        if 'MasstransportAnalysis' in analyses and not solution == 'TransientSolution' and not md.transient.ismasstransport:
            md = checkfield(md, 'fieldname', 'basalforcings.groundedice_melting_rate', 'NaN', 1, 'Inf', 1, 'timeseries', 1)
            md = checkfield(md, 'fieldname', 'basalforcings.deepwater_melting_rate', '>=', 0, 'singletimeseries', 1)
            md = checkfield(md, 'fieldname', 'basalforcings.upperwater_melting_rate', '>=', 0, 'singletimeseries', 1)
            md = checkfield(md, 'fieldname', 'basalforcings.deepwater_elevation', '<', 'basalforcings.upperwater_elevation', 'singletimeseries', 1)
            md = checkfield(md, 'fieldname', 'basalforcings.upperwater_elevation', '<=', 0, 'singletimeseries', 1)
        if 'BalancethicknessAnalysis' in analyses:
            md = checkfield(md, 'fieldname', 'basalforcings.groundedice_melting_rate', 'NaN', 1, 'Inf', 1, 'size', [md.mesh.numberofvertices])
            md = checkfield(md, 'fieldname', 'basalforcings.deepwater_melting_rate', '>=', 0, 'singletimeseries', 1)
            md = checkfield(md, 'fieldname', 'basalforcings.upperwater_melting_rate', '>=', 0, 'singletimeseries', 1)
            md = checkfield(md, 'fieldname', 'basalforcings.deepwater_elevation', '<', 'basalforcings.upperwater_elevation', 'singletimeseries', 1)
            md = checkfield(md, 'fieldname', 'basalforcings.upperwater_elevation', '<=', 0, 'singletimeseries', 1)
        if 'ThermalAnalysis' in analyses and not solution == 'TransientSolution' and not md.transient.isthermal:
            md = checkfield(md, 'fieldname', 'basalforcings.groundedice_melting_rate', 'NaN', 1, 'Inf', 1, 'timeseries', 1)
            md = checkfield(md, 'fieldname', 'basalforcings.deepwater_melting_rate', '>=', 0, 'singletimeseries', 1)
            md = checkfield(md, 'fieldname', 'basalforcings.upperwater_melting_rate', '>=', 0, 'singletimeseries', 1)
            md = checkfield(md, 'fieldname', 'basalforcings.deepwater_elevation', '<', 'basalforcings.upperwater_elevation', 'singletimeseries', 1)
            md = checkfield(md, 'fieldname', 'basalforcings.upperwater_elevation', '<=', 0, 'singletimeseries', 1)
            md = checkfield(md, 'fieldname', 'basalforcings.geothermalflux', 'NaN', 1, 'Inf', 1, 'timeseries', 1, '>=', 0)

        return md
    # }}}

    def marshall(self, prefix, md, fid):  # {{{
        yts = md.constants.yts

        WriteData(fid, prefix, 'name', 'md.basalforcings.model', 'data', 2, 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'perturbation_melting_rate', 'name', 'md.basalforcings.perturbation_melting_rate', 'format', 'DoubleMat', 'mattype', 1, 'scale', 1. / yts, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', yts)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'groundedice_melting_rate', 'name', 'md.basalforcings.groundedice_melting_rate', 'format', 'DoubleMat', 'mattype', 1, 'scale', 1. / yts, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', yts)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'geothermalflux', 'name', 'md.basalforcings.geothermalflux', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', yts)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'deepwater_melting_rate', 'format', 'DoubleMat', 'mattype', 3, 'timeserieslength', 2, 'name', 'md.basalforcings.deepwater_melting_rate', 'scale', 1. / yts, 'yts', yts)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'deepwater_elevation', 'format', 'DoubleMat', 'mattype', 3, 'name', 'md.basalforcings.deepwater_elevation', 'yts', yts)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'upperwater_melting_rate', 'format', 'DoubleMat', 'mattype', 3, 'timeserieslength', 2, 'name', 'md.basalforcings.upperwater_melting_rate', 'scale', 1. / yts, 'yts', yts)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'upperwater_elevation', 'format', 'DoubleMat', 'mattype', 3, 'name', 'md.basalforcings.upperwater_elevation', 'yts', yts)
    # }}}
