import numpy as np

from checkfield import checkfield
from fielddisplay import fielddisplay
from project3d import project3d
from WriteData import WriteData


class basalforcings(object):
    """BASAL FORCINGS class definition

    Usage:
        basalforcings = basalforcings()
    """

    def __init__(self):  # {{{
        self.groundedice_melting_rate = np.nan
        self.floatingice_melting_rate = np.nan
        self.perturbation_melting_rate = np.nan
        self.geothermalflux = np.nan

        self.setdefaultparameters()
    # }}}
    def __repr__(self):  # {{{
        s = '   basal forcings parameters:\n'
        s += '{}\n'.format(fielddisplay(self, 'groundedice_melting_rate', 'basal melting rate (positive if melting) [m/yr]'))
        s += '{}\n'.format(fielddisplay(self, 'floatingice_melting_rate', 'basal melting rate (positive if melting) [m/yr]'))
        s += '{}\n'.format(fielddisplay(self, 'perturbation_melting_rate', '(optional) perturbation in basal melting rate under floating ice [m/yr]'))
        s += '{}\n'.format(fielddisplay(self, 'geothermalflux', 'geothermal heat flux [W/m^2]'))
        return s
    # }}}
    def extrude(self, md):  # {{{
        self.groundedice_melting_rate = project3d(md, 'vector', self.groundedice_melting_rate, 'type', 'node', 'layer', 1)
        self.perturbation_melting_rate = project3d(md, 'vector', self.perturbation_melting_rate, 'type', 'node', 'layer', 1)
        self.floatingice_melting_rate = project3d(md, 'vector', self.floatingice_melting_rate, 'type', 'node', 'layer', 1)
        self.geothermalflux = project3d(md, 'vector', self.geothermalflux, 'type', 'node', 'layer', 1) # Bedrock only gets geothermal flux
        return self
    # }}}
    def initialize(self, md):  # {{{
        if np.all(np.isnan(self.groundedice_melting_rate)):
            self.groundedice_melting_rate = np.zeros((md.mesh.numberofvertices, ))
            print('      no basalforcings.groundedice_melting_rate specified: values set as zero')
        if np.all(np.isnan(self.floatingice_melting_rate)):
            self.floatingice_melting_rate = np.zeros((md.mesh.numberofvertices, ))
            print('      no basalforcings.floatingice_melting_rate specified: values set as zero')
        return self
    # }}}
    def setdefaultparameters(self):  # {{{
        return self
    # }}}
    def checkconsistency(self, md, solution, analyses):  # {{{
        if 'MasstransportAnalysis' in analyses and not solution == 'TransientSolution' and not md.transient.ismasstransport:
            md = checkfield(md, 'fieldname', 'basalforcings.groundedice_melting_rate', 'NaN', 1, 'Inf', 1, 'timeseries', 1)
            md = checkfield(md, 'fieldname', 'basalforcings.floatingice_melting_rate', 'NaN', 1, 'Inf', 1, 'timeseries', 1)
        if 'BalancethicknessAnalysis' in analyses:
            md = checkfield(md, 'fieldname', 'basalforcings.groundedice_melting_rate', 'NaN', 1, 'Inf', 1, 'size', [md.mesh.numberofvertices])
            md = checkfield(md, 'fieldname', 'basalforcings.floatingice_melting_rate', 'NaN', 1, 'Inf', 1, 'size', [md.mesh.numberofvertices])
        if 'ThermalAnalysis' in analyses and not solution == 'TransientSolution' and not md.transient.isthermal:
            md = checkfield(md, 'fieldname', 'basalforcings.groundedice_melting_rate', 'NaN', 1, 'Inf', 1, 'timeseries', 1)
            md = checkfield(md, 'fieldname', 'basalforcings.floatingice_melting_rate', 'NaN', 1, 'Inf', 1, 'timeseries', 1)
            md = checkfield(md, 'fieldname', 'basalforcings.geothermalflux', 'NaN', 1, 'Inf', 1, 'timeseries', 1, '>=', 0)
        return md
    # }}}
    def marshall(self, prefix, md, fid):  # {{{
        yts = md.constants.yts
        WriteData(fid, prefix, 'name', 'md.basalforcings.model', 'data', 1, 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'groundedice_melting_rate', 'format', 'DoubleMat', 'mattype', 1, 'scale', 1 / yts, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', yts)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'floatingice_melting_rate', 'format', 'DoubleMat', 'mattype', 1, 'scale', 1 / yts, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', yts)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'geothermalflux', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', yts)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'perturbation_melting_rate', 'format', 'DoubleMat', 'name', 'md.basalforcings.perturbation_melting_rate', 'mattype', 1, 'scale', 1 / yts, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', yts)
    # }}}
