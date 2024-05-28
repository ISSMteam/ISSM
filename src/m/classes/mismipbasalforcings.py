import numpy as np

from checkfield import checkfield
from fielddisplay import fielddisplay
from project3d import project3d
from WriteData import WriteData


class mismipbasalforcings(object):
    """MISMIP Basal Forcings class definition

    Usage:
        mismipbasalforcings = mismipbasalforcings()
    """

    def __init__(self):  # {{{
        self.groundedice_melting_rate = np.nan
        self.meltrate_factor          = np.nan
        self.threshold_thickness      = 0.
        self.upperdepth_melt          = 0.
        self.geothermalflux           = np.nan
        self.setdefaultparameters()
    # }}}
    def __repr__(self):  # {{{
        s = '   MISMIP + basal melt parameterization\n'
        s += '{}\n'.format(fielddisplay(self, "groundedice_melting_rate", "basal melting rate (positive if melting) [m / yr]"))
        s += '{}\n'.format(fielddisplay(self, "meltrate_factor", "Melt - rate rate factor [1 / yr] (sign is opposite to MISMIP + benchmark to remain consistent with ISSM convention of positive values for melting)"))
        s += '{}\n'.format(fielddisplay(self, "threshold_thickness", "Threshold thickness for saturation of basal melting [m]"))
        s += '{}\n'.format(fielddisplay(self, "upperdepth_melt", "Depth above which melt rate is zero [m]"))
        s += '{}\n'.format(fielddisplay(self, "geothermalflux", "Geothermal heat flux [W / m^2]"))
        return s
    # }}}
    def extrude(self, md):  # {{{
        self.groundedice_melting_rate = project3d(md, 'vector', self.groundedice_melting_rate, 'type', 'node', 'layer', 1)
        self.geothermalflux = project3d(md, 'vector', self.geothermalflux, 'type', 'node', 'layer', 1)  #bedrock only gets geothermal flux
        return self
    # }}}
    def initialize(self, md):  # {{{
        if np.all(np.isnan(self.groundedice_melting_rate)):
            self.groundedice_melting_rate = np.zeros((md.mesh.numberofvertices))
            print('no basalforcings.groundedice_melting_rate specified: values set as zero')
        if np.all(np.isnan(self.geothermalflux)):
            self.geothermalflux = np.zeros((md.mesh.numberofvertices))
            print("      no basalforcings.geothermalflux specified: values set as zero")
        return self
    # }}}
    def setdefaultparameters(self):  # {{{
        # default values for melting parameterization
        self.meltrate_factor = 0.2
        self.threshold_thickness = 75.
        self.upperdepth_melt = -100.
        return self
    # }}}
    def checkconsistency(self, md, solution, analyses):  # {{{
        # Early return
        if 'MasstransportAnalysis' in analyses and not solution == 'TransientSolution' and not md.transient.ismasstransport:
            md = checkfield(md, 'fieldname', 'basalforcings.groundedice_melting_rate', 'NaN', 1, 'Inf', 1, 'timeseries', 1)
            md = checkfield(md, 'fieldname', 'basalforcings.meltrate_factor', '>=', 0, 'numel', [1])
            md = checkfield(md, 'fieldname', 'basalforcings.threshold_thickness', '>=', 0, 'numel', [1])
            md = checkfield(md, 'fieldname', 'basalforcings.upperdepth_melt', '<=', 0, 'numel', [1])
        if 'BalancethicknessAnalysis' in analyses:
            md = checkfield(md, 'fieldname', 'basalforcings.groundedice_melting_rate', 'NaN', 1, 'Inf', 1, 'size', [md.mesh.numberofvertices])
            md = checkfield(md, 'fieldname', 'basalforcings.meltrate_factor', '>=', 0, 'numel', [1])
            md = checkfield(md, 'fieldname', 'basalforcings.threshold_thickness', '>=', 0, 'numel', [1])
            md = checkfield(md, 'fieldname', 'basalforcings.upperdepth_melt', '<=', 0, 'numel', [1])
        if 'ThermalAnalysis' in analyses and not (solution == 'TransientSolution' and not md.transient.isthermal):
            md = checkfield(md, 'fieldname', 'basalforcings.groundedice_melting_rate', 'NaN', 1, 'Inf', 1, 'timeseries', 1)
            md = checkfield(md, 'fieldname', 'basalforcings.meltrate_factor', '>=', 0, 'numel', [1])
            md = checkfield(md, 'fieldname', 'basalforcings.threshold_thickness', '>=', 0, 'numel', [1])
            md = checkfield(md, 'fieldname', 'basalforcings.upperdepth_melt', '<=', 0, 'numel', [1])
            md = checkfield(md, 'fieldname', 'basalforcings.geothermalflux', 'NaN', 1, 'Inf', 1, 'timeseries', 1, '>=', 0)
        return md
    # }}}
    def marshall(self, prefix, md, fid):  # {{{
        yts = md.constants.yts
        if yts != 365.2422 * 24. * 3600.:
            print('WARNING: value of yts for MISMIP + runs different from ISSM default!')
        WriteData(fid, prefix, 'name', 'md.basalforcings.model', 'data', 3, 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'groundedice_melting_rate', 'format', 'DoubleMat', 'name', 'md.basalforcings.groundedice_melting_rate', 'mattype', 1, 'scale', 1. / yts, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', yts)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'geothermalflux', 'name', 'md.basalforcings.geothermalflux', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', yts)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'meltrate_factor', 'format', 'Double', 'scale', 1. / yts)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'threshold_thickness', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'upperdepth_melt', 'format', 'Double')
    # }}}
