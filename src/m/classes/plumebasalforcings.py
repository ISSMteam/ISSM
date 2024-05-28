import numpy as np

from checkfield import checkfield
from fielddisplay import fielddisplay
from project3d import project3d
from structtoobj import structtoobj
from WriteData import WriteData


class plumebasalforcings(object):
    """PLUME BASAL FORCINGS class definition

    Usage:
        plumebasalforcings = plumebasalforcings()
    """

    def __init__(self, *args):  # {{{
        self.floatingice_melting_rate   = np.nan
        self.groundedice_melting_rate   = np.nan
        self.mantleconductivity         = np.nan
        self.nusselt                    = np.nan
        self.dtbg                       = np.nan
        self.plumeradius                = np.nan
        self.topplumedepth              = np.nan
        self.bottomplumedepth           = np.nan
        self.plumex                     = np.nan
        self.plumey                     = np.nan
        self.crustthickness             = np.nan
        self.uppercrustthickness        = np.nan
        self.uppercrustheat             = np.nan
        self.lowercrustheat             = np.nan

        nargs = len(args)
        if nargs == 0:
            self.setdefaultparameters()
        elif nargs == 1:
            # TODO: This has not been tested
            self = structtoobj(self, args[0])
        else:
            error('constuctor not supported')
    # }}}

    def __repr__(self):  # {{{
        s = '   mantle plume basal melt parameterization:\n'
        s += '{}\n'.format(fielddisplay(self, 'groundedice_melting_rate', 'basal melting rate (positive if melting) [m/yr]'))
        s += '{}\n'.format(fielddisplay(self, 'floatingice_melting_rate', 'basal melting rate (positive if melting) [m/yr]'))
        s += '{}\n'.format(fielddisplay(self, 'mantleconductivity', 'mantle heat conductivity [W/m^3]'))
        s += '{}\n'.format(fielddisplay(self, 'nusselt', 'nusselt number, ratio of mantle to plume [1]'))
        s += '{}\n'.format(fielddisplay(self, 'dtbg', 'background temperature gradient [degree/m]'))
        s += '{}\n'.format(fielddisplay(self, 'plumeradius', 'radius of the mantle plume [m]'))
        s += '{}\n'.format(fielddisplay(self, 'topplumedepth', 'depth of the mantle plume top below the crust [m]'))
        s += '{}\n'.format(fielddisplay(self, 'bottomplumedepth', 'depth of the mantle plume base below the crust [m]'))
        s += '{}\n'.format(fielddisplay(self, 'plumex', 'x coordinate of the center of the plume [m]'))
        s += '{}\n'.format(fielddisplay(self, 'plumey', 'y coordinate of the center of the plume [m]'))
        s += '{}\n'.format(fielddisplay(self, 'crustthickness', 'thickness of the crust [m]'))
        s += '{}\n'.format(fielddisplay(self, 'uppercrustthickness', 'thickness of the upper crust [m]'))
        s += '{}\n'.format(fielddisplay(self, 'uppercrustheat', 'volumic heat of the upper crust [w/m^3]'))
        s += '{}\n'.format(fielddisplay(self, 'lowercrustheat', 'volumic heat of the lowercrust [w/m^3]'))
        return s
    # }}}

    def initialize(self, md):  #{{{
        if np.all(np.isnan(self.groundedice_melting_rate)):
            self.groundedice_melting_rate = np.zeros((md.mesh.numberofvertices, ))
            print('      no basalforcings.groundedice_melting_rate specified: values set as zero')
        if np.all(np.isnan(self.floatingice_melting_rate)):
            self.floatingice_melting_rate = np.zeros((md.mesh.numberofvertices, ))
            print('      no basalforcings.floatingice_melting_rate specified: values set as zero')
        return self
    # }}}

    def extrude(self, md):  # {{{
        self.groundedice_melting_rate = project3d(md, 'vector', self.groundedice_melting_rate, 'type', 'node', 'layer', 1)
        self.floatingice_melting_rate = project3d(md, 'vector', self.floatingice_melting_rate, 'type', 'node', 'layer', 1)
        return self
    # }}}

    def setdefaultparameters(self):  # {{{
        # Default values for melting parameterization
        self.mantleconductivity = 2.2
        self.nusselt = 300
        self.dtbg = 11 / 1000.0
        self.plumeradius = 100000
        self.topplumedepth = 10000
        self.bottomplumedepth = 1050000
        self.crustthickness = 30000
        self.uppercrustthickness = 14000
        self.uppercrustheat = 1.7 * pow(10, -6)
        self.lowercrustheat = 0.4 * pow(10, -6)
        return self
    # }}}

    def checkconsistency(self, md, solution, analyses):  # {{{
        if 'MasstransportAnalysis' in analyses and not (solution == 'TransientSolution' and md.transient.ismasstransport == 0):
            md = checkfield(md, 'fieldname', 'basalforcings.groundedice_melting_rate', 'NaN', 1, 'timeseries', 1)
            md = checkfield(md, 'fieldname', 'basalforcings.floatingice_melting_rate', 'NaN', 1, 'timeseries', 1)
        if 'BalancethicknessAnalysis' in analyses:
            md = checkfield(md, 'fieldname', 'basalforcings.groundedice_melting_rate', 'NaN', 1, 'size', [md.mesh.numberofvertices])
            md = checkfield(md, 'fieldname', 'basalforcings.floatingice_melting_rate', 'NaN', 1, 'size', [md.mesh.numberofvertices])
        if 'ThermalAnalysis' in analyses and not (solution == 'TransientSolution' and md.transient.isthermal == 0):
            md = checkfield(md, 'fieldname', 'basalforcings.groundedice_melting_rate', 'NaN', 1, 'timeseries', 1)
            md = checkfield(md, 'fieldname', 'basalforcings.floatingice_melting_rate', 'NaN', 1, 'timeseries', 1)
            md = checkfield(md, 'fieldname', 'basalforcings.mantleconductivity', '>=', 0, 'numel', 1)
            md = checkfield(md, 'fieldname', 'basalforcings.nusselt', '>', 0, 'numel', 1)
            md = checkfield(md, 'fieldname', 'basalforcings.dtbg', '>', 0, 'numel', 1)
            md = checkfield(md, 'fieldname', 'basalforcings.topplumedepth', '>', 0, 'numel', 1)
            md = checkfield(md, 'fieldname', 'basalforcings.bottomplumedepth', '>', 0, 'numel', 1)
            md = checkfield(md, 'fieldname', 'basalforcings.plumex', 'numel', 1)
            md = checkfield(md, 'fieldname', 'basalforcings.plumey', 'numel', 1)
            md = checkfield(md, 'fieldname', 'basalforcings.crustthickness', '>', 0, 'numel', 1)
            md = checkfield(md, 'fieldname', 'basalforcings.uppercrustthickness', '>', 0, 'numel', 1)
            md = checkfield(md, 'fieldname', 'basalforcings.uppercrustheat', '>', 0, 'numel', 1)
            md = checkfield(md, 'fieldname', 'basalforcings.lowercrustheat', '>', 0, 'numel', 1)
        return md
    # }}}

    def marshall(self, prefix, md, fid):  # {{{
        yts = md.constants.yts

        WriteData(fid, prefix, 'name', 'md.basalforcings.model', 'data', 4, 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'floatingice_melting_rate', 'format', 'DoubleMat', 'name', 'md.basalforcings.floatingice_melting_rate', 'mattype', 1, 'scale', 1. / yts, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', yts)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'groundedice_melting_rate', 'format', 'DoubleMat', 'name', 'md.basalforcings.groundedice_melting_rate', 'mattype', 1, 'scale', 1. / yts, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', yts)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'mantleconductivity', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'nusselt', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'dtbg', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'plumeradius', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'topplumedepth', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'bottomplumedepth', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'plumex', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'plumey', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'crustthickness', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'uppercrustthickness', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'uppercrustheat', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'lowercrustheat', 'format', 'Double')
    # }}}
