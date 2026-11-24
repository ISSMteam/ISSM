from checkfield import checkfield
from fielddisplay import fielddisplay
#from project3d import project3d # Uncomment if/when extrude is implemented
from WriteData import WriteData


class calvingcrevassedepth(object):
    """CALVINGCREVASSEDEPTH class definition

    Usage:
        calvingcrevassedepth = calvingcrevassedepth()
    """

    def __init__(self):  # {{{
        self.crevasse_opening_stress = 1
        self.crevasse_threshold      = 1.
        self.surface_hydrology_type  = 1.
        self.water_height            = 0.
        self.hurst                   = 0.
        self.sigma                   = 0.
        self.melt_supply             = 0.

        #self.setdefaultparameters() # Uncomment if/when setdefaultparameters is used
    # }}}
    def __repr__(self):  # {{{
        s = '   Calving Pi parameters:'
        s += '{}\n'.format(fielddisplay(self, 'crevasse_opening_stress', '0: stress only in the ice-flow direction, 1: max principal'))
        s += '{}\n'.format(fielddisplay(self, 'crevasse_threshold', 'ratio of full thickness to calve (e.g. 0.75 is for 75% of the total ice thickness)'))
        s += '{}\n'.format(fielddisplay(self, 'surface_hydrology_type', '0: default crevasse depth calving, 1: use mean water depth param to prescibe water_height, 2: use mean lake depth param to prescribe as water_height'))
        s += '{}\n'.format(fielddisplay(self, 'water_height', 'water height in the crevasse [m]'))
        s += '{}\n'.format(fielddisplay(self, 'hurst', 'Hurst exponent'))
        s += '{}\n'.format(fielddisplay(self, 'sigma', 'Roughness Amplitude [m]'))
        s += '{}\n'.format(fielddisplay(self, 'melt_supply', 'Amount of supply meltwater able to be distributed across the surface [m]'))
        return s
    # }}}
    def setdefaultparameters(self):  # {{{
        return self
    # }}}
    def extrude(self, md):  # {{{
        return self
    # }}}
    def checkconsistency(self, md, solution, analyses):  # {{{
	#Early return
        if solution != 'TransientSolution' or not md.transient.ismovingfront:
            return md

        md = checkfield(md, 'fieldname', 'calving.crevasse_opening_stress', 'numel', [1], 'values', [0,1])
        md = checkfield(md, 'fieldname', 'calving.crevasse_threshold', 'numel', [1], '>', 0., '<=', 1.)
        md = checkfield(md, 'fieldname', 'calving.surface_hydrology_type', 'numel', [1], 'values', [0,1,2])
        if self.surface_hydrology_type ==0:
            md = checkfield(md, 'fieldname', 'calving.water_height', 'NaN', 1, 'Inf', 1, 'timeseries', 1, '>=', 0)
        if self.surface_hydrology_type ==1:
            md = checkfield(md, 'fieldname', 'calving.hurst', 'NaN', 1, 'size', [md.mesh.numberofvertices, 1],'>=',0,'<=',1)
            md = checkfield(md, 'fieldname', 'calving.sigma', 'NaN', 1, 'size', [md.mesh.numberofvertices, 1],'>=',0)
            md = checkfield(md, 'fieldname', 'calving.melt_supply', 'NaN', 1, 'Inf', 1, 'timeseries', 1, '>=',0)
        if self.surface_hydrology_type ==2:
            md = checkfield(md, 'fieldname', 'calving.hurst', 'NaN', 1, 'size', [md.md.mesh.numberofvertices, 1],'>=',0,'<=',1)
            md = checkfield(md, 'fieldname', 'calving.sigma', 'NaN', 1, 'size', [md.mesh.numberofvertices, 1],'>=',0)
            md = checkfield(md, 'fieldname', 'calving.melt_supply', 'NaN', 1, 'Inf', 1, 'timeseries', 1, '>=',0)

        return md
    # }}}
    def marshall(self, prefix, md, fid):  # {{{
        yts = md.constants.yts
        WriteData(fid, prefix, 'name', 'md.calving.law', 'data', 6, 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'crevasse_opening_stress', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'crevasse_threshold', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'surface_hydrology_type', 'format', 'Integer')
        if self.surface_hydrology_type ==0:
            WriteData(fid, prefix, 'object', self, 'fieldname', 'water_height', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', md.constants.yts)
        if self.surface_hydrology_type ==1:
            WriteData(fid, prefix, 'object', self, 'fieldname', 'hurst', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', md.constants.yts)
            WriteData(fid, prefix, 'object', self, 'fieldname', 'sigma', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', md.constants.yts)
            WriteData(fid, prefix, 'object', self, 'fieldname', 'melt_supply', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', md.constants.yts)
        if self.surface_hydrology_type ==2:
                WriteData(fid, prefix, 'object', self, 'fieldname', 'hurst', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', md.constants.yts)
                WriteData(fid, prefix, 'object', self, 'fieldname', 'sigma', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', md.constants.yts)
                WriteData(fid, prefix, 'object', self, 'fieldname', 'melt_supply', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', md.constants.yts)
    # }}}
