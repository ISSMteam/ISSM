from checkfield import checkfield
from fielddisplay import fielddisplay
#from project3d import project3d # Uncomment if/when extrude is implemented
from WriteData import WriteData


class calvingcrevassedepth(object):
    """CALVINCREVASSEDEPTH class definition

    Usage:
        calvingcrevassedepth = calvingcrevassedepth()
    """

    def __init__(self):  # {{{
        self.crevasse_opening_stress = 1
        self.crevasse_threshold      = 1.
        self.water_height            = 0.

        #self.setdefaultparameters() # Uncomment if/when setdefaultparameters is used
    # }}}
    def __repr__(self):  # {{{
        s = '   Calving Pi parameters:'
        s += '{}\n'.format(fielddisplay(self, 'crevasse_opening_stress', '0: stress only in the ice-flow direction, 1: max principal'))
        s += '{}\n'.format(fielddisplay(self, 'crevasse_threshold', 'ratio of full thickness to calve (e.g. 0.75 is for 75% of the total ice thickness)'))
        s += '{}\n'.format(fielddisplay(self, 'water_height', 'water height in the crevasse [m]'))
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
        md = checkfield(md, 'fieldname', 'calving.water_height', 'NaN', 1, 'Inf', 1, 'timeseries', 1, '>=', 0) 

        return md
    # }}}
    def marshall(self, prefix, md, fid):  # {{{
        yts = md.constants.yts
        WriteData(fid, prefix, 'name', 'md.calving.law', 'data', 6, 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'crevasse_opening_stress', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'crevasse_threshold', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'water_height', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', md.constants.yts)
    # }}}
