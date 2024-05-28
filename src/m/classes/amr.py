from fielddisplay import fielddisplay
from checkfield import checkfield
from WriteData import WriteData


class amr(object):
    """AMR Class definition

    Usage:
        amr = amr()
    """

    def __init__(self):  # {{{
        self.hmin = 0
        self.hmax = 0
        self.fieldname = ''
        self.err = 0
        self.keepmetric = 0
        self.gradation = 0
        self.groundingline_resolution = 0
        self.groundingline_distance = 0
        self.icefront_resolution = 0
        self.icefront_distance = 0
        self.thicknesserror_resolution = 0
        self.thicknesserror_threshold = 0
        self.thicknesserror_groupthreshold = 0
        self.thicknesserror_maximum = 0
        self.deviatoricerror_resolution = 0
        self.deviatoricerror_threshold = 0
        self.deviatoricerror_groupthreshold = 0
        self.deviatoricerror_maximum = 0
        self.restart = 0

        self.setdefaultparameters()
    # }}}

    def __repr__(self):  # {{{
        s = '   amr parameters:\n'
        s += '{}\n'.format(fielddisplay(self, 'hmin', 'minimum element length'))
        s += '{}\n'.format(fielddisplay(self, 'hmax', 'maximum element length'))
        s += '{}\n'.format(fielddisplay(self, 'fieldname', 'name of input that will be used to compute the metric (should be an input of FemModel)'))
        s += '{}\n'.format(fielddisplay(self, 'keepmetric', 'indicates whether the metric should be kept every remeshing time'))
        s += '{}\n'.format(fielddisplay(self, 'gradation', 'maximum ratio between two adjacent edges'))
        s += '{}\n'.format(fielddisplay(self, 'groundingline_resolution', 'element length near the grounding line'))
        s += '{}\n'.format(fielddisplay(self, 'groundingline_distance', 'distance around the grounding line which elements will be refined'))
        s += '{}\n'.format(fielddisplay(self, 'icefront_resolution', 'element length near the ice front'))
        s += '{}\n'.format(fielddisplay(self, 'icefront_distance', 'distance around the ice front which elements will be refined'))
        s += '{}\n'.format(fielddisplay(self, 'thicknesserror_resolution', 'element length when thickness error estimator is used'))
        s += '{}\n'.format(fielddisplay(self, 'thicknesserror_threshold', 'maximum threshold thickness error permitted'))
        s += '{}\n'.format(fielddisplay(self, 'thicknesserror_groupthreshold', 'maximum group threshold thickness error permitted'))
        s += '{}\n'.format(fielddisplay(self, 'thicknesserror_maximum', 'maximum thickness error permitted'))
        s += '{}\n'.format(fielddisplay(self, 'deviatoricerror_resolution', 'element length when deviatoric stress error estimator is used'))
        s += '{}\n'.format(fielddisplay(self, 'deviatoricerror_threshold', 'maximum threshold deviatoricstress error permitted'))
        s += '{}\n'.format(fielddisplay(self, 'deviatoricerror_groupthreshold', 'maximum group threshold deviatoric stress error permitted'))
        s += '{}\n'.format(fielddisplay(self, 'deviatoricerror_maximum', 'maximum deviatoricstress error permitted'))
        s += '{}\n'.format(fielddisplay(self, 'restart', 'indicates if ReMesh() will call before first time step'))
        return s
    # }}}

    def setdefaultparameters(self):  # {{{
        self.hmin = 100
        self.hmax = 100e3

        # Fields
        self.fieldname = 'Vel'
        self.err = 3

        # Keep metric?
        self.keepmetric = 1

        # Control of element lengths
        self.gradation = 1.5

        # Other criteria
        self.groundingline_resolution = 500
        self.groundingline_distance = 0
        self.icefront_resolution = 500
        self.icefront_distance = 0
        self.thicknesserror_resolution = 500
        self.thicknesserror_threshold = 0
        self.thicknesserror_groupthreshold = 0
        self.thicknesserror_maximum = 0
        self.deviatoricerror_resolution = 500
        self.deviatoricerror_threshold = 0
        self.deviatoricerror_groupthreshold = 0
        self.deviatoricerror_maximum = 0

        # Is restart? This calls femmodel->ReMesh() before first time step.
        self.restart = 0
        return self
    # }}}

    def checkconsistency(self, md, solution, analyses):  # {{{
        md = checkfield(md, 'fieldname', 'amr.hmax', 'numel', [1], '>', 0, 'NaN', 1)
        md = checkfield(md, 'fieldname', 'amr.hmin', 'numel', [1], '>', 0, '<', self.hmax, 'NaN', 1)
        md = checkfield(md, 'fieldname', 'amr.keepmetric', 'numel', [1], '>=', 0, '<=', 1, 'NaN', 1)
        md = checkfield(md, 'fieldname', 'amr.gradation', 'numel', [1], '>=', 1.1, '<=', 5, 'NaN', 1)
        md = checkfield(md, 'fieldname', 'amr.groundingline_resolution', 'numel', [1], '>', 0, '<', self.hmax, 'NaN', 1)
        md = checkfield(md, 'fieldname', 'amr.groundingline_distance', 'numel', [1], '>=', 0, 'NaN', 1, 'Inf', 1)
        md = checkfield(md, 'fieldname', 'amr.icefront_resolution', 'numel', [1], '>', 0, '<', self.hmax, 'NaN', 1)
        md = checkfield(md, 'fieldname', 'amr.icefront_distance', 'numel', [1], '>=', 0, 'NaN', 1, 'Inf', 1)
        md = checkfield(md, 'fieldname', 'amr.thicknesserror_resolution', 'numel', [1], '>', 0, '<', self.hmax, 'NaN', 1)
        md = checkfield(md, 'fieldname', 'amr.thicknesserror_threshold', 'numel', [1], '>=', 0, '<=', 1, 'NaN', 1)
        md = checkfield(md, 'fieldname', 'amr.thicknesserror_groupthreshold', 'numel', [1], '>=', 0, '<=', 1, 'NaN', 1)
        md = checkfield(md, 'fieldname', 'amr.thicknesserror_maximum', 'numel', [1], '>=', 0, 'NaN', 1, 'Inf', 1)
        md = checkfield(md, 'fieldname', 'amr.deviatoricerror_resolution', 'numel', [1], '>', 0, '<', self.hmax, 'NaN', 1)
        md = checkfield(md, 'fieldname', 'amr.deviatoricerror_threshold', 'numel', [1], '>=', 0, '<=', 1, 'NaN', 1)
        md = checkfield(md, 'fieldname', 'amr.deviatoricerror_groupthreshold', 'numel', [1], '>=', 0, '<=', 1, 'NaN', 1)
        md = checkfield(md, 'fieldname', 'amr.deviatoricerror_maximum', 'numel', [1], '>=', 0, 'NaN', 1, 'Inf', 1)
        md = checkfield(md, 'fieldname', 'amr.restart', 'numel', [1], '>=', 0, '<=', 1, 'NaN', 1)
        return md
   # }}}

    def marshall(self, prefix, md, fid):  # {{{
        WriteData(fid, prefix, 'name', 'md.amr.type', 'data', 1, 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'hmin', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'hmax', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'fieldname', 'format', 'String')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'err', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'keepmetric', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'gradation', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'groundingline_resolution', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'groundingline_distance', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'icefront_resolution', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'icefront_distance', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'thicknesserror_resolution', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'thicknesserror_threshold', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'thicknesserror_groupthreshold', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'thicknesserror_maximum', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'deviatoricerror_resolution', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'deviatoricerror_threshold', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'deviatoricerror_groupthreshold', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'deviatoricerror_maximum', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'class', 'amr', 'fieldname', 'restart', 'format', 'Integer')
    # }}}
