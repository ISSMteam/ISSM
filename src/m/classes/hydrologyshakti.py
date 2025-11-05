import numpy as np
from fielddisplay import fielddisplay
from checkfield import checkfield
from WriteData import WriteData


class hydrologyshakti(object):
    """hydrologyshakti class definition

    Usage:
        hydrologyshakti = hydrologyshakti()
    """

    def __init__(self):  # {{{
        self.head = np.nan
        self.gap_height = np.nan
        self.gap_height_min  = 1e-3;
        self.gap_height_max  = 1.;
        self.bump_spacing = np.nan
        self.bump_height = np.nan
        self.englacial_input = np.nan
        self.moulin_input = np.nan
        self.reynolds = np.nan
        self.spchead = np.nan
        self.neumannflux = np.nan
        self.relaxation = 0
        self.storage = np.nan
        self.requested_outputs = []

    #set defaults
        self.setdefaultparameters()

    # }}}
    def __repr__(self):  # {{{
        s = '   hydrologyshakti solution parameters:'
        s += '{}\n'.format(fielddisplay(self, 'head', 'subglacial hydrology water head (m)'))
        s += '{}\n'.format(fielddisplay(self, 'gap_height', 'height of gap separating ice to bed (m)'))
        s += '{}\n'.format(fielddisplay(self, 'gap_height_min', 'minimum allowed gap height (m)'))
        s += '{}\n'.format(fielddisplay(self, 'gap_height_max', 'minimum allowed gap height (m)'))
        s += '{}\n'.format(fielddisplay(self, 'bump_spacing', 'characteristic bedrock bump spacing (m)'))
        s += '{}\n'.format(fielddisplay(self, 'bump_height', 'characteristic bedrock bump height (m)'))
        s += '{}\n'.format(fielddisplay(self, 'englacial_input', 'liquid water input from englacial to subglacial system (m / yr)'))
        s += '{}\n'.format(fielddisplay(self, 'moulin_input', 'liquid water input from moulins (at the vertices) to subglacial system (m^3 / s)'))
        s += '{}\n'.format(fielddisplay(self, 'reynolds', 'Reynolds'' number'))
        s += '{}\n'.format(fielddisplay(self, 'neumannflux', 'water flux applied along the model boundary (m^2 / s)'))
        s += '{}\n'.format(fielddisplay(self, 'spchead', 'water head constraints (NaN means no constraint) (m)'))
        s += '{}\n'.format(fielddisplay(self, 'relaxation', 'under - relaxation coefficient for nonlinear iteration'))
        s += '{}\n'.format(fielddisplay(self, 'storage', 'englacial storage coefficient (void ratio)'))
        s += '{}\n'.format(fielddisplay(self, 'requested_outputs', 'additional outputs requested'))
        return s
    # }}}

    def extrude(self, md):  # {{{
        return self
    # }}}

    def setdefaultparameters(self):  # {{{
        # Set under-relaxation parameter to be 1 (no under-relaxation of nonlinear iteration)
        self.gap_height_min  = 1e-3;
        self.gap_height_max  = 1.;
        self.relaxation = 1
        self.storage = 0
        self.requested_outputs = ['default']
        return self
    # }}}

    def defaultoutputs(self, md):  # {{{
        list = ['HydrologyHead', 'HydrologyGapHeight', 'EffectivePressure', 'HydrologyBasalFlux', 'DegreeOfChannelization']
        return list
    # }}}

    def checkconsistency(self, md, solution, analyses):  # {{{
        # Early return
        if 'HydrologyShaktiAnalysis' not in analyses:
            return md

        md = checkfield(md, 'fieldname', 'hydrology.head', 'size', [md.mesh.numberofvertices], 'NaN', 1, 'Inf', 1)
        md = checkfield(md, 'fieldname', 'hydrology.gap_height', '>=', 0, 'size', [md.mesh.numberofelements], 'NaN', 1, 'Inf', 1)
        md = checkfield(md, 'fieldname', 'hydrology.gap_height_min', '>=', 0, 'numel', 1, 'NaN', 1, 'Inf', 1)
        md = checkfield(md, 'fieldname', 'hydrology.gap_height_max', '>=', 0, 'numel', 1, 'NaN', 1, 'Inf', 1)
        md = checkfield(md, 'fieldname', 'hydrology.bump_spacing', '>', 0, 'size', [md.mesh.numberofelements], 'NaN', 1, 'Inf', 1)
        md = checkfield(md, 'fieldname', 'hydrology.bump_height', '>=', 0, 'size', [md.mesh.numberofelements], 'NaN', 1, 'Inf', 1)
        md = checkfield(md, 'fieldname', 'hydrology.englacial_input', '>=', 0, 'NaN', 1, 'Inf', 1, 'timeseries', 1)
        md = checkfield(md, 'fieldname', 'hydrology.moulin_input', '>=', 0, 'NaN', 1, 'Inf', 1, 'timeseries', 1)
        md = checkfield(md, 'fieldname', 'hydrology.reynolds', '>', 0, 'size', [md.mesh.numberofelements], 'NaN', 1, 'Inf', 1)
        md = checkfield(md, 'fieldname', 'hydrology.neumannflux', 'timeseries', 1, 'NaN', 1, 'Inf', 1)
        md = checkfield(md, 'fieldname', 'hydrology.spchead', 'Inf', 1, 'timeseries', 1)
        md = checkfield(md, 'fieldname', 'hydrology.relaxation', '>=', 0)
        md = checkfield(md, 'fieldname', 'hydrology.storage', '>=', 0, 'universal', 1, 'NaN', 1, 'Inf', 1)
        md = checkfield(md, 'fieldname', 'hydrology.requested_outputs', 'stringrow', 1)
        return md
    # }}}

    def marshall(self, prefix, md, fid):  # {{{
        yts = md.constants.yts

        WriteData(fid, prefix, 'name', 'md.hydrology.model', 'data', 3, 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'class', 'hydrology', 'fieldname', 'head', 'format', 'DoubleMat', 'mattype', 1)
        WriteData(fid, prefix, 'object', self, 'class', 'hydrology', 'fieldname', 'gap_height', 'format', 'DoubleMat', 'mattype', 2)
        WriteData(fid, prefix, 'object', self, 'class', 'hydrology', 'fieldname', 'gap_height_min', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'class', 'hydrology', 'fieldname', 'gap_height_max', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'class', 'hydrology', 'fieldname', 'bump_spacing', 'format', 'DoubleMat', 'mattype', 2)
        WriteData(fid, prefix, 'object', self, 'class', 'hydrology', 'fieldname', 'bump_height', 'format', 'DoubleMat', 'mattype', 2)
        WriteData(fid, prefix, 'object', self, 'class', 'hydrology', 'fieldname', 'englacial_input', 'format', 'DoubleMat', 'mattype', 1, 'scale', 1. / yts, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', yts)
        WriteData(fid, prefix, 'object', self, 'class', 'hydrology', 'fieldname', 'moulin_input', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', yts)
        WriteData(fid, prefix, 'object', self, 'class', 'hydrology', 'fieldname', 'reynolds', 'format', 'DoubleMat', 'mattype', 2)
        WriteData(fid, prefix, 'object', self, 'class', 'hydrology', 'fieldname', 'neumannflux', 'format', 'DoubleMat', 'mattype', 2, 'timeserieslength', md.mesh.numberofelements + 1, 'yts', yts)
        WriteData(fid, prefix, 'object', self, 'class', 'hydrology', 'fieldname', 'spchead', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', md.constants.yts)
        WriteData(fid, prefix, 'object', self, 'class', 'hydrology', 'fieldname', 'relaxation', 'format', 'Double')
        # NOTE: We first have to check if we have a NumPy array here
        if np.size(self.storage)==1 or (((np.shape(self.storage)[0] == md.mesh.numberofvertices) or (np.shape(self.storage)[0] == md.mesh.numberofvertices + 1)) or ((len(np.shape(self.storage)) == 2) and (np.shape(self.storage)[0] == md.mesh.numberofelements) and (np.shape(self.storage)[1] > 1))):
            mattype = 1
            tsl = md.mesh.numberofvertices
        else:
            mattype = 2
            tsl = md.mesh.numberofelements
        WriteData(fid, prefix, 'object', self, 'class', 'hydrology', 'fieldname', 'storage', 'format', 'DoubleMat', 'mattype', mattype, 'timeserieslength', tsl + 1, 'yts', md.constants.yts)

        # Process requested outputs
        outputs = self.requested_outputs
        indices = [i for i, x in enumerate(outputs) if x == 'default']
        if len(indices) > 0:
            outputscopy = outputs[0:max(0, indices[0] - 1)] + self.defaultoutputs(md) + outputs[indices[0] + 1:]
            outputs = outputscopy
        WriteData(fid, prefix, 'data', outputs, 'name', 'md.hydrology.requested_outputs', 'format', 'StringArray')
    # }}}
