import numpy as np
from checkfield import checkfield
from fielddisplay import fielddisplay
from project3d import project3d
from WriteData import WriteData


class hydrologyglads(object):
    """hydrologyglads class definition

    Usage:
        hydrologyglads = hydrologyglads()
    """

    def __init__(self, *args):  # {{{
       # Sheet
        self.pressure_melt_coefficient = 0.
        self.sheet_conductivity = np.nan
        self.cavity_spacing = 0.
        self.bump_height = np.nan
        self.omega = 0.
        self.sheet_alpha = np.nan 
        self.sheet_beta = np.nan 
        self.rheology_B_base = np.nan
        self.isincludesheetthickness = 0
        self.creep_open_flag = 1
        self.rheology_B_base = np.nan

        # Channels
        self.ischannels = 0
        self.channel_conductivity = np.nan
        self.channel_sheet_width = 0.
        self.channel_alpha = np.nan 
        self.channel_beta = np.nan 

        # Other
        self.spcphi = np.nan
        self.moulin_input = np.nan
        self.neumannflux = np.nan
        self.englacial_void_ratio = 0.
        self.requested_outputs = []
        self.melt_flag = 0
        self.istransition = 0

        nargs = len(args)
        if nargs == 0:
            self.setdefaultparameters()
        elif nargs == 1:
            # TODO: Replace the following with constructor
            self.setdefaultparameters()
        else:
            raise Exception('constructor not supported')
        # }}}

    def __repr__(self):  # {{{
        s = '   GlaDS (hydrologyglads) solution parameters:\n'
        s += '\t--SHEET\n'
        s += '{}\n'.format(fielddisplay(self, 'pressure_melt_coefficient', 'Pressure melt coefficient (c_t) [K Pa^ - 1]'))
        s += '{}\n'.format(fielddisplay(self, 'sheet_conductivity', 'sheet conductivity (k) [m^(7 / 4) kg^(- 1 / 2)]'))
        s += '{}\n'.format(fielddisplay(self, 'sheet_alpha', 'First sheet-flow exponent (alpha_s) []')) #TH
        s += '{}\n'.format(fielddisplay(self, 'sheet_beta', 'Second sheet-flow exponent (beta_s) []')) #TH
        s += '{}\n'.format(fielddisplay(self, 'cavity_spacing', 'cavity spacing (l_r) [m]'))
        s += '{}\n'.format(fielddisplay(self, 'bump_height', 'typical bump height (h_r) [m]'))
        s += '{}\n'.format(fielddisplay(self, 'omega', 'transition parameter (omega) []')) #TH
        s += '{}\n'.format(fielddisplay(self, 'rheology_B_base', 'ice rheology factor B at base of ice (B) [Pa s^(-1/3)]')) #SE
        s += '{}\n'.format(fielddisplay(self, 'isincludesheetthickness', 'Do we add rho_w*g*h in effective pressure calculation? 1: yes, 0: no'))
        s += '{}\n'.format(fielddisplay(self, 'creep_open_flag', 'Do we allow cavities to open by creep when N<0? 1: yes, 0: no'))
        s += '\t--CHANNELS\n'
        s += '{}\n'.format(fielddisplay(self, 'ischannels', 'Do we allow for channels? 1: yes, 0: no'))
        s += '{}\n'.format(fielddisplay(self, 'channel_conductivity', 'channel conductivity (k_c) [m^(3 / 2) kg^(- 1 / 2)]'))
        s += '{}\n'.format(fielddisplay(self, 'channel_sheet_width', 'channel sheet width [m]'))
        s += '{}\n'.format(fielddisplay(self, 'channel_alpha', 'First channel-flow exponent (alpha_s) []')) #TH
        s += '{}\n'.format(fielddisplay(self, 'channel_beta', 'Second channel-flow exponent (beta_s) []')) #TH
        s += '\t--OTHER\n'
        s += '{}\n'.format(fielddisplay(self, 'spcphi', 'Hydraulic potential Dirichlet constraints [Pa]'))
        s += '{}\n'.format(fielddisplay(self, 'neumannflux', 'water flux applied along the model boundary (m^2 / s)'))
        s += '{}\n'.format(fielddisplay(self, 'moulin_input', 'moulin input (Q_s) [m^3 / s]'))
        s += '{}\n'.format(fielddisplay(self, 'englacial_void_ratio', 'englacial void ratio (e_v)'))
        s += '{}\n'.format(fielddisplay(self, 'requested_outputs', 'additional outputs requested'))
        s += '{}\n'.format(fielddisplay(self, 'melt_flag', 'User specified basal melt? 0: no (default), 1: use md.basalforcings.groundedice_melting_rate'))
        s += '{}\n'.format(fielddisplay(self, 'istransition','do we use standard [0, default] or transition model [1]'))
        return s
    # }}}

    def defaultoutputs(self, md):  # {{{
        list = ['EffectivePressure', 'HydraulicPotential', 'HydrologySheetThickness', 'ChannelArea', 'ChannelDischarge']
        return list
    # }}}

    def extrude(self, md):  # {{{
        self.sheet_conductivity = project3d(md, 'vector', self.sheet_conductivity, 'type', 'node', 'layer', 1)
        self.bump_height = project3d(md, 'vector', self.bump_height, 'type', 'node', 'layer', 1)

        # Other
        self.spcphi = project3d(md, 'vector', self.spcphi, 'type', 'node', 'layer', 1)
        self.moulin_input = project3d(md, 'vector', self.moulin_input, 'type', 'node', 'layer', 1)
        self.neumannflux = project3d(md, 'vector', self.neumannflux, 'type', 'node', 'layer', 1)
        return self
    # }}}

    def setdefaultparameters(self):  # {{{
        # Sheet parameters
        self.pressure_melt_coefficient = 7.5e-8  #K / Pa (See table 1 in Werder et al. 2013)
        self.cavity_spacing = 2.  #m
        self.sheet_alpha = 5.0/4.0
        self.sheet_beta = 3.0/2.0
        self.omega = 1./2000. 
        self.creep_open_flag = 1

        # Channel parameters
        self.ischannels = False
        self.channel_conductivity = 5.e-2  #Dow's default, Table uses 0.1
        self.channel_sheet_width = 2.  #m
        self.channel_alpha = 5.0/4.0
        self.channel_beta = 3.0/2.0

        # Other
        self.englacial_void_ratio = 1.e-5  #Dow's default, Table from Werder et al. uses 1e-3
        self.requested_outputs = ['default']
        self.melt_flag = 0
        self.istransition = 0  #by default use turbulent physics

        return self
    # }}}

    def checkconsistency(self, md, solution, analyses):  # {{{
        # Early return
        if 'HydrologyGladsAnalysis' not in analyses:
            return md

        # Sheet
        md = checkfield(md, 'fieldname', 'hydrology.pressure_melt_coefficient', 'numel', [1], '>=', 0)
        md = checkfield(md, 'fieldname', 'hydrology.sheet_conductivity', 'size', [md.mesh.numberofvertices], '>', 0, 'np.nan', 1, 'Inf', 1)
        md = checkfield(md, 'fieldname', 'hydrology.cavity_spacing', 'numel', [1], '>', 0)
        md = checkfield(md, 'fieldname', 'hydrology.bump_height', 'size', [md.mesh.numberofvertices], '>=', 0, 'np.nan', 1, 'Inf', 1)
        md = checkfield(md,'fieldname','hydrology.omega', 'numel', [1], '>=', 0) 
        md = checkfield(md,'fieldname','hydrology.sheet_alpha', 'numel', [1], '>', 0) 
        md = checkfield(md,'fieldname','hydrology.sheet_beta', 'numel', [1], '>', 0) 
        md = checkfield(md,'fieldname','hydrology.rheology_B_base', 'size', [md.mesh.numberofvertices], '>=', 0, 'np.nan', 1, 'Inf', 1)
        md = checkfield(md, 'fieldname', 'hydrology.isincludesheetthickness', 'numel', [1], 'values', [0, 1])
        md = checkfield(md, 'fieldname', 'hydrology.creep_open_flag', 'numel', [1], 'values', [0, 1])
        md = checkfield(md,'fieldname','hydrology.rheology_B_base', 'size', [md.mesh.numberofvertices], '>=', 0, 'np.nan', 1, 'Inf', 1)

        # Channels
        md = checkfield(md, 'fieldname', 'hydrology.ischannels', 'numel', [1], 'values', [0, 1])
        md = checkfield(md, 'fieldname', 'hydrology.channel_conductivity', 'size', [md.mesh.numberofvertices], '>', 0)
        md = checkfield(md, 'fieldname', 'hydrology.channel_sheet_width', 'numel', [1], '>=', 0)
        md = checkfield(md,'fieldname','hydrology.channel_alpha', 'numel', [1], '>', 0) 
        md = checkfield(md,'fieldname','hydrology.channel_beta', 'numel', [1], '>', 0) 

        # Other
        md = checkfield(md, 'fieldname', 'hydrology.spcphi', 'Inf', 1, 'timeseries', 1)
        md = checkfield(md, 'fieldname', 'hydrology.englacial_void_ratio', 'numel', [1], '>=', 0)
        md = checkfield(md, 'fieldname', 'hydrology.moulin_input', '>=', 0, 'timeseries', 1, 'NaN', 1, 'Inf', 1)
        md = checkfield(md, 'fieldname', 'hydrology.neumannflux', 'timeseries', 1, 'NaN', 1, 'Inf', 1)
        md = checkfield(md, 'fieldname', 'hydrology.requested_outputs', 'stringrow', 1)
        md = checkfield(md, 'fieldname', 'hydrology.melt_flag', 'numel', [1], 'values', [0, 1])
        md = checkfield(md, 'fieldname', 'hydrology.istransition', 'numel', [1], 'values', [0, 1])
        if self.melt_flag == 1 or self.melt_flag == 2:
            md = checkfield(md, 'fieldname', 'basalforcings.groundedice_melting_rate', 'NaN', 1, 'Inf', 1, 'timeseries', 1)
    # }}}

    def marshall(self, prefix, md, fid):  # {{{
        yts = md.constants.yts
        # Marshall model code first
        WriteData(fid, prefix, 'name', 'md.hydrology.model', 'data', 5, 'format', 'Integer')

       # Sheet
        WriteData(fid, prefix, 'object', self, 'class', 'hydrology', 'fieldname', 'pressure_melt_coefficient', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'class', 'hydrology', 'fieldname', 'sheet_conductivity', 'format', 'DoubleMat', 'mattype', 1)
        WriteData(fid, prefix, 'object', self, 'class', 'hydrology', 'fieldname', 'cavity_spacing', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'class', 'hydrology', 'fieldname', 'bump_height', 'format', 'DoubleMat', 'mattype', 1)
        WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','omega','format','Double') 
        WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','sheet_alpha','format','Double') 
        WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','sheet_beta','format','Double') 
        WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','rheology_B_base','format','DoubleMat', 'mattype', 1)
        WriteData(fid, prefix, 'object', self, 'class', 'hydrology', 'fieldname', 'isincludesheetthickness', 'format', 'Boolean')
        WriteData(fid, prefix, 'object', self, 'class', 'hydrology', 'fieldname', 'creep_open_flag', 'format', 'Boolean')
        WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','rheology_B_base','format','DoubleMat', 'mattype', 1);

        # Channels
        WriteData(fid, prefix, 'object', self, 'class', 'hydrology', 'fieldname', 'ischannels', 'format', 'Boolean')
        WriteData(fid, prefix, 'object', self, 'class', 'hydrology', 'fieldname', 'channel_conductivity', 'format', 'DoubleMat', 'mattype', 1)
        WriteData(fid, prefix, 'object', self, 'class', 'hydrology', 'fieldname', 'channel_sheet_width', 'format', 'Double')
        WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','channel_alpha','format','Double') 
        WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','channel_beta','format','Double') 

        # Others
        WriteData(fid, prefix, 'object', self, 'class', 'hydrology', 'fieldname', 'spcphi', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', yts)
        WriteData(fid, prefix, 'object', self, 'class', 'hydrology', 'fieldname', 'neumannflux', 'format', 'DoubleMat', 'mattype', 2, 'timeserieslength', md.mesh.numberofelements + 1, 'yts', yts)
        WriteData(fid, prefix, 'object', self, 'class', 'hydrology', 'fieldname', 'moulin_input', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', md.constants.yts)
        WriteData(fid, prefix, 'object', self, 'class', 'hydrology', 'fieldname', 'englacial_void_ratio', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'class', 'hydrology', 'fieldname', 'melt_flag', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'class', 'hydrology', 'fieldname', 'istransition', 'format', 'Boolean')

        outputs = self.requested_outputs
        indices = [i for i, x in enumerate(outputs) if x == 'default']
        if len(indices) > 0:
            outputscopy = outputs[0:max(0, indices[0] - 1)] + self.defaultoutputs(md) + outputs[indices[0] + 1:]
            outputs = outputscopy
        WriteData(fid, prefix, 'data', outputs, 'name', 'md.hydrology.requested_outputs', 'format', 'StringArray')
    # }}}
