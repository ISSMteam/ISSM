import numpy as np

from checkfield import checkfield
from fielddisplay import fielddisplay
import MatlabFuncs as m
from project3d import project3d
from WriteData import WriteData


class flowequation(object):
    """FLOWEQUATION class definition

    Usage:
        flowequation = flowequation()
    """

    def __init__(self):  # {{{
        self.isSIA = 0
        self.isSSA = 0
        self.isL1L2 = 0
        self.isMOLHO = 0
        self.isHO = 0
        self.isFS = 0
        self.isNitscheBC = 0
        self.FSNitscheGamma = 1e6
        self.fe_SSA = ''
        self.fe_HO = ''
        self.fe_FS = ''
        self.augmented_lagrangian_r = 1
        self.augmented_lagrangian_rhop = 1
        self.augmented_lagrangian_rlambda = 1
        self.augmented_lagrangian_rholambda = 1
        self.XTH_theta = 0
        self.vertex_equation = np.nan
        self.element_equation = np.nan
        self.borderSSA = np.nan
        self.borderHO = np.nan
        self.borderFS = np.nan
        
        self.setdefaultparameters()
    # }}}

    def __repr__(self):  # {{{
        s = '   flow equation parameters:\n'
        s += '{}\n'.format(fielddisplay(self, 'isSIA', "is the Shallow Ice Approximation (SIA) used?"))
        s += '{}\n'.format(fielddisplay(self, 'isSSA', "is the Shelfy-Stream Approximation (SSA) used?"))
        s += '{}\n'.format(fielddisplay(self, 'isL1L2', "are L1L2 equations used?"))
        s += '{}\n'.format(fielddisplay(self, 'isMOLHO', "are MOno-layer Higher-Order (MOLHO) equations used?"))
        s += '{}\n'.format(fielddisplay(self, 'isHO', "is the Higher-Order (HO) approximation used?"))
        s += '{}\n'.format(fielddisplay(self, 'isFS', "are the Full-FS (FS) equations used?"))
        s += '{}\n'.format(fielddisplay(self, 'isNitscheBC', "is weakly imposed condition used?"))
        s += '{}\n'.format(fielddisplay(self, 'FSNitscheGamma', "Gamma value for the Nitsche term (default: 1e6)"))
        s += '{}\n'.format(fielddisplay(self, 'fe_SSA', "Finite Element for SSA: 'P1', 'P1bubble' 'P1bubblecondensed' 'P2'"))
        s += '{}\n'.format(fielddisplay(self, 'fe_HO', "Finite Element for HO:  'P1', 'P1bubble', 'P1bubblecondensed', 'P1xP2', 'P2xP1', 'P2', 'P2bubble', 'P1xP3', 'P2xP4'"))
        s += '{}\n'.format(fielddisplay(self, 'fe_FS', "Finite Element for FS:  'P1P1' (debugging only) 'P1P1GLS' 'MINIcondensed' 'MINI' 'TaylorHood' 'LATaylorHood' 'XTaylorHood'"))
        s += '{}\n'.format(fielddisplay(self, 'vertex_equation', "flow equation for each vertex"))
        s += '{}\n'.format(fielddisplay(self, 'element_equation', "flow equation for each element"))
        s += '{}\n'.format(fielddisplay(self, 'borderSSA', "vertices on SSA's border (for tiling)"))
        s += '{}\n'.format(fielddisplay(self, 'borderHO', "vertices on HO's border (for tiling)"))
        s += '{}\n'.format(fielddisplay(self, 'borderFS', "vertices on FS' border (for tiling)"))
        return s
    # }}}

    def setdefaultparameters(self):  # {{{
        # P1 for SSA
        self.fe_SSA = 'P1'

        # P1 for HO
        self.fe_HO = 'P1'

        # MINI condensed element for FS by default
        self.fe_FS = 'MINIcondensed'
        return self
    # }}}

    def extrude(self, md):  # {{{
        self.element_equation = project3d(md, 'vector', self.element_equation, 'type', 'element')
        self.vertex_equation = project3d(md, 'vector', self.vertex_equation, 'type', 'node')
        self.borderSSA = project3d(md, 'vector', self.borderSSA, 'type', 'node')
        self.borderHO = project3d(md, 'vector', self.borderHO, 'type', 'node')
        self.borderFS = project3d(md, 'vector', self.borderFS, 'type', 'node')
        return self
    # }}}

    def checkconsistency(self, md, solution, analyses):  # {{{
        # Early return
        if ('StressbalanceAnalysis' not in analyses and 'StressbalanceSIAAnalysis' not in analyses) or (solution == 'TransientSolution' and not md.transient.isstressbalance):
            return md
        md = checkfield(md, 'fieldname', 'flowequation.isSIA', 'numel', [1], 'values', [0, 1])
        md = checkfield(md, 'fieldname', 'flowequation.isSSA', 'numel', [1], 'values', [0, 1])
        md = checkfield(md, 'fieldname', 'flowequation.isL1L2', 'numel', [1], 'values', [0, 1])
        md = checkfield(md, 'fieldname', 'flowequation.isMOLHO', 'numel', [1], 'values', [0, 1])
        md = checkfield(md, 'fieldname', 'flowequation.isHO', 'numel', [1], 'values', [0, 1])
        md = checkfield(md, 'fieldname', 'flowequation.isFS', 'numel', [1], 'values', [0, 1])
        md = checkfield(md, 'fieldname', 'flowequation.isNitscheBC', 'numel', [1], 'values', [0, 1])
        md = checkfield(md, 'fieldname', 'flowequation.FSNitscheGamma', 'numel', [1], '>=', 0.)
        md = checkfield(md, 'fieldname', 'flowequation.fe_SSA', 'values', ['P1', 'P1bubble', 'P1bubblecondensed', 'P2', 'P2bubble'])
        md = checkfield(md, 'fieldname', 'flowequation.fe_HO', 'values', ['P1', 'P1bubble', 'P1bubblecondensed', 'P1xP2', 'P2xP1', 'P2', 'P2bubble', 'P1xP3', 'P2xP4'])
        md = checkfield(md, 'fieldname', 'flowequation.fe_FS', 'values', ['P1P1', 'P1P1GLS', 'MINIcondensed', 'MINI', 'TaylorHood', 'LATaylorHood', 'XTaylorHood', 'OneLayerP4z', 'CrouzeixRaviart', 'LACrouzeixRaviart'])
        md = checkfield(md, 'fieldname', 'flowequation.borderSSA', 'size', [md.mesh.numberofvertices], 'values', [0, 1])
        md = checkfield(md, 'fieldname', 'flowequation.borderHO', 'size', [md.mesh.numberofvertices], 'values', [0, 1])
        md = checkfield(md, 'fieldname', 'flowequation.borderFS', 'size', [md.mesh.numberofvertices], 'values', [0, 1])
        md = checkfield(md, 'fieldname', 'flowequation.augmented_lagrangian_r', 'numel', [1], '>', 0.)
        md = checkfield(md, 'fieldname', 'flowequation.augmented_lagrangian_rhop', 'numel', [1], '>', 0.)
        md = checkfield(md, 'fieldname', 'flowequation.augmented_lagrangian_rlambda', 'numel', [1], '>', 0.)
        md = checkfield(md, 'fieldname', 'flowequation.augmented_lagrangian_rholambda', 'numel', [1], '>', 0.)
        md = checkfield(md, 'fieldname', 'flowequation.XTH_theta', 'numel', [1], '>=', 0., '<', .5)
        if m.strcmp(md.mesh.domaintype(), '2Dhorizontal'):
            md = checkfield(md, 'fieldname', 'flowequation.vertex_equation', 'size', [md.mesh.numberofvertices], 'values', [1, 2, 4])
            md = checkfield(md, 'fieldname', 'flowequation.element_equation', 'size', [md.mesh.numberofelements], 'values', [1, 2, 4])
        elif m.strcmp(md.mesh.domaintype(), '3Dsurface'):
            md = checkfield(md, 'fieldname', 'flowequation.vertex_equation', 'size', [md.mesh.numberofvertices], 'values', np.arange(1, 2 + 1))
            md = checkfield(md, 'fieldname', 'flowequation.element_equation', 'size', [md.mesh.numberofelements], 'values', np.arange(1, 2 + 1))
        elif m.strcmp(md.mesh.domaintype(), '2Dvertical'):
            md = checkfield(md, 'fieldname', 'flowequation.vertex_equation', 'size', [md.mesh.numberofvertices], 'values', [2, 5, 6])
            md = checkfield(md, 'fieldname', 'flowequation.element_equation', 'size', [md.mesh.numberofelements], 'values', [2, 5, 6])
        elif m.strcmp(md.mesh.domaintype(), '3D'):
            md = checkfield(md, 'fieldname', 'flowequation.vertex_equation', 'size', [md.mesh.numberofvertices], 'values', np.arange(0, 9 + 1))
            md = checkfield(md, 'fieldname', 'flowequation.element_equation', 'size', [md.mesh.numberofelements], 'values', np.arange(0, 9 + 1))
        else:
            raise RuntimeError('Case not supported yet')

        if not (self.isSIA or self.isSSA or self.isL1L2 or self.isMOLHO or self.isHO or self.isFS):
            md.checkmessage("no element types set for this model")
        if 'StressbalanceSIAAnalysis' in analyses:
            if any(self.element_equation == 1):
                if np.any(np.logical_and(self.vertex_equation, md.mask.ocean_levelset)):
                    print("\n !!! Warning: SIA's model is not consistent on ice shelves !!!\n")
        return md
    # }}}

    def marshall(self, prefix, md, fid):  # {{{
        WriteData(fid, prefix, 'object', self, 'fieldname', 'isSIA', 'format', 'Boolean')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'isSSA', 'format', 'Boolean')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'isL1L2', 'format', 'Boolean')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'isMOLHO', 'format', 'Boolean')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'isHO', 'format', 'Boolean')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'isFS', 'format', 'Boolean')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'isNitscheBC', 'format', 'Boolean')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'FSNitscheGamma', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'fe_SSA', 'data', self.fe_SSA, 'format', 'String')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'fe_HO', 'data', self.fe_HO, 'format', 'String')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'fe_FS', 'data', self.fe_FS, 'format', 'String')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'augmented_lagrangian_r', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'augmented_lagrangian_rhop', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'augmented_lagrangian_rlambda', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'augmented_lagrangian_rholambda', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'XTH_theta', 'data', self.XTH_theta, 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'borderSSA', 'format', 'DoubleMat', 'mattype', 1)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'borderHO', 'format', 'DoubleMat', 'mattype', 1)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'borderFS', 'format', 'DoubleMat', 'mattype', 1)
        # Convert approximations to enums
        WriteData(fid, prefix, 'data', self.vertex_equation, 'name', 'md.flowequation.vertex_equation', 'format', 'DoubleMat', 'mattype', 1)
        WriteData(fid, prefix, 'data', self.element_equation, 'name', 'md.flowequation.element_equation', 'format', 'DoubleMat', 'mattype', 2)
    # }}}
