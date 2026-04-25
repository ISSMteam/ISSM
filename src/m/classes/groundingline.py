import numpy as np
from checkfield import checkfield
from fielddisplay import fielddisplay
import MatlabFuncs as m
from WriteData import WriteData


class groundingline(object):
    """
    GROUNDINGLINE class definition

       Usage:
          groundingline = groundingline()
    """

    def __init__(self):  # {{{
        self.migration = ''
        self.friction_interpolation = ''
        self.melt_interpolation = ''
        self.intrusion_distance = 0
        self.requested_outptuts = []

        # Set defaults
        self.setdefaultparameters()

    # }}}

    def __repr__(self):  # {{{
        s = '   grounding line migration parameters:\n'
        s += '{}\n'.format(fielddisplay(self, 'migration', 'type of grounding line migration: \'SoftMigration\', \'SubelementMigration\', \'AggressiveMigration\', \'Contact\', \'None\''))
        s += '{}\n'.format(fielddisplay(self, 'migration', 'type of friction interpolation on partially floating elements: ''SubelementFriction1'', ''SubelementFriction2'', ''NoFrictionOnPartiallyFloating'''))
        s += '{}\n'.format(fielddisplay(self, 'migration', 'type of melt interpolation on partially floating elements: \'NoMeltOnPartiallyFloating\', \'FullMeltOnPartiallyFloating\', \'SubelementMelt1\', \'SubelementMelt2\', \'IntrusionMelt\''))
        s += '{}\n'.format(fielddisplay(self, 'requested_outputs', 'additional outputs requested'))
        return s
    # }}}

    def defaultoutputs(self, md):  # {{{
        return ['Surface', 'Base','MaskOceanLevelset']

    # }}}

    def setdefaultparameters(self):  # {{{
        # Type of migration
        self.migration = 'SubelementMigration'
        self.friction_interpolation = 'SubelementFriction1'
        self.melt_interpolation = 'NoMeltOnPartiallyFloating'
        self.intrusion_distance =  0
        # Default output
        self.requested_outputs = ['default']

        return self
    # }}}

    def checkconsistency(self, md, solution, analyses):  # {{{
        md = checkfield(md, 'fieldname', 'groundingline.migration', 'values', ['None', 'SubelementMigration', 'AggressiveMigration', 'SoftMigration', 'Contact', 'GroundingOnly'])
        md = checkfield(md, 'fieldname', 'groundingline.friction_interpolation', 'values', ['SubelementFriction1', 'SubelementFriction2', 'NoFrictionOnPartiallyFloating'])
        md = checkfield(md, 'fieldname', 'groundingline.melt_interpolation', 'values', ['NoMeltOnPartiallyFloating', 'FullMeltOnPartiallyFloating', 'SubelementMelt1', 'SubelementMelt2', 'IntrusionMelt'])
        md = checkfield(md, 'fieldname', 'groundingline.intrusion_distance', 'NaN', 1, 'Inf', 1, '>=', 0)
        md = checkfield(md, 'fieldname', 'groundingline.requested_outputs', 'stringrow', 1)

        if(not m.strcmp(self.migration, 'None') and md.transient.isgroundingline and solution == 'TransientSolution'):
            if np.any(np.isnan(md.geometry.bed)):
                md.checkmessage("requesting grounding line migration, but bathymetry is absent!")
            pos = np.nonzero(md.mask.ocean_levelset > 0.)[0]
            if any(np.abs(md.geometry.base[pos] - md.geometry.bed[pos]) > pow(10, -10)):
                md.checkmessage("base not equal to bed on grounded ice!")
            if any(md.geometry.bed - md.geometry.base > pow(10, -9)):
                md.checkmessage("bed superior to base on floating ice!")

        return md
    # }}}

    def marshall(self, prefix, md, fid):  # {{{
        WriteData(fid, prefix, 'data', self.migration, 'name', 'md.groundingline.migration', 'format', 'String')
        WriteData(fid, prefix, 'data', self.friction_interpolation, 'name', 'md.groundingline.friction_interpolation', 'format', 'String')
        WriteData(fid, prefix, 'data', self.melt_interpolation, 'name', 'md.groundingline.melt_interpolation', 'format', 'String')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'intrusion_distance', 'format', 'DoubleMat', 'mattype', 1)
        
        # Process requested outputs
        outputs = self.requested_outputs
        indices = [i for i, x in enumerate(outputs) if x == 'default']
        if len(indices) > 0:
            outputscopy = outputs[0:max(0, indices[0] - 1)] + self.defaultoutputs(md) + outputs[indices[0] + 1:]
            outputs = outputscopy
        WriteData(fid, prefix, 'data', outputs, 'name', 'md.groundingline.requested_outputs', 'format', 'StringArray')
    # }}}
