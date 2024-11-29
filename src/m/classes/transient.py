from checkfield import checkfield
from fielddisplay import fielddisplay
from WriteData import WriteData


class transient(object):
    """TRANSIENT class definition

    Usage:
        transient = transient()
    """

    def __init__(self, *args):  # {{{
        self.isage = 0
        self.issmb = 0
        self.ismasstransport = 0
        self.ismmemasstransport = 0
        self.isoceantransport = 0
        self.isstressbalance = 0
        self.isthermal = 0
        self.isgroundingline = 0
        self.isesa = 0
        self.isdamageevolution = 0
        self.ismovingfront = 0
        self.ishydrology = 0
        self.isdebris = 0
        self.issampling = 0
        self.isslc = 0
        self.amr_frequency = 0
        self.isoceancoupling = 0
        self.requested_outputs = []

        if len(args) == 0:
            self.setdefaultparameters()
        else:
            raise Exception('constructor not supported')
    # }}}

    def __repr__(self):  # {{{
        s = '   transient solution parameters:\n'
        s += '{}\n'.format(fielddisplay(self, 'isage', 'indicates if age model is requested in the transient'))
        s += '{}\n'.format(fielddisplay(self, 'issmb', 'indicates if a surface mass balance solution is used in the transient'))
        s += '{}\n'.format(fielddisplay(self, 'ismasstransport', 'indicates if a masstransport solution is used in the transient'))
        s += '{}\n'.format(fielddisplay(self, 'ismmemasstransport', 'indicates whether an MME masstransport solution is used in the transient'))
        s += '{}\n'.format(fielddisplay(self, 'isoceantransport', 'indicates whether an ocean masstransport solution is used in the transient'))
        s += '{}\n'.format(fielddisplay(self, 'isstressbalance', 'indicates if a stressbalance solution is used in the transient'))
        s += '{}\n'.format(fielddisplay(self, 'isthermal', 'indicates if a thermal solution is used in the transient'))
        s += '{}\n'.format(fielddisplay(self, 'isgroundingline', 'indicates if a groundingline migration is used in the transient'))
        s += '{}\n'.format(fielddisplay(self, 'isesa', 'indicates whether an elastic adjustment model is used in the transient'))
        s += '{}\n'.format(fielddisplay(self, 'isdamageevolution', 'indicates whether damage evolution is used in the transient'))
        s += '{}\n'.format(fielddisplay(self, 'ismovingfront', 'indicates whether a moving front capability is used in the transient'))
        s += '{}\n'.format(fielddisplay(self, 'ishydrology', 'indicates whether an hydrology model is used'))
        s += '{}\n'.format(fielddisplay(self, 'isdebris', 'indicates whether a debris model is used'))
        s += '{}\n'.format(fielddisplay(self, 'issampling', 'indicates whether sampling is used in the transient'))
        s += '{}\n'.format(fielddisplay(self, 'isslc', 'indicates if a sea level change solution is used in the transient'))
        s += '{}\n'.format(fielddisplay(self, 'isoceancoupling', 'indicates whether coupling with an ocean model is used in the transient (1 for cartesian coordinates, 2 for lat/long coordinates'))
        s += '{}\n'.format(fielddisplay(self, 'amr_frequency', 'frequency at which mesh is refined in simulations with multiple time_steps'))
        s += '{}\n'.format(fielddisplay(self, 'requested_outputs', 'list of additional outputs requested'))
        return s
    # }}}

    def defaultoutputs(self, md):  # {{{
        return []
    # }}}

    def deactivateall(self):  #{{{
        self.isage = 0
        self.issmb = 0
        self.ismasstransport = 0
        self.ismmemasstransport = 0
        self.isoceantransport = 0
        self.isstressbalance = 0
        self.isthermal = 0
        self.isgroundingline = 0
        self.isesa = 0
        self.isdamageevolution = 0
        self.ismovingfront = 0
        self.ishydrology = 0
        self.isdebris = 0
        self.issampling = 0
        self.isslc = 0
        self.isoceancoupling = 0
        self.amr_frequency = 0

        # Default output
        self.requested_outputs = []
        return self
    # }}}

    def setdefaultparameters(self):  # {{{
        # Full analysis: Stressbalance, Masstransport and Thermal but no 
        # groundingline migration for now
        self.isage = 0
        self.issmb = 1
        self.ismasstransport = 1
        self.ismmemasstransport = 0
        self.isoceantransport = 0
        self.isstressbalance = 1
        self.isthermal = 1
        self.isgroundingline = 0
        self.isesa = 0
        self.isdamageevolution = 0
        self.ismovingfront = 0
        self.ishydrology = 0
        self.isdebris = 0
        self.issampling = 0
        self.isslc = 0
        self.isoceancoupling = 0
        self.amr_frequency = 0

        # Default output
        self.requested_outputs = ['default']
        return self
    # }}}

    def checkconsistency(self, md, solution, analyses):  # {{{
        #Early return
        if not solution == 'TransientSolution':
            return md

        md = checkfield(md, 'fieldname', 'transient.isage', 'numel', [1], 'values', [0, 1])
        md = checkfield(md, 'fieldname', 'transient.issmb', 'numel', [1], 'values', [0, 1])
        md = checkfield(md, 'fieldname', 'transient.ismasstransport', 'numel', [1], 'values', [0, 1])
        md = checkfield(md, 'fieldname', 'transient.ismmemasstransport', 'numel', [1], 'values', [0, 1])
        md = checkfield(md, 'fieldname', 'transient.isoceantransport', 'numel', [1], 'values', [0, 1])
        md = checkfield(md, 'fieldname', 'transient.isstressbalance', 'numel', [1], 'values', [0, 1])
        md = checkfield(md, 'fieldname', 'transient.isthermal', 'numel', [1], 'values', [0, 1])
        md = checkfield(md, 'fieldname', 'transient.isgroundingline', 'numel', [1], 'values', [0, 1])
        md = checkfield(md, 'fieldname', 'transient.isesa', 'numel', [1], 'values', [0, 1])
        md = checkfield(md, 'fieldname', 'transient.isdamageevolution', 'numel', [1], 'values', [0, 1])
        md = checkfield(md, 'fieldname', 'transient.ishydrology', 'numel', [1], 'values', [0, 1])
        md = checkfield(md, 'fieldname', 'transient.isdebris', 'numel', [1], 'values', [0, 1])
        md = checkfield(md, 'fieldname', 'transient.issampling', 'numel', [1], 'values', [0, 1])
        md = checkfield(md, 'fieldname', 'transient.ismovingfront', 'numel', [1], 'values', [0, 1])
        md = checkfield(md, 'fieldname', 'transient.isslc', 'numel', [1], 'values', [0, 1])
        md = checkfield(md, 'fieldname', 'transient.isoceancoupling', 'numel', [1], 'values', [0, 1, 2])
        md = checkfield(md, 'fieldname', 'transient.amr_frequency', 'numel', [1], '>=', 0, 'NaN', 1, 'Inf', 1)

        if solution != 'TransientSolution' and md.transient.iscoupling:
            md.checkmessage("Coupling with ocean can only be done in transient simulations!")
        if md.transient.isdamageevolution and not hasattr(md.materials, 'matdamageice'):
            md.checkmessage("requesting damage evolution but md.materials is not of class matdamageice")
        return md
    # }}}

    def marshall(self, prefix, md, fid):  # {{{
        WriteData(fid, prefix, 'object', self, 'fieldname', 'isage', 'format', 'Boolean')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'issmb', 'format', 'Boolean')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'ismasstransport', 'format', 'Boolean')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'ismmemasstransport', 'format', 'Boolean')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'isoceantransport', 'format', 'Boolean')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'isstressbalance', 'format', 'Boolean')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'isthermal', 'format', 'Boolean')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'isgroundingline', 'format', 'Boolean')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'isesa', 'format', 'Boolean')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'isdamageevolution', 'format', 'Boolean')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'ismovingfront', 'format', 'Boolean')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'ishydrology', 'format', 'Boolean')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'isdebris', 'format', 'Boolean')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'issampling', 'format', 'Boolean')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'isslc', 'format', 'Boolean')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'isoceancoupling', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'amr_frequency', 'format', 'Integer')

        # Process requested outputs
        outputs = self.requested_outputs
        indices = [i for i, x in enumerate(outputs) if x == 'default']
        if len(indices) > 0:
            outputscopy = outputs[0:max(0, indices[0] - 1)] + self.defaultoutputs(md) + outputs[indices[0] + 1:]
            outputs = outputscopy
        WriteData(fid, prefix, 'data', outputs, 'name', 'md.transient.requested_outputs', 'format', 'StringArray')
    # }}}
