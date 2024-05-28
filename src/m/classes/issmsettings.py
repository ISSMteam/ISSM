from checkfield import checkfield
from fielddisplay import fielddisplay
from WriteData import WriteData


class issmsettings(object):
    """ISSMSETTINGS class definition

    Usage:
        issmsettings = issmsettings()
    """

    def __init__(self):  # {{{
        self.results_on_nodes = []
        self.io_gather = 0
        self.lowmem = 0
        self.output_frequency = 0
        self.coupling_frequency = 0
        self.checkpoint_frequency = 0
        self.waitonlock = 0
        self.solver_residue_threshold = 0

        # Set defaults
        self.setdefaultparameters()
    # }}}

    def __repr__(self):  # {{{
        s = "   general issmsettings parameters:\n"
        s += '{}\n'.format(fielddisplay(self, "results_on_nodes", "list of output for which results will be output for all the nodes of each element, Use 'all' for all output on nodes."))
        s += '{}\n'.format(fielddisplay(self, "io_gather", "I / O gathering strategy for result outputs (default 1)"))
        s += '{}\n'.format(fielddisplay(self, "lowmem", "is the memory limited ? (0 or 1)"))
        s += '{}\n'.format(fielddisplay(self, "output_frequency", "number of time steps between two saves (e.g., 5 means that results are only saved every 5 time steps)"))
        s += '{}\n'.format(fielddisplay(self, "sb_coupling_frequency", "frequency at which StressBalance solver is coupled (default 1)"))
        s += '{}\n'.format(fielddisplay(self, "checkpoint_frequency", "frequency at which the runs are being recorded, allowing for a restart"))
        s += '{}\n'.format(fielddisplay(self, "waitonlock", "maximum number of minutes to wait for batch results, or return 0"))
        s += '{}\n'.format(fielddisplay(self, "solver_residue_threshold", "throw an error if solver residue exceeds this value (NaN to deactivate)"))
        return s
    # }}}

    def setdefaultparameters(self):  # {{{
        # Are we short in memory? (0 faster but requires more memory)
        self.lowmem = 0
        # I/O:
        self.io_gather = 1
        # Results frequency by default every step
        self.output_frequency = 1
        # Coupling frequency of the stress balance solver by default every step
        self.sb_coupling_frequency = 1
        # Checkpoints frequency, by default never:
        self.checkpoint_frequency = 0
        # This option can be activated to load automatically the results onto 
        # the model after a parallel run by waiting for the lock file N minutes 
        # that is generated once the solution has converged
        # Set to 0 to deactivate
        self.waitonlock = pow(2, 31) - 1
        # Throw an error if solver residue exceeds this value
        self.solver_residue_threshold = 1e-6

        return self
    # }}}

    def checkconsistency(self, md, solution, analyses):  # {{{
        md = checkfield(md, 'fieldname', 'settings.results_on_nodes', 'stringrow', 1)
        md = checkfield(md, 'fieldname', 'settings.io_gather', 'numel', [1], 'values', [0, 1])
        md = checkfield(md, 'fieldname', 'settings.lowmem', 'numel', [1], 'values', [0, 1])
        md = checkfield(md, 'fieldname', 'settings.output_frequency', 'numel', [1], '>=', 1)
        md = checkfield(md, 'fieldname', 'settings.sb_coupling_frequency', 'numel', [1], '>=', 1)
        md = checkfield(md, 'fieldname', 'settings.checkpoint_frequency', 'numel', [1], '>=', 0)
        md = checkfield(md, 'fieldname', 'settings.waitonlock', 'numel', [1])
        md = checkfield(md, 'fieldname', 'settings.solver_residue_threshold', 'numel', [1], '>', 0)

        return md
    # }}}

    def marshall(self, prefix, md, fid):  # {{{
        WriteData(fid, prefix, 'data', self.results_on_nodes, 'name', 'md.settings.results_on_nodes', 'format', 'StringArray')
        WriteData(fid, prefix, 'object', self, 'class', 'settings', 'fieldname', 'io_gather', 'format', 'Boolean')
        WriteData(fid, prefix, 'object', self, 'class', 'settings', 'fieldname', 'lowmem', 'format', 'Boolean')
        WriteData(fid, prefix, 'object', self, 'class', 'settings', 'fieldname', 'output_frequency', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'class', 'settings', 'fieldname', 'sb_coupling_frequency', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'class', 'settings', 'fieldname', 'checkpoint_frequency', 'format', 'Integer')

        if self.waitonlock > 0:
            WriteData(fid, prefix, 'name', 'md.settings.waitonlock', 'data', True, 'format', 'Boolean')
        else:
            WriteData(fid, prefix, 'name', 'md.settings.waitonlock', 'data', False, 'format', 'Boolean')

        WriteData(fid, prefix, 'object', self, 'fieldname', 'solver_residue_threshold', 'format', 'Double')
    # }}}
