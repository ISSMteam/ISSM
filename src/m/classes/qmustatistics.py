import numpy as np

from helpers import *
from WriteData import WriteData


class qmustatistics(object):
    """QMUSTATISTICS class definition

    Usage:
        stats = qmustatistics()
    """

    def __init__(self, *args):  #{{{
        self.nfiles_per_directory = 5  # Number of files per output directory
        self.ndirectories  = 50  # Number of output directories; should be < numcpus

        self.method = [{}]
        self.method[0]['name'] = 'None'
        self.method[0]['fields'] = []
        self.method[0]['steps'] = []
        self.method[0]['nbins'] = np.nan
        self.method[0]['indices'] = []

        # name : name of method, one of 'None', 'Histogram', 'SampleSeries', or 'MeanVariance'
        # fields : fields for the statistics being requested, ex: 'Sealevel', 'BslrIce', 'BsrlHydro'
        # steps : time steps at which each field statistic is computed, ex: [1, 2, 5, 20] or [range(1:100)]
        # nbins : number of bins for 'Histogram' statistics
        # indices : vertex indices at which to retrieve samples

        nargs = len(args)
        if nargs == 0:
            # Create a default object
            self.setdefaultparameters()
        elif nargs == 1:
            # NOTE: The following has not been tested. Remove this note when it has
            inputstruct = args[0]
            list1 = properties('qmustatistics')
            list2 = fieldnames(inputstruct)
            for i in range(len(list1)):
                fieldname = list1[i]
                if fieldname in list2:
                    setattr(self, fieldname, getattr(inputstruct, fieldname))
        else:
            raise Exception('constructor not supported')
    # }}}

    def __repr__(self):  # {{{
        s = 'qmustatistics: post-Dakota run processing of QMU statistics:\n'
        if self.method[0]['name'] == 'None':
            return s
        s += '{}\n'.format(fielddisplay(self, 'nfiles_per_directory', 'number of files per output directory'))
        s += '{}\n'.format(fielddisplay(self, 'ndirectories', 'number of output directories; should be < numcpus'))
        for i in range(len(self.method)):
            s += '{}\n'.format('   method # {}'.format(i))
            s += '{}\n'.format(self.method[i])
        return s
    # }}}

    def setdefaultparameters(self):  # {{{
        self.method[0]['name'] = 'None'
        self.nfiles_per_directory = 5  # Number of files per output directory
        self.ndirectories = 50  # Number of output directories; should be < numcpus
        return self
    # }}}

    @staticmethod
    def checkconsistency(self, md, solution, analyses):  # {{{
        # Early return
        if md.qmu.isdakota == 0:
            return
        if self.method[0]['name'] == 'None':
            return

        # Checks
        md = checkfield(md, 'fieldname', 'qmu.statistics.nfiles_per_directory', '>=', 1)
        if self.ndirectories > md.cluster.np:
            raise Exception('qmustatistics consistency check: number of cluster CPUs should be > number of output directories')
        if self.ndirectories * self.nfiles_per_directory != md.qmu.method.params.sample:
            raise Exception('qmustatistics consistency check: number of directories x number of files per directory should be == to number of samples requested!')
        for i in range(len(self.method)):
            m = self.method[i]
            if m.name == 'Histogram':
                md = checkfield(md, 'fieldname', 'qmu.statistics.method[{}].nbins'.format(i), '>=', 1, '<=', md.qmu.method.params.samples)
            for f in range(len(m['fields'])):
                if not isinstance(m['fields'][f], str):
                    raise Exception('qmustatistics consistency check error: qmu.statistics.method[{}][\'fields\'][{}] is not a string!'.format(i, f))
            for s in range(len(m['steps'])):
                if m['steps'][s] <= 0:
                    raise Exception('qmustatistics consistency check error: qmu.statistics.method[{}][\'steps\'][{}] should be > 0!'.format(i, s))
                if m['steps'][s] > md.mesh.numberofvertices:
                    raise Exception('qmustatistics consistency check error: qmu.statistics.method[{}][\'steps\'][{}] should be < md.mesh.numberofvertices!'.format(i, s))
    # }}}

    def defaultoutputs(self, md):  # {{{
        outputs = []
        return outputs
    # }}}

    def marshall(self, prefix, md, fid):  # {{{
        if self.method[0]['name'] == 'None':
            WriteData(fid, prefix, 'name', 'md.qmu.statistics', 'data', 0, 'format', 'Boolean')
            statistics = 0
            return
        else:
            WriteData(fid, prefix, 'name', 'md.qmu.statistics', 'data', 1, 'format', 'Boolean')
            statistics = 1

        if statistics:
            WriteData(fid, prefix, 'name', 'md.qmu.statistics.nfiles_per_directory', 'data', self.nfiles_per_directory, 'format', 'Integer')
            WriteData(fid, prefix, 'name', 'md.qmu.statistics.ndirectories', 'data', self.ndirectories, 'format', 'Integer')
            WriteData(fid, prefix, 'name', 'md.qmu.statistics.numstatistics', 'data', len(self.method), 'format', 'Integer')
            for i in range(1, len(self.method) + 1):
                m = self.method[i - 1]
                WriteData(fid, prefix, 'name', 'md.qmu.statistics.method({}).name'.format(i), 'data', m['name'], 'format', 'String')
                WriteData(fid, prefix, 'data', m['fields'], 'name', 'md.qmu.statistics.method({}).fields'.format(i), 'format', 'StringArray')
                WriteData(fid, prefix, 'data', m['steps'], 'name', 'md.qmu.statistics.method({}).steps'.format(i), 'format', 'IntMat', 'mattype', 3)
                if m['name'] == 'Histogram':
                    WriteData(fid, prefix, 'name', 'md.qmu.statistics.method({}).nbins'.format(i), 'data', m['nbins'], 'format', 'Integer')
                elif m['name'] == 'MeanVariance':
                    pass # do nothing
                elif m['name'] == 'SampleSeries':
                    WriteData(fid, prefix, 'data', m['indices'], 'name', 'md.qmu.statistics.method({}).indices'.format(i), 'format', 'IntMat', 'mattype', 3)
                else:
                    raise Exception('qmustatistics marshall error message: unknown type ''{}'' for qmu.statistics.method[{}]'.format(m['name'], i))
    # }}}

    def extrude(self, md):  # {{{
        return self
    # }}}

    def addmethod(self, *args):  # {{{
        """ADDMETHOD - Add new, empty method or passed dict to self.method
        """
        nargs = len(args)
        if nargs == 0:
            self.method.append({})
        elif nargs == 1:
            self.method.append(args[0])
        else:
            raise Exception('Number of args should be 0 (appends empty dict to methods member) or 1 (appends passed dict to methods member)')
    # }}}
