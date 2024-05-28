import numpy as np

from collections import OrderedDict
from dakota_method import *
from fielddisplay import fielddisplay
from helpers import *
from IssmConfig import *
from MatlabFuncs import *
from qmustatistics import qmustatistics
from WriteData import WriteData


class qmu(object):
    """QMU class definition

    Usage:
        qmu = qmu()
    """

    def __init__(self):  # {{{
        self.isdakota = 0
        self.output = 0
        self.variables = OrderedStruct()
        self.correlation_matrix = []
        self.responses = OrderedStruct()
        self.method = OrderedDict()
        self.params = OrderedStruct()
        self.statistics = qmustatistics()
        self.results = OrderedDict()
        self.numberofresponses = 0
        self.variabledescriptors = []
        self.variablepartitions = []
        self.variablepartitions_npart = []
        self.variablepartitions_nt = []
        self.responsedescriptors = []
        self.responsepartitions = []
        self.responsepartitions_npart = []
        self.responsepartitions_nt = []
        self.mass_flux_profile_directory = float('NaN')
        self.mass_flux_profiles = float('NaN')
        self.mass_flux_segments = []
        self.adjacency = float('NaN')
        self.vertex_weight = float('NaN')

        self.setdefaultparameters()
    # }}}
    def __repr__(self):  # {{{
        s = '   qmu parameters:\n'
        s += '{}\n'.format(fielddisplay(self, 'isdakota', 'is QMU analysis activated?'))
        s += '{}\n'.format(fielddisplay(self, 'output', 'are we outputting ISSM results, default is 0'))
        maxlen = 0
        s += '         variables:  (arrays of each variable class)\n'

    # OrderedStruct's iterator returns individual name / array - of - functions pairs
        for variable in self.variables:
            fname = variable[0]
            maxlen = max(maxlen, len(fname))
            size = np.shape(variable[1])
            a = size[0]
            b = 1 if len(size) < 2 else size[1]
            s += "            %-*s:    [%ix%i]    '%s'\n" % (maxlen + 1, fname, a, b, type(variable[1][0]))

        s += "         responses:  (arrays of each response class)\n"
        for response in self.responses:
            fname = response[0]
            maxlen = max(maxlen, len(fname))
            size = np.shape(response[1])
            a = size[0]
            b = 1 if len(size) < 2 else size[1]
            s += "            %-*s:    [%ix%i]    '%s'\n" % (maxlen + 1, fname, a, b, type(response[1][0]))

        s += '{}\n'.format(fielddisplay(self, 'numberofresponses', 'number of responses'))

        if type(self.method) != OrderedDict:
            self.method = [self.method]
    # self.method must be iterable
        for method in self.method:
            if isinstance(method, dakota_method):
                s += "            method :    '%s'\n" % (method.method)

    # params could have a number of forms (mainly 1 struct or many)
        if type(self.params) == OrderedStruct:
            params = [self.params]
        else:
            params = np.hstack(np.atleast_1d(np.array(self.params)))
        for param in params:
            print(type(param))
            print(param)
            s += "         params:  (array of method-independent parameters)\n"
            fnames = vars(param)
            maxlen = 0
            for fname in fnames:
                maxlen = max(maxlen, len(fname))

            for fname in fnames:
                s += "            %-*s: %s\n" % (maxlen + 1, fname, str(getattr(param, fname)))

    # results could be have a number of forms (mainly 1 struct or many)
        results = np.hstack(np.atleast_1d(np.array(self.results)))
        for result in results:
            s += "         results:  (information from dakota files)\n"
            fnames = vars(result)
            maxlen = 0
            for fname in fnames:
                maxlen = max(maxlen, len(fname))

            for fname in fnames:
                size = np.shape(response[1])
                a = size[0]
                b = 0 if len(size) < 2 else size[1]
                size = np.shape(getattr(result, fname))
                s += "            %-*s:    [%ix%i]    '%s'\n" % (maxlen + 1, fname, a, b, type(getattr(result, fname)))

        s += '{}\n'.format(fielddisplay(self, 'variablepartitions', ''))
        s += '{}\n'.format(fielddisplay(self, 'variablepartitions_npart', ''))
        s += '{}\n'.format(fielddisplay(self, 'variablepartitions_nt', ''))
        s += '{}\n'.format(fielddisplay(self, 'variabledescriptors', ''))
        s += '{}\n'.format(fielddisplay(self, 'responsedescriptors', ''))
        s += '{}\n'.format(fielddisplay(self, 'method', 'array of dakota_method class'))
        s += '{}\n'.format(fielddisplay(self, 'mass_flux_profile_directory', 'directory for mass flux profiles'))
        s += '{}\n'.format(fielddisplay(self, 'mass_flux_profiles', 'list of mass_flux profiles'))
        s += '{}\n'.format(fielddisplay(self, 'mass_flux_segments', ''))
        s += '{}\n'.format(fielddisplay(self, 'adjacency', ''))
        s += '{}\n'.format(fielddisplay(self, 'vertex_weight', 'weight applied to each mesh vertex'))
        return s
    # }}}

    def extrude(self, md):  # {{{
        return self
    # }}}

    def setdefaultparameters(self):  # {{{
        return self
    # }}}

    def checkconsistency(self, md, solution, analyses):  # {{{
        #Early return
        if not md.qmu.isdakota:
            return

        version = IssmConfig('_DAKOTA_VERSION_')
        version = float(version[0])

        if version < 6:
            if not md.qmu.params.evaluation_concurrency == 1:
                md.checkmessage("concurrency should be set to 1 when running dakota in library mode")
        else:
            if not strcmpi(self.params.evaluation_scheduling, 'master'):
                md.checkmessage('evaluation_scheduling in qmu.params should be set to "master"')

            if md.cluster.nprocs() <= 1:
                md.checkmessage('in parallel library mode, Dakota needs to run on at least 2 cpus, 1 cpu for the master, 1 cpu for the slave. Modify md.cluster.np accordingly.')

            if self.params.processors_per_evaluation < 1:
                md.checkmessage('in parallel library mode, Dakota needs to run at least one slave on one cpu (md.qmu.params.processors_per_evaluation >= 1)!')

            if np.mod(md.cluster.nprocs() - 1, self.params.processors_per_evaluation):
                #md.checkmessage('in parallel library mode, the requirement is for md.cluster.np = md.qmu.params.processors_per_evaluation * number_of_slaves, where number_of_slaves will automatically be determined by Dakota. Modify md.cluster.np accordingly')
                pass

        # Go through variables and check for consistency
        fv = fieldnames(self.variables)
        for i in range(len(fv)):
            variable = getattr(self.variables, fv[i])
            if hasattr(variable, 'checkconsistency'):
                variable.checkconsistency(md, solution, analyses)

        # Go through variables and check that we have normal uncertains first,
        # then uniform uncertains and finally histogram_bin_uncertain. Indeed,
        # Dakota will order them this waym, and when we send partitions for
        # scaled variables, they better show up in the order Dakota is feeding
        # them to us in InputUpdateFromDakotax!
        fv = fieldnames(self.variables)
        classlist = []
        for i in range(len(fv)):
            classlist.append(self.variables[fv[i]].__class__.__name__)
        n = 0
        u = 0
        h = 0
        for i in range(len(classlist)):
            if classlist[i] == 'normal_uncertain':
                if u != 0 or h != 0:
                    raise Exception('normal uncertain variables should be declared before uniform and hhistogram_bin uncertain variables')
                else:
                    n = 1
            if classlist[i] == 'uniform_uncertain':
                if h != 0:
                    raise Exception('uniform_uncertain variables should be declared before histogram_bin uncertain variables')
                else:
                    u = 1
            if classlist[i] == 'histogram_bin_uncertain':
                h = 1

        return md
    # }}}

    def marshall(self, prefix, md, fid):  # {{{
        WriteData(fid, prefix, 'object', self, 'fieldname', 'isdakota', 'format', 'Boolean')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'output', 'format', 'Boolean')
        if not self.isdakota:
            WriteData(fid, prefix, 'data', False, 'name', 'md.qmu.mass_flux_segments_present', 'format', 'Boolean')
            return
        WriteData(fid, prefix, 'data', self.method.params.samples, 'name', 'md.qmu.method.params.samples', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'numberofresponses', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'variabledescriptors', 'format', 'StringArray')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'variablepartitions', 'format', 'MatArray')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'variablepartitions_npart', 'format', 'IntMat', 'mattype', 3)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'variablepartitions_nt', 'format', 'IntMat', 'mattype', 3)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'responsedescriptors', 'format', 'StringArray')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'responsepartitions', 'format', 'MatArray')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'responsepartitions_npart', 'format', 'IntMat', 'mattype', 3)
        if not isempty(self.mass_flux_segments):
            WriteData(fid, prefix, 'data', self.mass_flux_segments, 'name', 'md.qmu.mass_flux_segments', 'format', 'MatArray')
            flag = True
        else:
            flag = False
        WriteData(fid, prefix, 'data', flag, 'name', 'md.qmu.mass_flux_segments_present', 'format', 'Boolean')
        self.statistics.marshall(prefix, md, fid)
    # }}}
