import numpy as np
from dependent import dependent
from independent import independent
from fielddisplay import fielddisplay
from checkfield import checkfield
from WriteData import WriteData
from MatlabArray import *


class autodiff(object):
    """autodiff class definition

    Usage:
        autodiff = autodiff()
    """
    def __init__(self, *args):  # {{{
        self.isautodiff = False
        self.dependents = []
        self.independents = []
        self.driver = 'fos_forward'
        self.obufsize = np.nan
        self.lbufsize = np.nan
        self.cbufsize = np.nan
        self.tbufsize = np.nan
        self.gcTriggerMaxSize = np.nan
        self.gcTriggerRatio = np.nan
        self.tapeAlloc = np.nan
        self.outputTapeMemory = 0
        self.outputTime = 0
        self.enablePreaccumulation = 0
        if not len(args):
            self.setdefaultparameters()
        else:
            raise RuntimeError("constructor not supported")
    # }}}

    def __repr__(self):  # {{{
        s = '      automatic differentiation parameters:\n'
        s += '{}\n'.format(fielddisplay(self, 'isautodiff', "indicates if the automatic differentiation is activated"))
        s += '{}\n'.format(fielddisplay(self, 'dependents', "list of dependent variables"))
        s += '{}\n'.format(fielddisplay(self, 'independents', "list of independent variables"))
        s += '{}\n'.format(fielddisplay(self, 'driver', "ADOLC driver ('fos_forward' or 'fov_forward')"))
        s += '{}\n'.format(fielddisplay(self, 'obufsize', "Number of operations per buffer (== OBUFSIZE in usrparms.h)"))
        s += '{}\n'.format(fielddisplay(self, 'lbufsize', "Number of locations per buffer (== LBUFSIZE in usrparms.h)"))
        s += '{}\n'.format(fielddisplay(self, 'cbufsize', "Number of values per buffer (== CBUFSIZE in usrparms.h)"))
        s += '{}\n'.format(fielddisplay(self, 'tbufsize', "Number of taylors per buffer (<=TBUFSIZE in usrparms.h)"))
        s += '{}\n'.format(fielddisplay(self, 'gcTriggerRatio', "free location block sorting / consolidation triggered if the ratio between allocated and used locations exceeds gcTriggerRatio"))
        s += '{}\n'.format(fielddisplay(self, 'gcTriggerMaxSize', "free location block sorting / consolidation triggered if the allocated locations exceed gcTriggerMaxSize)"))
        s += '{}\n'.format(fielddisplay(self, 'tapeAlloc', 'Iteration count of a priori memory allocation of the AD tape'))
        s += '{}\n'.format(fielddisplay(self, 'outputTapeMemory', 'Write AD tape memory statistics to file ad_mem.dat'))
        s += '{}\n'.format(fielddisplay(self, 'outputTime', 'Write AD recording and evaluation times to file ad_time.dat'))
        s += '{}\n'.format(fielddisplay(self, 'enablePreaccumulation', 'Enable CoDiPack preaccumulation in augmented places'))

        return s
    # }}}

    def setdefaultparameters(self):  # {{{
        self.obufsize = 524288
        self.lbufsize = 524288
        self.cbufsize = 524288
        self.tbufsize = 524288
        self.gcTriggerRatio = 2.0
        self.gcTriggerMaxSize = 65536
        self.tapeAlloc = 15000000
        return self
    # }}}

    def checkconsistency(self, md, solution, analyses):  # {{{
        # Early return
        if not self.isautodiff:
            return md

        md = checkfield(md, 'fieldname', 'autodiff.obufsize', '>=', 524288)
        md = checkfield(md, 'fieldname', 'autodiff.lbufsize', '>=', 524288)
        md = checkfield(md, 'fieldname', 'autodiff.cbufsize', '>=', 524288)
        md = checkfield(md, 'fieldname', 'autodiff.tbufsize', '>=', 524288)
        md = checkfield(md, 'fieldname', 'autodiff.gcTriggerRatio', '>=', 2.0)
        md = checkfield(md, 'fieldname', 'autodiff.gcTriggerMaxSize', '>=', 65536)
        md = checkfield(md, 'fieldname', 'autodiff.tapeAlloc', '>=', 0)

        # Memory and time output
        md = checkfield(md, 'fieldname', 'autodiff.outputTapeMemory', 'numel', [1], 'values', [0, 1])
        md = checkfield(md, 'fieldname', 'autodiff.outputTime', 'numel', [1], 'values', [0, 1])

        # Memory reduction options
        md = checkfield(md, 'fieldname', 'autodiff.enablePreaccumulation', '>=', 0)

        # Driver value
        md = checkfield(md, 'fieldname', 'autodiff.driver', 'values', ['fos_forward', 'fov_forward', 'fov_forward_all', 'fos_reverse', 'fov_reverse', 'fov_reverse_all'])

        # Go through our dependents and independents and check consistency
        for dep in self.dependents:
            dep.checkconsistency(md, solution, analyses)
        for i, indep in enumerate(self.independents):
            indep.checkconsistency(md, i, solution, analyses, self.driver)

        return md
    # }}}

    def marshall(self, prefix, md, fid):  # {{{
        WriteData(fid, prefix, 'object', self, 'fieldname', 'isautodiff', 'format', 'Boolean')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'driver', 'format', 'String')

        # Early return
        if not self.isautodiff:
            WriteData(fid, prefix, 'data', False, 'name', 'md.autodiff.mass_flux_segments_present', 'format', 'Boolean')
            WriteData(fid, prefix, 'data', False, 'name', 'md.autodiff.keep', 'format', 'Boolean')
            return

        # Buffer sizes
        WriteData(fid, prefix, 'object', self, 'fieldname', 'obufsize', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'lbufsize', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'cbufsize', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'tbufsize', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'gcTriggerRatio', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'gcTriggerMaxSize', 'format', 'Double')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'tapeAlloc', 'format', 'Integer')

        # Output of memory and time
        WriteData(fid, prefix, 'object', self, 'fieldname', 'outputTapeMemory', 'format', 'Boolean')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'outputTime', 'format', 'Boolean')

        # Memory reduction options
        WriteData(fid, prefix, 'object', self, 'fieldname', 'enablePreaccumulation', 'format', 'Boolean')

        # Process dependent variables
        num_dependent_objects = len(self.dependents)
        WriteData(fid, prefix, 'data', num_dependent_objects, 'name', 'md.autodiff.num_dependent_objects', 'format', 'Integer')

        if num_dependent_objects:
            names = []
            for i, dep in enumerate(self.dependents):
                names.append(dep.name)

            WriteData(fid, prefix, 'data', names, 'name', 'md.autodiff.dependent_object_names', 'format', 'StringArray')

        # Process independent variables
        num_independent_objects = len(self.independents)
        WriteData(fid, prefix, 'data', num_independent_objects, 'name', 'md.autodiff.num_independent_objects', 'format', 'Integer')

        for indep in self.independents:
            WriteData(fid, prefix, 'data', indep.name, 'name', 'md.autodiff.independent_name', 'format', 'String')
            WriteData(fid, prefix, 'data', indep.min_parameters, 'name','md.autodiff.independent_min_parameters','format', 'DoubleMat', 'mattype', 3)
            WriteData(fid, prefix, 'data', indep.max_parameters, 'name', 'md.autodiff.independent_max_parameters', 'format', 'DoubleMat', 'mattype', 3)
            WriteData(fid, prefix, 'data', indep.control_scaling_factor, 'name', 'md.autodiff.independent_scaling_factor', 'format', 'Double')
            WriteData(fid, prefix, 'data', indep.control_size, 'name', 'md.autodiff.independent_control_size', 'format', 'Integer')

        # If driver is fos_forward, build index
        if strcmpi(self.driver, 'fos_forward'):
            index = 0

            for indep in self.independents:
                if not np.isnan(indep.fos_forward_index):
                    index += indep.fos_forward_index
                    break
                else:
                    if strcmpi(indep.type, 'scalar'):
                        index += 1
                    else:
                        index += indep.nods

            index -= 1  # get c-index numbering going
            WriteData(fid, prefix, 'data', index, 'name', 'md.autodiff.fos_forward_index', 'format', 'Integer')

        # If driver is fos_reverse, build index
        if strcmpi(self.driver, 'fos_reverse'):
            index = 0

            for dep in self.dependents:
                if not np.isnan(dep.fos_reverse_index):
                    index += dep.fos_reverse_index
                    break
                else:
                    index += 1

            index -= 1  # get c-index numbering going
            WriteData(fid, prefix, 'data', index, 'name', 'md.autodiff.fos_reverse_index', 'format', 'Integer')

        # If driver is fov_forward, build indices
        if strcmpi(self.driver, 'fov_forward'):
            indices = 0

            for indep in self.independents:
                if indep.fos_forward_index:
                    indices += indep.fov_forward_indices
                    break
                else:
                    if strcmpi(indep.type, 'scalar'):
                        indices += 1
                    else:
                        indices += indep.nods

            indices -= 1  # get c-indices numbering going
            WriteData(fid, prefix, 'data', indices, 'name', 'md.autodiff.fov_forward_indices', 'format', 'IntMat', 'mattype', 3)

        # Deal with mass fluxes
        mass_flux_segments = [dep.segments for dep in self.dependents if strcmpi(dep.name, 'MassFlux')]

        if mass_flux_segments:
            WriteData(fid, prefix, 'data', mass_flux_segments, 'name', 'md.autodiff.mass_flux_segments', 'format', 'MatArray')
            flag = True
        else:
            flag = False
        WriteData(fid, prefix, 'data', flag, 'name', 'md.autodiff.mass_flux_segments_present', 'format', 'Boolean')

        # Deal with trace keep on
        keep = False

        # From ADOLC userdoc:
        # The optional integer argument keep of trace on determines whether the 
        # numerical values of all active variables are recorded in a buffered 
        # temporary array or file called the taylor stack. This option takes 
        # effect if keep = 1 and prepares the scene for an immediately 
        # following gradient evaluation by a call to a routine implementing the 
        # reverse mode as described in the Section 4 and Section 5.
        #
        if len(self.driver) <= 3:
            keep = False  # there is no "_reverse" string within the driver string
        else:
            if strncmpi(self.driver[3:], '_reverse', 8):
                keep = True
            else:
                keep = False
        WriteData(fid, prefix, 'data', keep, 'name', 'md.autodiff.keep', 'format', 'Boolean')
    # }}}

        return
    # }}}
