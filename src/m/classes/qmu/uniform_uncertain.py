import numpy as np

from MatlabArray import *
from MatlabFuncs import *
from fielddisplay import fielddisplay
from pairoptions import pairoptions
from partition_npart import *
from qmupart2npart import qmupart2npart


class uniform_uncertain(object):
    '''
    UNIFORM_UNCERTAIN class definition

    Usage:
        [uuv] = uniform_uncertain(
            'descriptor', descriptor,
            'lower', lower,
            'upper', upper,
            'partition', partition
            )

        where uuv is the uniform_uncertain object returned by the constructor, 
        lower and upper are the pdf distribution bounds, and partition is the 
        partition vector for distributed variables. Can be a partition vector 
        over elements or vertices.

    Example:
        md.qmu.variables.rheology = uniform_uncertain(
            'descriptor', 'RheologyBBar',
            'lower', 1e8,
            'upper', 1e9
            )
        md.qmu.variables.rheology = uniform_uncertain(
            'descriptor', 'RheologyBBar',
            'lower', 1e8,
            'upper', 1e9,
            'partition', vpartition
            )
    '''
    def __init__(self):  # {{{
        self.descriptor = ''
        self.lower      = -np.inf
        self.upper      = np.inf
        self.partition  = []
        self.nsteps     = 0
    # }}}

    @staticmethod
    def uniform_uncertain(*args):  # {{{
        nargin = len(args)

        # create a default object
        if nargin == 0:
            return uniform_uncertain()

        # copy the object
        elif nargin == 1:
            if isinstance(args[0], uniform_uncertain):
                uuv = args[0]
            else:
                raise Exception('Object ' + str(args[0]) + ' is a ' + str(type(args[0])) + ' class object, not "uniform_uncertain".')

        # create the object from the input
        else:
            uuv = uniform_uncertain()

            #recover options:
            options = pairoptions(*args)

            #initialize fields:
            uuv.descriptor = options.getfieldvalue('descriptor')
            uuv.lower      = options.getfieldvalue('lower')
            uuv.upper      = options.getfieldvalue('upper')

            #if the variable is scaled, a partition vector should have been 
            #supplied, and  that partition vector should have as many 
            #partitions as the lower and upper vectors:
            if uuv.isscaled():
                uuv.partition = options.getfieldvalue('partition')
                uuv.nsteps = options.getfieldvalue('nsteps', 1)
                npart = qmupart2npart(uuv.partition)
                if npart != uuv.upper.shape[0]:
                    raise Exception("uniform_uncertain constructor: for the scaled variable %s the upper field is not currently a vector of values for all the partitions described in the partition vector" % uuv.descriptor)
                if npart != uuv.lower.shape[0]:
                    raise Exception("uniform_uncertain constructor: for the scaled variable %s the lower field is not currently a vector of values for all the partitions described in the partition vector" % uuv.descriptor)
                if uuv.nsteps != uuv.upper.shape[1]:
                    raise Exception("uniform_uncertain constructor: for the scaled variable %s the col size of the upper field should be identical to the number of time steps" % uuv.descriptor)
                if uuv.nsteps != uuv.lower.shape[1]:
                    raise Exception("uniform_uncertain constructor: for the scaled variable %s the col size of the lower field should be identical to the number of time steps" % uuv.descriptor)

        return [uuv] # Always return a list, so we have something akin to a MATLAB single row matrix
    # }}}

    def __repr__(self):  # {{{
        string = '   uniform uncertain variable: '
        string = "%s\n%s" % (string, fielddisplay(self, 'descriptor', 'name tag'))
        string = "%s\n%s" % (string, fielddisplay(self, 'lower', 'pdf lower bound'))
        string = "%s\n%s" % (string, fielddisplay(self, 'upper', 'pdf upper bound'))
        if self.partition != []:
            string = "%s\n%s" % (string, fielddisplay(self, 'partition', 'partition vector defining where sampling will occur'))
        string = "%s\n%s" % (string, fielddisplay(self, 'nsteps', 'number of time steps'))

        return string
    # }}}

    def __len__(self):  # {{{
        if type(self.lower) in [list, np.ndarray]:
            return len(self.lower)
        else:
            return 1
    # }}}
    
    def checkconsistency(self, md, solution, analyses):  # {{{
        md = checkfield(md, 'field', self.upper, 'fieldname', 'uniform_uncertain.upper', 'NaN', 1, 'Inf', 1, '>', self.lower, 'numel', len(self.lower))
        md = checkfield(md, 'field', self.lower, 'fieldname', 'uniform_uncertain.upper', 'NaN', 1, 'Inf', 1, '<', self.upper, 'numel', len(self.upper))
        if self.isscaled():
            if self.partition == []:
                raise Exception("uniform_uncertain is a scaled variable, but it's missing a partition vector")
            #better have a partition vector that has as many partitions as 
            #upper and lower's size:
            if self.upper.shape[0] != partition_npart(self.partititon):
                raise Exception("uniform_uncertain error message: row size of upper and partition size should be identical")
            if self.lower.shape[0] != partition_npart(self.partition):
                raise Exception("uniform_uncertain error message: row size of lower and partition size should be identical")
            #we need as steps in upper and lower as there are time steps
            if self.stddev.shape[1] != self.nsteps:
                raise Exception("uniform_uncertain error message: col size of upper and partition size should be identical")
            if self.mean.shape[1] != self.nsteps:
                raise Exception("uniform_uncertain error message: col size of lower and partition size should be identical")
            md = checkfield(md, 'field', self.partition, 'fieldname', 'uniform_uncertain.partition', 'NaN', 1, 'Inf', 1, '>=', -1, 'numel', [md.mesh.numberofvertices, md.mesh.numberofvertices])
            if self.partition.shape[1] > 1:
                raise Exception("uniform_uncertain error message: partition should be a column vector")
            partcheck = np.unique(self.partition)
            partmin = min(partcheck)
            partmax = max(partcheck)
            if partmax < -1:
                raise Exception("uniform_uncertain error message: partition vector's min value should be -1 (for no partition), or start at 0")
            nmax = max(md.mesh.numberofelements, md.mesh.numberofvertices)
            if partmax > nmax:
                raise Exception("uniform_uncertain error message: partition vector's values cannot go over the number of vertices or elements")
    # }}}

    #virtual functions needed by qmu processing algorithms:
    #implemented:

    @staticmethod
    def prop_desc(uuv, dstr):  # {{{
        desc = ['' for i in range(np.size(uuv))]
        for i in range(np.size(uuv)):
            if uuv[i].descriptor != '' or type(uuv[i].descriptor) != str:
                desc[i] = str(uuv[i].descriptor)
            elif dstr != '':
                desc[i] = str(dstr) + str(string_dim(uuv, i, 'vector'))
            else:
                desc[i] = 'uuv' + str(string_dim(uuv, i, 'vector'))

            desc = allempty(desc)

        return desc
    # }}}

    @staticmethod
    def prop_stddev(uuv):  # {{{
        stddev = []
        return stddev
    # }}}

    @staticmethod
    def prop_mean(uuv):  # {{{
        mean = []
        return mean
    # }}}

    @staticmethod
    def prop_lower(uuv):  # {{{
        lower = np.zeros(np.size(uuv))
        for i in range(np.size(uuv)):
            lower[i] = uuv[i].lower

        lower = allequal(lower, -np.inf)

        return lower
    # }}}

    @staticmethod
    def prop_upper(uuv):  # {{{
        upper = np.zeros(np.size(uuv))
        for i in range(np.size(uuv)):
            upper[i] = uuv[i].upper

        #upper = allequal(upper, np.inf)

        return upper
    # }}}

    @staticmethod
    def prop_abscissas(hbu):  # {{{
        abscissas = []
        return abscissas
    # }}}

    @staticmethod
    def prop_pairs_per_variable(hbu):  # {{{
        pairs_per_variable = []
        return pairs_per_variable
    # }}}

    @staticmethod
    def prop_counts(hbu):  # {{{
        counts = []
        return counts
    # }}}

    @staticmethod
    def prop_initpt(uuv):  # {{{
        initpt = []
        return initpt
    # }}}

    @staticmethod
    def prop_initst(uuv):  # {{{
        initst = []
        return initst
    # }}}

    @staticmethod
    def prop_stype(uuv):  # {{{
        stype = []
        return stype
    # }}}

    @staticmethod
    def prop_scale(uuv):  # {{{
        scale = []
        return scale
    # }}}

    #new methods:
    def isscaled(self):  # {{{
        if strncmp(self.descriptor, 'scaled_', 7):
            return True
        else:
            return False
    # }}}

    @staticmethod
    def dakota_write(fidi, dvar):  # {{{
        # possible namespace pollution, the above import seems not to work
        from vlist_write import vlist_write
        # collect only the variables of the appropriate class
        uuv = deepcopy(dvar)
        fields = fieldnames(uuv)
        for field in fields:
            if getattr(uuv, field)[0].__class__.__name__ != 'uniform_uncertain':
                delattr(uuv, field)
        if len(uuv) > 0:
            vlist_write(fidi, 'uniform_uncertain', 'uuv', uuv)
    # }}}
