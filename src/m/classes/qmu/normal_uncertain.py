import numpy as np

from MatlabArray import *
from MatlabFuncs import *
from fielddisplay import fielddisplay
from pairoptions import pairoptions
from partition_npart import *
from qmupart2npart import qmupart2npart


class normal_uncertain(object):
    """NORMAL_UNCERTAIN class definition

    Usage:
        [nuv] = normal_uncertain(
            'descriptor', descriptor,
            'mean', mean,
            'stddev', stddev,
            'partition', partition
            )

        where nuv is the normal_uncertain object returned by the constructor, 
        mean and stddev are self explanatory, and partition is the partition 
        vector for distributed variables. Can be a partition vector over 
        elements or vertices.

    Example:
        md.qmu.variables.rheology=normal_uncertain(
            'descriptor','RheologyBBar',
            'mean',1,
            'stddev',.05
            )
        md.qmu.variables.rheology=normal_uncertain(
            'descriptor','scaled_RheologyBBar',
            'mean',1,
            'stddev',.05,
            'partition',vpartition
            )
    """

    def __init__(self):  # {{{
        self.descriptor = ''
        self.mean       = np.nan
        self.stddev     = np.nan
        self.partition  = []
        self.nsteps     = 0
    # }}}

    @staticmethod
    def normal_uncertain(*args):  # {{{
        nargin = len(args)

        # create a default object
        if nargin == 0:
            return normal_uncertain()

        # copy the object
        elif nargin == 1:
            if isinstance(args[0], normal_uncertain):
                nuv = args[0]
            else:
                raise Exception('Object ' + str(args[0]) + ' is a ' + str(type(args[0])) + ' class object, not "normal_uncertain".')

        # create the object from the input
        else:
            nuv = normal_uncertain()

            #recover options:
            options = pairoptions(*args)

            #initialize fields:
            nuv.descriptor = options.getfieldvalue('descriptor')
            nuv.mean       = options.getfieldvalue('mean')
            nuv.stddev     = options.getfieldvalue('stddev')

            #if the variable is scaled, a partition vector should have been 
            #supplied, and that partition vector should have as many partitions 
            #as the mean and stddev vectors:
            if nuv.isscaled():
                nuv.partition = options.getfieldvalue('partition')
                nuv.nsteps = options.getfieldvalue('nsteps', 1)
                npart = qmupart2npart(nuv.partition)
                if npart != nuv.mean.shape[0]:
                    raise Exception("normal_uncertain constructor: for the scaled variable %s the row size of the mean field should be identical to the number of partitions" % nuv.descriptor)
                if npart != nuv.stddev.shape[0]:
                    raise Exception("normal_uncertain constructor: for the scaled variable %s the row size of the stddev field should be identical to the number of partitions" % nuv.descriptor)
                if nuv.nsteps != nuv.mean.shape[1]:
                    raise Exception("normal_uncertain constructor: for the scaled variable %s the col size of the mean field should be identical to the number of time steps" % nuv.descriptor)
                if nuv.nsteps != nuv.stddev.shape[1]:
                    raise Exception("normal_uncertain constructor: for the scaled variable %s the col size of the stddev field should be identical to the number of time steps" % nuv.descriptor)

        return [nuv] # Always return a list, so we have something akin to a MATLAB single row matrix
    # }}}

    def __repr__(self):  # {{{
        string = '   normal uncertain variable: '
        string = "%s\n%s" % (string, fielddisplay(self, 'descriptor', 'name tag'))
        string = "%s\n%s" % (string, fielddisplay(self, 'mean', 'pdf mean'))
        string = "%s\n%s" % (string, fielddisplay(self, 'stddev', 'pdf standard deviation'))
        if self.partition != []:
            string = "%s\n%s" % (string, fielddisplay(self, 'partition', 'partition vector defining where sampling will occur'))
        string = "%s\n%s" % (string, fielddisplay(self, 'nsteps', 'number of time steps'))

        return string
    # }}}

    def __len__(self):  # {{{
        if type(self.mean) in [list, np.ndarray]:
            return len(self.mean)
        else:
            return 1
    # }}}

    def checkconsistency(self, md, solution, analyses):  # {{{
        md = checkfield(md, 'field', self.mean, 'fieldname', 'normal_uncertain.mean', 'NaN', 1, 'Inf', 1, '>=', 0)
        md = checkfield(md, 'field', self.stddev, 'fieldname', 'normal_uncertain.stddev', 'NaN', 1, 'Inf', 1, '>=', 0)
        if self.isscaled():
            if self.partition == []:
                raise Exception("normal_uncertain is a scaled variable, but it's missing a partition vector")
            #better have a partition vector that has as many partitions as stddev's size:
            if self.stddev.shape[0] != partition_npart(self.partititon):
                raise Exception("normal_uncertain error message: row size of stddev and partition size should be identical")
            if self.mean.shape[0] != partition_npart(self.partition):
                raise Exception("normal_uncertain error message: row size of mean and partition size should be identical")
            #we need as many steps in stddev and mean as there are in time steps
            if self.stddev.shape[1] != self.nsteps:
                raise Exception("normal_uncertain error message: col size of stddev and partition size should be identical")
            if self.mean.shape[1] != self.nsteps:
                raise Exception("normal_uncertain error message: col size of mean and partition size should be identical")
            md = checkfield(md, 'field', self.partition, 'fieldname', 'normal_uncertain.partition', 'NaN', 1, 'Inf', 1, '>=', -1, 'numel', [md.mesh.numberofvertices, md.mesh.numberofvertices])
            if self.partition.shape[1] > 1:
                raise Exception("normal_uncertain error message: partition should be a column vector")
            partcheck = np.unique(self.partition)
            partmin = min(partcheck)
            partmax = max(partcheck)
            if partmax < -1:
                raise Exception("normal_uncertain error message: partition vector's min value should be -1 (for no partition), or start at 0")
            nmax = max(md.mesh.numberofelements, md.mesh.numberofvertices)
            if partmax > nmax:
                raise Exception("normal_uncertain error message: partition vector's values cannot go over the number of vertices or elements")
    # }}}

    #virtual functions needed by qmu processing algorithms
    #implemented:

    @staticmethod
    def prop_desc(nuv, dstr):  # {{{
        desc = ['' for i in range(np.size(nuv))]
        for i in range(np.size(nuv)):
            if nuv[i].descriptor != '' or type(nuv[i].descriptor) != str:
                desc[i] = str(nuv[i].descriptor)
            elif dstr != '':
                desc[i] = str(dstr) + str(string_dim(nuv, i, 'vector'))
            else:
                desc[i] = 'nuv' + str(string_dim(nuv, i, 'vector'))

        desc = allempty(desc)

        return desc
    # }}}

    @staticmethod
    def prop_mean(nuv):  # {{{
        mean = np.zeros(np.size(nuv))
        for i in range(np.size(nuv)):
            mean[i] = nuv[i].mean
        return mean
    # }}}

    @staticmethod
    def prop_stddev(nuv):  # {{{
        stddev = np.zeros(np.size(nuv))
        for i in range(np.size(nuv)):
            stddev[i] = nuv[i].stddev
        return stddev
    # }}}

    @staticmethod
    def prop_lower(nuv):  # {{{
        lower = []
        return lower
    # }}}

    @staticmethod
    def prop_upper(nuv):  # {{{
        upper = []
        return upper
    # }}}

    #default
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
    def prop_initpt(nuv):  # {{{
        initpt = []
        return initpt
    # }}}

    @staticmethod
    def prop_initst(nuv):  # {{{
        inist = []
        return inist
    # }}}

    @staticmethod
    def prop_stype(nuv):  # {{{
        stype = []
        return stype
    # }}}

    @staticmethod
    def prop_scale(nuv):  # {{{
        scale = []
        return scale
    # }}}

    #new methods:
    def isdistributed(self):  # {{{
        if strncmp(self.descriptor, 'distributed_', 12):
            return True
        else:
            return False
    # }}}
    
    def isscaled(self):  # {{{
        if strncmp(self.descriptor, 'scaled_', 7):
            return True
        else:
            return False
    # }}}

    @staticmethod
    def dakota_write(fidi, dvar):
        # possible namespace pollution, the above import seems not to work
        from vlist_write import vlist_write
        # collect only the variables of the appropriate class
        nuv = deepcopy(dvar)
        fields = fieldnames(nuv)
        for field in fields:
            if getattr(nuv, field)[0].__class__.__name__ != 'normal_uncertain':
                delattr(nuv, field)
        if len(nuv) > 0:
            vlist_write(fidi, 'normal_uncertain', 'nuv', nuv)
    # }}}

