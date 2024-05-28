import numpy as np
from vlist_write import *
from MatlabArray import *


class continuous_design(object):
    '''
  definition for the continuous_design class.

  [cdv] = continuous_design.continuous_design(args)
   cdv = continuous_design()

  where the required args are:
    descriptor    (char, description, '')
    initpt        (double, initial point, 0.)
  and the optional args and defaults are:
    lower         (double, lower bound, -Inf)
    upper         (double, upper bound, Inf)
    scale_type    (char, scaling type, 'none')
    scale         (double, scaling factor, 1.)

  note that zero arguments constructs a default instance, one
  argument of the class copies the instance, and two or more
  arguments constructs a new instance from the arguments.
'''

    def __init__(self):
        self.descriptor = ''
        self.initpt = 0.
        self.lower = -np.inf
        self.upper = np.inf
        self.scale_type = 'none'
        self.scale = 1.

    @staticmethod
    def continuous_design(*args):
        nargin = len(args)

    #  create a default object
        if nargin == 0:
            return continuous_design()

    #  copy the object
        if nargin == 1:
            if isinstance(args[0], continuous_design):
                cdv = args[0]
            else:
                raise RuntimeError('Object is a ' + str(type(args[0])) + ' class object, not "continuous_design".')

    #  create the object from the input
        else:
            shapec = array_size(*args[0:min(nargin, 6)])
            cdv = [continuous_design() for i in range(shapec[0]) for j in range(shapec[1])]
    # in case cdv doesn't use array - like args
    #cdv = continuous_design()

            for i in range(np.size(cdv)):
                if (np.size(args[0]) > 1):
                    cdv[i].descriptor = args[0][i]
                else:
                    cdv[i].descriptor = str(args[0]) + string_dim(cdv, i, 'vector')

            if (nargin >= 2):
                for i in range(np.size(cdv)):
                    if (np.size(args[1]) > 1):
                        cdv[i].initpt = args[1][i]
                    else:
                        cdv[i].initpt = args[1]

            if (nargin >= 3):
                for i in range(np.size(cdv)):
                    if (np.size(args[2]) > 1):
                        cdv[i].lower = args[2][i]
                    else:
                        cdv[i].lower = args[2]

            if (nargin >= 4):
                for i in range(np.size(cdv)):
                    if (np.size(args[3]) > 1):
                        cdv[i].upper = args[3][i]
                    else:
                        cdv[i].upper = args[3]

            if (nargin >= 5):
                for i in range(np.size(cdv)):
                    if (np.size(args[4]) > 1):
                        cdv[i].scale_type = args[4][i]
                    else:
                        cdv[i].scale_type = str(args[4])

            if (nargin >= 6):
                for i in range(np.size(cdv)):
                    if (np.size(args[5]) > 1):
                        cdv[i].scale = args[5][i]
                    else:
                        cdv[i].scale = args[5]

            if (nargin > 6):
                print('WARNING: continuous_design:extra_arg: Extra arguments for object of class ' + str(type(cdv)) + '.')

        return cdv

    def __repr__(self):
        #  display the object
        string = '\n'
        string += 'class "continuous_design" object = \n'
        string += '    descriptor: ' + str(self.descriptor) + '\n'
        string += '        initpt: ' + str(self.initpt) + '\n'
        string += '         lower: ' + str(self.lower) + '\n'
        string += '         upper: ' + str(self.upper) + '\n'
        string += '    scale_type: ' + str(self.scale_type) + '\n'
        string += '         scale: ' + str(self.scale) + '\n'

        return string

    @staticmethod
    def prop_desc(cdv, dstr):
        if type(cdv) not in [list, np.ndarray]:
            cdv = [cdv]
    # in case cdv doesn't use array - like args
    #if cdv.descriptor != '' or type(cdv.descriptor) != str:
    #desc = str(cdv.descriptor)
    #elif dstr != '':
    #desc = str(dstr)
    #else:
    #desc = 'cdv'
    #return desc

        desc = ['' for i in range(np.size(cdv))]
        for i in range(np.size(cdv)):
            if cdv[i].descriptor != '' or type(cdv[i].descriptor) != str:
                desc[i] = str(cdv[i].descriptor)
            elif dstr != '':
                desc[i] = str(dstr) + str(string_dim(cdv, i, 'vector'))
            else:
                desc[i] = 'cdv' + str(string_dim(cdv, i, 'vector'))

        desc = allempty(desc)

        return desc

    @staticmethod
    def prop_initpt(cdv):
        if type(cdv) not in [list, np.ndarray]:
            return cdv.initpt

        initpt = np.zeros(np.size(cdv))
        for i in range(np.size(cdv)):
            initpt[i] = cdv[i].initpt

        initpt = allequal(initpt, 0.)

        return initpt

    @staticmethod
    def prop_lower(cdv):
        if type(cdv) not in [list, np.ndarray]:
            return cdv.lower

        lower = np.zeros(np.size(cdv))
        for i in range(np.size(cdv)):
            lower[i] = cdv[i].lower

        lower = allequal(lower, -np.inf)

        return lower

    @staticmethod
    def prop_upper(cdv):
        if type(cdv) not in [list, np.ndarray]:
            return cdv.upper

        upper = np.zeros(np.size(cdv))
        for i in range(np.size(cdv)):
            upper[i] = cdv[i].upper

        upper = allequal(upper, np.inf)

        return upper

    @staticmethod
    def prop_mean(cdv):
        mean = []
        return mean

    @staticmethod
    def prop_stddev(cdv):
        stddev = []
        return stddev

    @staticmethod
    def prop_initst(cdv):
        initst = []
        return initst

    @staticmethod
    def prop_stype(cdv):
        if type(cdv) not in [list, np.ndarray]:
            return str(cdv.scale_type)

        stype = np.empty(np.size(cdv))
        stype.fill(0.0)
        for i in range(np.size(cdv)):
            stype[i] = str(cdv[i].scale_type)

        stype = allequal(stype, 'none')

        return stype

    @staticmethod
    def prop_scale(cdv):
        if type(cdv) not in [list, np.ndarray]:
            return cdv.scale

        scale = np.zeros(np.size(cdv))
        for i in range(np.size(cdv)):
            scale[i] = cdv[i].scale

        scale = allequal(scale, 1.)

        return scale

    @staticmethod
    def dakota_write(fidi, dvar):
        #  collect only the variables of the appropriate class
        cdv = [struc_class(i, 'continuous_design', 'cdv') for i in dvar]

    #  write variables
        vlist_write(fidi, 'continuous_design', 'cdv', cdv)
