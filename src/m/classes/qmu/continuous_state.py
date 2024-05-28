import numpy as np
from vlist_write import *
from MatlabArray import *


class continuous_state(object):
    '''
  definition for the continuous_state class.

  [csv] = continuous_state.continuous_state(args)
   csv = continuous_state()

  where the required args are:
    descriptor    (char, description, '')
    initst        (double, initial state, 0.)
  and the optional args and defaults are:
    lower         (double, lower bound, -Inf)
    upper         (double, upper bound, Inf)

  note that zero arguments constructs a default instance, one
  argument of the class copies the instance, and two or more
  arguments constructs a new instance from the arguments.
'''

    def __init__(self):
        self.descriptor = ''
        self.initst = 0.
        self.lower = -np.inf
        self.upper = np.inf

    @staticmethod
    def continuous_state(*args):
        nargin = len(args)

        #  create a default object
        if nargin == 0:
            return continuous_state()

        #  copy the object
        if nargin == 1:
            if isinstance(args[0], continuous_state):
                csv = args[0]
            else:
                raise RuntimeError('Object is a ' + str(type(args[0])) + ' class object, not "continuous_state".')

        #  create the object from the input
        else:
            shapec = np.shape(*args[0:min(nargin, 4)])
            csv = [continuous_state() for i in range(shapec[0]) for j in range(shapec[1])]

            for i in range(np.size(csv)):
                if (np.size(args[0]) > 1):
                    csv[i].descriptor = args[0][i]
                else:
                    csv[i].descriptor = str(args[0]) + string_dim(csv, i, 'vector')

            if (nargin >= 2):
                for i in range(np.size(csv)):
                    if (np.size(args[1]) > 1):
                        csv[i].initst = args[1][i]
                    else:
                        csv[i].initst = args[1]

            if (nargin >= 3):
                for i in range(np.size(csv)):
                    if (np.size(args[2]) > 1):
                        csv[i].lower = args[2][i]
                    else:
                        csv[i].lower = args[2]

            if (nargin >= 4):
                for i in range(np.size(csv)):
                    if (np.size(args[3]) > 1):
                        csv[i].upper = args[3][i]
                    else:
                        csv[i].upper = args[3]

            if (nargin > 4):
                print('continuous_state:extra_arg', 'Extra arguments for object of class ' + str(type(csv)) + '.')

        return csv

    def __repr__(self):
        #  display the object
        string = '\n'
        string += 'class "continuous_state" object = \n'
        string += '    descriptor: ' + str(self.descriptor) + '\n'
        string += '        initst: ' + str(self.initst) + '\n'
        string += '         lower: ' + str(self.lower) + '\n'
        string += '         upper: ' + str(self.upper) + '\n'

        return string

    @staticmethod
    def prop_desc(csv, dstr):
        if type(csv) not in [list, np.ndarray]:
            csv = [csv]

        desc = ['' for i in range(np.size(csv))]
        for i in range(np.size(csv)):
            if csv[i].descriptor != '' or type(cdv[i].descriptor) != str:
                desc[i] = str(csv[i].descriptor)
            elif dstr != '':
                desc[i] = str(dstr) + str(string_dim(csv, i, 'vector'))
            else:
                desc[i] = 'csv' + str(string_dim(csv, i, 'vector'))

        desc = allempty(desc)

        return desc

    @staticmethod
    def prop_initpt(csv):
        initpt = []
        return initpt

    @staticmethod
    def prop_lower(csv):
        if type(csv) not in [list, np.ndarray]:
            return csv.lower

        lower = np.zeros(np.size(csv))
        for i in range(np.size(csv)):
            lower[i] = csv[i].lower

        lower = allequal(lower, -np.inf)

        return lower

    @staticmethod
    def prop_upper(csv):
        if type(csv) not in [list, np.ndarray]:
            return csv.upper

        upper = np.zeros(np.size(csv))
        for i in range(np.size(csv)):
            upper[i] = csv[i].upper

        upper = allequal(upper, np.inf)

        return upper

    @staticmethod
    def prop_mean(csv):
        mean = []
        return mean

    @staticmethod
    def prop_stddev(csv):
        stddev = []
        return stddev

    @staticmethod
    def prop_initst(csv):
        if type(csv) not in [list, np.ndarray]:
            return csv.initst

        initst = np.zeros(np.size(csv))
        for i in range(np.size(csv)):
            initst[i] = csv[i].initst

        initst = allequal(initst, 0.)

        return initst

    @staticmethod
    def prop_stype(csv):
        stype = ''
        return stype

    @staticmethod
    def prop_scale(csv):
        scale = []
        return scale

    @staticmethod
    def dakota_write(fidi, dvar):
        #  collect only the variables of the appropriate class
        csv = [struc_class(i, 'continuous_state', 'csv') for i in dvar]

        #  write variables
        vlist_write(fidi, 'continuous_state', 'csv', csv)
