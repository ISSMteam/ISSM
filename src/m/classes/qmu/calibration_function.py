import numpy as np
from rlist_write import *
from MatlabArray import *


class calibration_function(object):
    '''
  definition for the calibration_function class.

  [cf] = calibration_function.calibration_function(args)
   cf = calibration_function()

  where the required args are:
    descriptor    (char, description, '')
  and the optional args and defaults are:
    scale_type    (char, scaling type, 'none')
    scale         (double, scaling factor, 1.)
    weight        (double, weighting factor, 1.)

  note that zero arguments constructs a default instance, one
  argument of the class copies the instance, and one or more
  arguments constructs a new instance from the arguments.
'''
    def __init__(self):
        self.descriptor = ''
        self.scale_type = 'none'
        self.scale = 1.
        self.weight = 1.

    @staticmethod
    def calibration_function(*args):
        nargin = len(args)

    # create a default object
        if nargin == 0:
            return calibration_function()

    # copy the object or create the object from the input
        else:
            if (nargin == 1) and isinstance(args[0], calibration_function):
                cf = args[0]
            else:
                asizec = array_size(*args[0:min(nargin, 4)])
                cf = [calibration_function() for i in range(asizec[0]) for j in range(asizec[1])]

            for i in range(np.size(cf)):
                if (np.size(args[0]) > 1):
                    cf[i].descriptor = args[0][i]
                else:
                    cf[i].descriptor = str(args[0]) + string_dim(cf, i, 'vector')

            if nargin >= 2:
                for i in range(np.size(cf)):
                    cf[i].scale_type = str(args[1])
            if nargin >= 3:
                for i in range(np.size(cf)):
                    cf[i].scale = args[2]
            if nargin >= 4:
                for i in range(np.size(cf)):
                    cf[i].weight = args[3]
            if nargin > 4:
                print('WARNING: calibration_function:extra_arg: Extra arguments for object of class ' + str(type(cf)) + '.')

        return cf

    def __repr__(self):
        # display the object
        string = '\n'
        string += 'class "calibration_function" object = \n'
        string += '    descriptor: ' + str(self.descriptor) + '\n'
        string += '    scale_type: ' + str(self.scale_type) + '\n'
        string += '         scale: ' + str(self.scale) + '\n'
        string += '        weight: ' + str(self.weight) + '\n'
        return string

    # from here on, cf is either a single, or a 1d vector of, calibration_function
    @staticmethod
    def prop_desc(cf, dstr):
        if type(cf) not in [list, np.ndarray]:
            if cf.descriptor != '' or type(cf.descriptor) != str:
                desc = str(cf.descriptor)
            elif dstr != '':
                desc = str(dstr)
            else:
                desc = 'cf'
            return desc

        desc = ['' for i in range(np.size(cf))]
        for i in range(np.size(cf)):
            if cf[i].descriptor != '' or type(cf[i].descriptor) != str:
                desc[i] = str(cf[i].descriptor)
            elif dstr != '':
                desc[i] = str(dstr) + str(string_dim(cf, i, 'vector'))
            else:
                desc[i] = 'cf' + str(string_dim(cf, i, 'vector'))

        desc = allempty(desc)
        return desc

    @staticmethod
    def prop_stype(cf):
        stype = ''
        return stype

    @staticmethod
    def prop_weight(cf):
        weight = []
        return weight

    @staticmethod
    def prop_lower(cf):
        lower = []
        return lower

    @staticmethod
    def prop_upper(cf):
        upper = []
        return upper

    @staticmethod
    def prop_target(cf):
        target = []
        return target

    @staticmethod
    def prop_scale(cf):
        scale = []
        return scale

    @staticmethod
    def dakota_write(fidi, dresp, rdesc):
        # collect only the responses of the appropriate class
        cf = [struc_class(i, 'calibration_function', 'cf') for i in dresp]

        # write responses
        rdesc = rlist_write(fidi, 'calibration_terms', 'calibration_function', cf, rdesc)
        return rdesc

    @staticmethod
    def dakota_rlev_write(fidi, dresp, params):
        return
