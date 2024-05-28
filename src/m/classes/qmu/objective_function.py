import numpy as np
from rlist_write import *
from MatlabArray import *


class objective_function(object):
    '''
  definition for the objective_function class.

  [of] = objective_function.objective_function(args)
   of = objective_function()

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
    def objective_function(*args):
        nargin = len(args)

    #  create a default object
        if nargin == 0:
            return objective_function()

    #  copy the object or create the object from the input
        else:
            if (nargin == 1) and isinstance(args[0], objective_function):
                of = args[0]
            else:
                shapec = array_size(*args[0:min(nargin, 4)])
                of = [objective_function() for i in range(shapec[0]) for j in range(shapec[1])]

                for i in range(np.size(of)):
                    if (np.size(args[0]) > 1):
                        of[i].descriptor = args[0][i]
                    else:
                        of[i].descriptor = str(args[0]) + string_dim(of, i, 'vector')

                if (nargin >= 2):
                    for i in range(np.size(of)):
                        if (np.size(args[1]) > 1):
                            of[i].scale_type = args[1][i]
                        else:
                            of[i].scale_type = str(args[1])

                if (nargin >= 3):
                    for i in range(np.size(of)):
                        if (np.size(args[2]) > 1):
                            of[i].scale = args[2][i]
                        else:
                            of[i].scale = args[2]

                if (nargin >= 4):
                    for i in range(np.size(of)):
                        if (np.size(args[3]) > 1):
                            of[i].weight = args[3][i]
                        else:
                            of[i].weight = args[3]

                if (nargin > 4):
                    print('WARNING: objective_function:extra_arg Extra arguments for object of class ' + str(type(of)) + '.')

        return of

    def __repr__(self):
        #  display the object
        string = '\n'
        string += 'class "objective_function" object = \n'
        string += '    descriptor: ' + str(self.descriptor) + '\n'
        string += '    scale_type: ' + str(self.scale_type) + '\n'
        string += '         scale: ' + str(self.scale) + '\n'
        string += '        weight: ' + str(self.weight) + '\n'
        return string

    @staticmethod
    def prop_desc(of, dstr):
        if type(of) not in [list, np.ndarray]:
            if of.descriptor != '' or type(of.descriptor) != str:
                desc = str(of.descriptor)
            elif dstr != '':
                desc = str(dstr)
            else:
                desc = 'of'
            return desc

        desc = ['' for i in range(np.size(of))]
        for i in range(np.size(of)):
            if of[i].descriptor != '' or type(of[i].descriptor) != str:
                desc[i] = str(of[i].descriptor)
            elif dstr != '':
                desc[i] = str(dstr) + str(string_dim(of, i, 'vector'))
            else:
                desc[i] = 'of' + str(string_dim(of, i, 'vector'))

        desc = allempty(desc)
        return desc

    @staticmethod
    def prop_lower(of):
        lower = []
        return lower

    @staticmethod
    def prop_upper(of):
        upper = []
        return upper

    @staticmethod
    def prop_target(of):
        target = []
        return target

    @staticmethod
    def prop_weight(of):
        if type(of) not in [list, np.ndarray]:
            return of.weight

        weight = np.zeros(np.shape(of))
        for i in range(np.size(of)):
            weight[i] = of[i].weight

        weight = allequal(weight, 1.)
        return weight

    @staticmethod
    def prop_stype(of):
        if type(of) not in [list, np.ndarray]:
            return of.scale_type

        stype = ['' for i in range(np.size(of))]
        for i in range(np.size(of)):
            stype[i] = str(of[i].scale_type)

        stype = allequal(stype, 'none')
        return stype

    @staticmethod
    def prop_scale(of):
        if type(of) not in [list, np.ndarray]:
            return of.scale

        scale = np.zeros(np.shape(of))
        for i in range(np.size(of)):
            scale[i] = of[i].scale

        scale = allequal(scale, 1.)
        return scale

    @staticmethod
    def dakota_write(fidi, dresp, rdesc):
        # coloft only the variables of the appropriate class
        of = [struc_class(i, 'objective_functions', 'of') for i in dresp]

        # write constraints
        rdesc = rlist_write(fidi, 'objective_functions', 'objective_function', of, rdesc)
        return rdesc

    @staticmethod
    def dakota_rlev_write(fidi, dresp, params):
        return
