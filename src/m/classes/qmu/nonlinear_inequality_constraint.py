import numpy as np
from rlist_write import *
from MatlabArray import *


class nonlinear_inequality_constraint:
    '''
  constructor for the nonlinear_inequality_constraint class.

  [nic] = nonlinear_inequality_constraint.nonlinear_inequality_constraint(args)
   nic = nonlinear_inequality_constraint()

  where the required args are:
    descriptor    (char, description, '')
    lower         (double, lower bound, -np.inf)
    upper         (double, upper bound, 0.)
  and the optional args and defaults are:
    scale_type    (char, scaling type, 'none')
    scale         (double, scaling factor, 1.)

  note that zero arguments constructs a default instance, one
  argument of the class copies the instance, and three or more
  arguments constructs a new instance from the arguments.
'''
    def __init__(self):
        self.descriptor = ''
        self.lower = -np.inf
        self.upper = 0.
        self.scale_type = 'none'
        self.scale = 1.

    @staticmethod
    def nonlinear_inequality_constraint(*args):
        nargin = len(args)

        # create a default object
        if nargin == 0:
            return nonlinear_inequality_constraint()

        # copy the object
        if nargin == 1:
            if isinstance(args[0], nonlinear_inequality_constraint):
                nic = args[0]
            else:
                raise RuntimeError('Object is a ' + str(type(args[0])) + ' class object, not "nonlinear_inequality_constraint".')

        # not enough arguments
        if nargin == 2:
            raise RuntimeError('Construction of nonlinear_inequality_constraint class object requires at least 3 inputs.')

        # create the object from the input
        else:
            asizec = array_size(*args[0:min(nargin, 3)])
            nic = [nonlinear_inequality_constraint() for i in range(asizec[0]) for j in range(asizec[1])]

            for i in range(np.size(nic)):
                if (np.size(args[0]) > 1):
                    nic[i].descriptor = args[0][i]
                else:
                    nic[i].descriptor = str(args[0]) + string_dim(nic, i, 'vector')
                if (np.size(args[1]) > 1):
                    nic[i].lower = args[1][i]
                else:
                    nic[i].lower = args[1]

                if (np.size(args[2]) > 1):
                    nic[i].upper = args[2][i]
                else:
                    nic[i].upper = args[2]

            if (nargin >= 4):
                for i in range(np.size(nic)):
                    if (np.size(args[3]) > 1):
                        nic[i].scale_type = args[3][i]
                    else:
                        nic[i].scale_type = str(args[3])

            if (nargin >= 5):
                for i in range(np.size(nic)):
                    if (np.size(args[4]) > 1):
                        nic[i].upper = args[4][i]
                    else:
                        nic[i].upper = args[4]

            if (nargin > 5):
                print('WARNING: nonlinear_inequality_constraint:extra_arg: Extra arguments for object of class ' + str(type(nic)) + '.')

        return nic

    def __repr__(self):
        # display the object
        string = '\n'
        string += 'class "nonlinear_inequality_constraint" object = \n'
        string += '    descriptor: ' + str(self.descriptor) + '\n'
        string += '         lower: ' + str(self.lower) + '\n'
        string += '         upper: ' + str(self.upper) + '\n'
        string += '    scale_type: ' + str(self.scale_type) + '\n'
        string += '         scale: ' + str(self.scale) + '\n'

        return string

    @staticmethod
    def prop_desc(nic, dstr):
        if type(nic) not in [list, np.ndarray]:
            if nic.descriptor != '' or type(nic.descriptor) != str:
                desc = str(nic.descriptor)
            elif dstr != '':
                desc = str(dstr)
            else:
                desc = 'nic'
            return desc

        desc = ['' for i in range(np.size(nic))]
        for i in range(np.size(nic)):
            if nic[i].descriptor != '' or type(nic[i].descriptor) != str:
                desc[i] = str(nic[i].descriptor)
            elif dstr != '':
                desc[i] = str(dstr) + str(string_dim(nic, i, 'vector'))
            else:
                desc[i] = 'nic' + str(string_dim(nic, i, 'vector'))

        desc = allempty(desc)

        return desc

    @staticmethod
    def prop_stype(nic):
        if type(nic) not in [list, np.ndarray]:
            return nic.scale_type

        stype = ['' for i in range(np.size(nic))]
        for i in range(np.size(nic)):
            stype[i] = str(nic[i].scale_type)

        stype = allequal(stype, 'none')

        return stype

    @staticmethod
    def prop_scale(nic):
        if type(nic) not in [list, np.ndarray]:
            return nic.scale

        scale = np.zeros(np.shape(nic))
        for i in range(np.size(nic)):
            scale[i] = nic[i].scale

        scale = allequal(scale, 1.)

        return scale

    @staticmethod
    def prop_weight(nic):
        weight = []
        return weight

    @staticmethod
    def prop_lower(nic):
        if type(nic) not in [list, np.ndarray]:
            return nic.lower

        lower = np.zeros(np.shape(nic))
        for i in range(np.size(nic)):
            lower[i] = nic[i].lower

        lower = allequal(lower, -np.inf)

        return lower

    @staticmethod
    def prop_upper(nic):
        if type(nic) not in [list, np.ndarray]:
            return nic.upper

        upper = np.zeros(np.shape(nic))
        for i in range(np.size(nic)):
            upper[i] = nic[i].upper

        upper = allequal(upper, 0.)

        return upper

    @staticmethod
    def prop_target(nic):
        target = []
        return target

    @staticmethod
    def dakota_write(fidi, dresp, rdesc):
        # collect only the variables of the appropriate class
        nic = [struc_class(i, 'nonlinear_inequality_constraint', 'nic') for i in dresp]

        # write constraints
        rdesc = rlist_write(fidi, 'nonlinear_inequality_constrants', 'nonlinear_inequality', nic, rdesc)
        return rdesc

    @staticmethod
    def dakota_rlev_write(fidi, dresp, params):
        return
