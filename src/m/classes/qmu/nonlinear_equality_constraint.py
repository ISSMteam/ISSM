import numpy as np
from rlist_write import *
from MatlabArray import *


class nonlinear_equality_constraint:
    '''
  constructor for the nonlinear_equality_constraint class.

  [nec] = nonlinear_equality_constraint.nonlinear_equality_constraint(args)
   nec = nonlinear_equality_constraint()

  where the required args are:
    descriptor    (char, description, '')
    target        (double vector, target values, 0.)
  and the optional args and defaults are:
    scale_type    (char, scaling type, 'none')
    scale         (double, scaling factor, 1.)

  note that zero arguments constructs a default instance, one
  argument of the class copies the instance, and two or more
  arguments constructs a new instance from the arguments.
'''
    def __init__(self):
        self.descriptor = ''
        self.target = 0.
        self.scale_type = 'none'
        self.scale = 1.

    @staticmethod
    def nonlinear_equality_constraint(*args):
        nargin = len(args)

    # create a default object
        if nargin == 0:
            return nonlinear_equality_constraint()

    # copy the object
        elif nargin == 1:
            if isinstance(args[0], nonlinear_equality_constraint):
                nec = args[0]
            else:
                raise RuntimeError('Object is a ' + str(type(args[0])) + ' class object, not "nonlinear_equality_constraint"')

    # create the object from the input
        else:
            asizec = array_size(*args[0:min(nargin, 4)])
            nec = [nonlinear_equality_constraint() for i in range(asizec[0]) for j in range(asizec[1])]

            for i in range(np.size(nec)):
                if (np.shape(args[0])[0] > 1):
                    nec[i].descriptor = args[0][i]
                else:
                    nec[i].descriptor = str(args[0])
                if (np.size(args[1]) > 1):
                    nec[i].target = args[1][i]
                else:
                    nec[i].target = args[1]

            if (nargin >= 3):
                for i in range(np.size(nec)):
                    if (np.size(args[2]) > 1):
                        nec[i].scale_type = args[2][i]
                    else:
                        nec[i].scale_type = str(args[2])

            if (nargin >= 4):
                for i in range(np.size(nec)):
                    if (np.size(args[3]) > 1):
                        nec[i].scale = args[3][i]
                    else:
                        nec[i].scale = args[3]

            if (nargin > 4):
                print('WARNING: nonlinear_equality_constraint:extra_arg: Extra arguments for object of class ' + str(type(nec)) + '.')

        return nec

    def __repr__(self):
        # display the object
        string = '\n'
        string += 'class "nonlinear_equality_constraint" object = \n'
        string += '    descriptor: ' + str(self.descriptor) + '\n'
        string += '        target: ' + str(self.target) + '\n'
        string += '    scale_type: ' + str(self.scale_type) + '\n'
        string += '         scale: ' + str(self.scale) + '\n'

        return string

    @staticmethod
    def prop_desc(nec, dstr):
        if type(nec) not in [list, np.ndarray]:
            if nec.descriptor != '' or type(nec.descriptor) != str:
                desc = str(nec.descriptor)
            elif dstr != '':
                desc = str(dstr)
            else:
                desc = 'nec'
            return desc

        desc = ['' for i in range(np.size(nec))]
        for i in range(np.size(nec)):
            if nec[i].descriptor != '' or type(nec[i].descriptor) != str:
                desc[i] = str(nec[i].descriptor)
            elif dstr != '':
                desc[i] = str(dstr) + str(string_dim(nec, i, 'vector'))
            else:
                desc[i] = 'nec' + str(string_dim(nec, i, 'vector'))

        desc = allempty(desc)

        return desc

    @staticmethod
    def prop_lower(nec):
        lower = []
        return lower

    @staticmethod
    def prop_upper(nec):
        upper = []
        return upper

    @staticmethod
    def prop_weight(nec):
        weight = []
        return weight

    @staticmethod
    def prop_target(nec):
        if type(nec) not in [list, np.ndarray]:
            return nec.target

        target = np.zeros(np.shape(nec))
        for i in range(np.size(nec)):
            target[i] = nec[i].target

        target = allequal(target, 0.)

        return target

    @staticmethod
    def prop_stype(nec):
        if type(nec) not in [list, np.ndarray]:
            return nec.scale_type

        stype = ['' for i in range(np.size(nec))]
        for i in range(np.size(nec)):
            stype[i] = str(nec[i].scale_type)

        stype = allequal(stype, 'none')

        return stype

    @staticmethod
    def prop_scale(nec):
        if type(nec) not in [list, np.ndarray]:
            return nec.scale

        scale = np.zeros(np.shape(nec))
        for i in range(np.size(nec)):
            scale[i] = nec[i].scale

        scale = allequal(scale, 1.)

        return scale

    @staticmethod
    def dakota_write(fidi, dresp, rdesc):
        #  colnect only the variables of the appropriate class
        nec = [struc_type(i, 'nonlinear_equality_constraint', 'nec') for i in dresp]

        #  write constraints
        rdesc = rlist_write(fidi, 'nonlinear_equality_constraints', 'nonlinear_equality', nec, rdesc)
        return rdesc

    @staticmethod
    def dakota_rlev_write(fidi, dresp, params):
        return
