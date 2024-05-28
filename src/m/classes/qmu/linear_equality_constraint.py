import numpy as np
from lclist_write import *
from MatlabArray import *


class linear_equality_constraint:
    '''
  constructor for the linear_equality_constraint class.

  [lec] = linear_equality_constraint.linear_equality_constraint(args)
   lec = linear_equality_constraint()

  where the required args are:
    matrix        (double row, variable coefficients, float('NaN'))
    target        (double vector, target values, 0.)
  and the optional args and defaults are:
    scale_type    (char, scaling type, 'none')
    scale         (double, scaling factor, 1.)

  note that zero arguments constructs a default instance, one
  argument of the class copies the instance, and two or more
  arguments constructs a new instance from the arguments.
'''
    def __init__(self):
        self.matrix = np.array([[float('NaN')]])
        self.target = 0.
        self.scale_type = 'none'
        self.scale = 1.

    @staticmethod
    def linear_equality_constraint(*args):
        nargin = len(args)

        #  create a default object
        if nargin == 0:
            return linear_equality_constraint()

        #  copy the object
        elif nargin == 1:
            if isinstance(args[0], linear_equality_constraint):
                lec = args[0]
            else:
                raise RuntimeError('Object is a ' + str(type(args[0])) + ' class object, not "linear_equality_constraint"')

        #  create the object from the input
        else:
            if (np.shape(args[0], 1) == array_np.size(args[1:min(nargin, 4)]) or np.shape(args[0], 1) == 1):
                asizec = np.shape(args[1:min(nargin, 4)])
            elif (array_np.size(args[1:min(nargin, 4)]) == 1):
                asizec = [np.shape(args[0], 1), 1]
            else:
                raise RuntimeError('Matrix for object of class ' + str(type(lec)) + ' has inconsistent number of rows.')

            lec = [linear_equality_constraint() for i in range(asizec[0]) for j in range(asizec[1])]

            for i in range(np.size(lec)):
                if (np.shape(args[0])[0] > 1):
                    lec[i].matrix = args[0][i, :]
                else:
                    lec[i].matrix = args[0]

            if (nargin >= 2):
                for i in range(np.size(lec)):
                    if (np.size(args[1]) > 1):
                        lec[i].target = args[1][i]
                    else:
                        lec[i].target = args[1]

            if (nargin >= 3):
                for i in range(np.size(lec)):
                    if (np.size(args[2]) > 1):
                        lec[i].scale_type = args[2][i]
                    else:
                        lec[i].scale_type = str(args[2])

            if (nargin >= 4):
                for i in range(np.size(lec)):
                    if (np.size(args[3]) > 1):
                        lec[i].scale = args[3][i]
                    else:
                        lec[i].scale = args[3]

            if (nargin > 4):
                print('WARNING: linear_equality_constraint:extra_arg: Extra arguments for object of class ' + str(type(lec)) + '.')

        return lec

    def __repr__(self):
        # display the object
        string = '\n'
        string += 'class "linear_equality_constraint" object = \n'
        string += '        matrix: ' + str(self.matrix) + '\n'
        string += '        target: ' + str(self.target) + '\n'
        string += '    scale_type: ' + str(self.scale_type) + '\n'
        string += '         scale: ' + str(self.scale) + '\n'

        return string

    @staticmethod
    def prop_matrix(lec):
        if type(lec) not in [list, np.ndarray]:
            return lec.matrix

        matrix = np.zeros(np.size(lec))
        for i in range(np.size(lec)):
            matrix[i, 0:np.shape(lec[i].matrix)[1]] = lec[i].matrix[0, :]

        return matrix

    @staticmethod
    def prop_lower(lec):
        lower = []
        return lower

    @staticmethod
    def prop_upper(lec):
        upper = []
        return upper

    @staticmethod
    def prop_target(lec):
        if type(lec) not in [list, np.ndarray]:
            return lec.target

        target = np.zeros(np.shape(lec))
        for i in range(np.size(lec)):
            target[i] = lec[i].target

        target = allequal(target, 0.)

        return target

    @staticmethod
    def prop_stype(lec):
        if type(lec) not in [list, np.ndarray]:
            return lec.scale_type

        stype = ['' for i in range(np.size(lec))]
        for i in range(np.size(lec)):
            stype[i] = str(lec[i].scale_type)

        stype = allequal(stype, 'none')

        return stype

    @staticmethod
    def prop_scale(lec):
        if type(lec) not in [list, np.ndarray]:
            return lec.scale

        scale = np.zeros(np.shape(lec))
        for i in range(np.size(lec)):
            scale[i] = lec[i].scale

        scale = allequal(scale, 1.)

        return scale

    @staticmethod
    def dakota_write(fidi, dvar):
        # collect only the variables of the appropriate class
        lec = [struc_type(i, 'linear_equality_constraint', 'lec') for i in dvar]

        # write constraints
        lclist_write(fidi, 'linear_equality_constraints', 'linear_equality', lec)
