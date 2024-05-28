import numpy as np
from lclist_write import *
from MatlabArray import *


class linear_inequality_constraint:
    '''
  constructor for the linear_inequality_constraint class.

  [lic] = linear_inequality_constraint.linear_inequality_constraint(args)
   lic = linear_inequality_constraint()

  where the required args are:
    matrix        (double row, variable coefficients, float('NaN'))
    lower         (double vector, lower bounds, -np.inf)
    upper         (double vector, upper bounds, 0.)
  and the optional args and defaults are:
    scale_type    (char, scaling type, 'none')
    scale         (double, scaling factor, 1.)

  note that zero arguments constructs a default instance, one
  argument of the class copies the instance, and three or more
  arguments constructs a new instance from the arguments.
'''
    def __init__(self):
        self.matrix = np.array([[float('NaN')]])
        self.lower = -np.inf
        self.upper = 0.
        self.scale_type = 'none'
        self.scale = 1.

    @staticmethod
    def linear_inequality_constraint(*args):
        nargin = len(args)

        # create a default object
        if nargin == 0:
            return linear_inequality_constraint()

        # copy the object
        if nargin == 1:
            if isinstance(args[0], linear_inequality_constraint):
                lic = args[0]
            else:
                raise RuntimeError('Object is a ' + str(type(args[0])) + ' class object, not "linear_inequality_constraint".')

        # not enough arguments
        if nargin == 2:
            raise RuntimeError('Construction of linear_inequality_constraint class object requires at least 3 inputs.')

        # create the object from the input
        else:
            if (np.shape(args[0], 1) == array_numel(args[1:min(nargin, 5)]) or np.shape(args[0], 1) == 1):
                asizec = array_size(args[1:min(nargin, 5)])
            elif (array_numel(args[1:min(nargin, 5)]) == 1):
                asizec = [array_size(args[0], 1), 1]
            else:
                raise RuntimeError('Matrix for object of class ' + str(type(lic)) + ' has inconsistent number of rows.')

            lic = [linear_inequality_constraint() for i in range(asizec[0]) for j in range(asizec[1])]

            for i in range(np.size(lic)):
                if (np.shape(args[0], 1) > 1):
                    lic[i].matrix = args[0][i, :]
                else:
                    lic[i].matrix = args[0]

            if (nargin >= 2):
                for i in range(np.size(lic)):
                    if (np.size(args[1]) > 1):
                        lic[i].lower = args[1][i]
                    else:
                        lic[i].lower = args[1]

            if (nargin >= 3):
                for i in range(np.size(lic)):
                    if (np.size(args[2]) > 1):
                        lic[i].upper = args[2][i]
                    else:
                        lic[i].upper = args[2]

            if (nargin >= 4):
                for i in range(np.size(lic)):
                    if (np.size(args[3]) > 1):
                        lic[i].scale_type = args[3][i]
                    else:
                        lic[i].scale_type = str(args[3])

            if (nargin >= 5):
                for i in range(np.size(lic)):
                    if (np.size(args[4]) > 1):
                        lic[i].scale = args[4][i]
                    else:
                        lic[i].scale = args[4]

            if (nargin > 5):
                print('WARNING: linear_inequality_constraint:extra_arg: Extra arguments for object of class ' + str(type(lic)) + '.')

        return lic

    def __repr__(self):
        # display the object
        string = '\n'
        string += 'class "linear_inequality_constraint" object = \n'
        string += '        matrix: ' + str(string_vec(self.matrix)) + '\n'
        string += '         lower: ' + str(self.lower) + '\n'
        string += '         upper: ' + str(self.upper) + '\n'
        string += '    scale_type: ' + str(self.scale_type) + '\n'
        string += '         scale: ' + str(self.scale) + '\n'

        return string

    @staticmethod
    def prop_matrix(lic):
        if type(lic) not in [list, np.ndarray]:
            return lic.matrix

        matrix = np.zeros(np.size(lic))
        for i in range(np.size(lic)):
            matrix[i, 0:np.shape(lic[i].matrix)[1]] = lic[i].matrix[0, :]

        return matrix

    @staticmethod
    def prop_lower(lic):
        if type(lic) not in [list, np.ndarray]:
            return lic.lower

        lower = np.zeros(np.shape(lic))
        for i in range(np.size(lic)):
            lower[i] = lic[i].lower

        lower = allequal(lower, -np.inf)

        return lower

    @staticmethod
    def prop_upper(lic):
        if type(lic) not in [list, np.ndarray]:
            return lic.upper

        upper = np.zeros(np.shape(lic))
        for i in range(np.size(lic)):
            upper[i] = lic[i].upper

        upper = allequal(upper, 0.)

        return upper

    @staticmethod
    def prop_target(lic):
        target = []
        return target

    @staticmethod
    def prop_stype(lic):
        if type(lic) not in [list, np.ndarray]:
            return lic.scale_type

        stype = ['' for i in range(np.size(lic))]
        for i in range(np.size(lic)):
            stype[i] = str(lic[i].scale_type)

        stype = allequal(stype, 'none')

        return stype

    @staticmethod
    def prop_scale(lic):
        if type(lic) not in [list, np.ndarray]:
            return lic.scale

        scale = np.zeros(np.shape(lic))
        for i in range(np.size(lic)):
            scale[i] = lic[i].scale

        scale = allequal(scale, 1.)

        return scale

    @staticmethod
    def dakota_write(fidi, dvar):
        # collect only the variables of the appropriate class
        lic = [struc_class(i, 'linear_inequality_constraint', 'lic') for i in dvar]

        # write constraints
        lclist_write(fidi, 'linear_inequality_constraints', 'linear_inequality', lic)
