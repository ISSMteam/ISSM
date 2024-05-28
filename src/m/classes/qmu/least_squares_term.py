import numpy as np
from rlist_write import *
from MatlabArray import *


class least_squares_term(object):
    '''
  definition for the least_squares_term class.

  [lst] = least_squares_term.least_squares_term(args)
   lst = least_squares_term()

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
    def least_squares_term(*args):
        nargin = len(args)

        #create a default object
        if nargin == 0:
            return least_squares_term()

        #copy the object or create the object from the input
        else:
            if (nargin == 1) and isinstance(args[0], least_squares_term):
                lst = args[0]
            else:
                asizec = np.shape(*args[0:min(nargin, 4)])
                lst = [least_squares_term() for i in range(asizec[0]) for j in range(asizec[1])]

                for i in range(np.size(lst)):
                    if (np.size(args[0]) > 1):
                        lst[i].descriptor = args[0][i]
                    else:
                        lst[i].descriptor = str(args[0]) + string_dim(lst, i, 'vector')

                if (nargin >= 2):
                    for i in range(np.size(lst)):
                        if (np.size(args[1]) > 1):
                            lst[i].scale_type = args[1][i]
                        else:
                            lst[i].scale_type = str(args[1])

                if (nargin >= 3):
                    for i in range(np.size(lst)):
                        if (np.size(args[2]) > 1):
                            lst[i].scale = args[2][i]
                        else:
                            lst[i].scale = args[2]

                if (nargin >= 4):
                    for i in range(np.size(lst)):
                        if (np.size(args[3]) > 1):
                            lst[i].weight = args[3][i]
                        else:
                            lst[i].weight = args[3]

                if (nargin > 4):
                    print('WARNING: least_squares_term:extra_arg Extra arguments for object of class ' + str(type(lst)) + '.')

        return lst

    def __repr__(self):
        # display the object
        string = '\n'
        string += 'class "least_squares_term" object = \n'
        string += '    descriptor: ' + str(self.descriptor) + '\n'
        string += '    scale_type: ' + str(self.scale_type) + '\n'
        string += '         scale: ' + str(self.scale) + '\n'
        string += '        weight: ' + str(self.weight) + '\n'
        return string

    @staticmethod
    def prop_desc(lst, dstr):
        if type(lst) not in [list, np.ndarray]:
            lst = [lst]

        desc = ['' for i in range(np.size(lst))]
        for i in range(np.size(lst)):
            if lst[i].descriptor != '' or type(cdv[i].descriptor) != str:
                desc[i] = str(lst[i].descriptor)
            elif dstr != '':
                desc[i] = str(dstr) + str(string_dim(lst, i, 'vector'))
            else:
                desc[i] = 'lst' + str(string_dim(lst, i, 'vector'))

        desc = allempty(desc)
        return desc

    @staticmethod
    def prop_stype(lst):
        if type(lst) not in [list, np.ndarray]:
            return str(lst.scale_type)

        stype = ['' for i in range(np.size(lst))]
        for i in range(np.size(lst)):
            stype[i] = str(lst[i].scale_type)

        stype = allequal(stype, 'none')
        return stype

    @staticmethod
    def prop_scale(lst):
        if type(lst) not in [list, np.ndarray]:
            return lst.scale

        scale = np.zeros(np.size(lst))
        for i in range(np.size(lst)):
            scale[i] = lst[i].scale

        scale = allequal(scale, 1.)
        return scale

    @staticmethod
    def prop_weight(lst):
        if type(lst) not in [list, np.ndarray]:
            return lst.weight

        weight = np.zeros(np.size(lst))
        for i in range(np.size(lst)):
            weight[i] = lst[i].weight

        weight = allequal(weight, 1.)
        return weight

    @staticmethod
    def prop_lower(lst):
        lower = []
        return lower

    @staticmethod
    def prop_upper(lst):
        upper = []
        return upper

    @staticmethod
    def prop_target(lst):
        target = []
        return target

    @staticmethod
    def dakota_write(fidi, dresp, rdesc):
        #collect only the responses of the appropriate class
        lst = [struc_class(i, 'least_squares_term', 'lst') for i in dresp]

        #write responses
        rdesc = rlist_write(fidi, 'least_squares_terms', 'least_squares_term', lst, rdesc)
        return rdesc

    @staticmethod
    def dakota_rlev_write(fidi, dresp, params):
        return
