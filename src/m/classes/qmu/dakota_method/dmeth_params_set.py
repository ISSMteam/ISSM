from helpers import *
from dakota_method import *


def dmeth_params_set(dm, *args):
    #
    #  set parameters of a dakota_method object.
    #
    #  dm = dmeth_params_set(dm, *args)
    #

    if not isinstance(dm, dakota_method):
        raise RuntimeError('Provided object is a \'' + str(type(dm)) + '\' class object, not \'dakota_method\'')

    #  loop through each parameter field in the input list
    for i in range(0, len(args), 2):
        if isfield(dm.params, args[i]):
            #vars(dresp)[fnames[i]]
            exec(('dm.params.%s = args[i + 1]') % (args[i]))
    #vars(dm.params)[args[i]] = args[i + 1]
        else:
            print('WARNING: dmeth_params_set:unknown_param No parameter \'' + str(args[i]) + '\' for dakota_method \'' + str(dm.method) + '\'.')

    return dm
