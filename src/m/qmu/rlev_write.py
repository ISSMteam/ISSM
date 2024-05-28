import numpy as np
#move this later
from helpers import *
from vector_write import *
from param_write import *
#import relevent qmu classes
from MatlabArray import *


def rlevi_write(fidi, ltype, levels):
    '''
  function to each type of response level
'''

    fidi.write('\t  num_' + str(ltype) + ' = ')
    levels = np.array(levels)

    if len(levels) > 0 and type(levels[0]) in [list, np.ndarray]:
        for i in range(len(levels)):
            fidi.write(' ' + str(len(levels[i])))
    else:
        fidi.write(' ' + str(len(levels)))

    fidi.write('\n')
    fidi.write('\t  ' + str(ltype) + ' =\n')

    # check if we have a vector of vectors, or just 1 vector
    if np.size(levels) > 0 and type(levels[0]) in [list, np.ndarray]:
        for i in range(len(levels)):
            if len(levels[i]) != 0:
                vector_write(fidi, '\t    ', levels[i], 8, 76)
    else:
        vector_write(fidi, '\t    ', levels, 8, 76)

    return


def rlev_write(fidi, dresp, cstring, params):
    '''
  function to write response levels
'''
    from response_function import response_function

    if len(dresp) == 0 or len(fieldnames(dresp[0])) == 0:
        return

    if type(dresp) in [list, np.ndarray]:
        if len(dresp) > 0 and type(dresp[0]) == struct:
            func = eval(cstring)
        else:
            func = type(dresp[0])
    elif type(dresp) == struct:
        # type is defined within struct's contents
        func = None
        dresp = [dresp]
    else:
        func = type(dresp)
        dresp = [dresp]

    # put responses into lists for writing

    nresp = 0
    respl = []
    probl = []
    rell = []
    grell = []

    # assume all fields in dvar[0:n] are consistent (ex. all are normal_uncertain)
    #   which will always be true since this is called per field
    fnames = fieldnames(dresp[0])
    for j in range(len(dresp)):
        for i in range(np.size(fnames)):
            if func is None:
                func = type(vars(dresp[j])[fnames[i]])

            nresp += 1
            [respli, probli, relli, grelli] = func.prop_levels([vars(dresp[j])[fnames[i]]])
            respl.extend(respli)
            probl.extend(probli)
            rell.extend(relli)
            grell.extend(grelli)

    # write response levels
    respl = allempty(respl)
    probl = allempty(probl)
    rell = allempty(rell)
    grell = allempty(grell)

    param_write(fidi, '\t  ', 'distribution', ' ', '\n', params)
    if len(respl) != 0:
        rlevi_write(fidi, 'response_levels', respl)
        param_write(fidi, '\t  ', 'compute', ' ', '\n', params)

    if len(probl) != 0:
        rlevi_write(fidi, 'probability_levels', probl)

    if len(rell) != 0:
        rlevi_write(fidi, 'reliability_levels', rell)

    if len(grell) != 0:
        rlevi_write(fidi, 'gen_reliability_levels', grell)

    return
