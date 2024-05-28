#move this later:
from helpers import *


def param_write(fidi, sbeg, pname, smid, s, params):
    '''
function to write a parameter
'''
    if not isfield(params, pname):
        print('WARNING: param_write:param_not_found: Parameter {} not found in structure.'.format(pname))
        return

    params_pname = vars(params)[pname]

    if type(params_pname) == bool and (not params_pname):
        return

    if type(params_pname) == bool:
        fidi.write(sbeg + str(pname) + s)

    elif type(params_pname) in [str]:
        fidi.write(sbeg + pname + smid + params_pname + s)

    elif type(params_pname) in [int, float]:
        fidi.write(sbeg + str(pname) + smid + str(params_pname) + s)
