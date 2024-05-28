import numpy as np
#move this later
from helpers import *

from vector_write import *
from response_function import *


def rlist_write(fidi, cstring, cstring2, dresp, rdesc):
    '''
    function to write response list
    '''
    if dresp is None:
        return

    func = eval(cstring)

    # if type(dresp) not in [list, np.ndarray]:
    #     dresp = [dresp]

    fnames = fieldnames(dresp)

    nresp = 0
    pstype = []
    pscale = []
    pweight = []
    plower = []
    pupper = []
    ptarget = []

    for i in range(len(fnames)):
        field = getattr(dresp, format(fnames[i]))
        nresp = nresp + len(field)
        pstype.extend(func.prop_stype(field))
        pscale.extend(func.prop_scale(field))
        pweight.extend(func.prop_weight(field))
        plower.extend(func.prop_lower(field))
        pupper.extend(func.prop_upper(field))
        ptarget.extend(func.prop_target(field))
        rdesc.extend(func.prop_desc(field, field[0].descriptor))

    # write responses
    print('  Writing ' + str(nresp) + ' ' + cstring + ' responses.')

    if strcmp(cstring, 'calibration_terms') == 1:
        fidi.write('\t' + cstring + ' = ' + str(nresp) + '\n')

    else:
        fidi.write('\tnum_' + cstring + 's = ' + str(nresp) + '\n')

    if not isempty(pstype):
        fidi.write('\t  ' + cstring2 + '_scale_types =\n')
        vector_write(fidi, '\t    ', pstype, 6, 76)

    if not isempty(pscale):
        fidi.write('\t  ' + cstring2 + '_scales =\n')
        vector_write(fidi, '\t    ', pscale, 6, 76)

    if not isempty(pweight):
        if cstring2 == 'objective_function':
            fidi.write('\t  multi_objective_weights =\n')
            vector_write(fidi, '\t    ', pweight, 6, 76)
        elif cstring2 == 'least_squares_term':
            fidi.write('\t  least_squares_weights =\n')
            vector_write(fidi, '\t    ', pweight, 6, 76)

    if not isempty(plower):
        fidi.write('\t  ' + cstring2 + '_lower_bounds =\n')
        vector_write(fidi, '\t    ', plower, 6, 76)

    if not isempty(pupper):
        fidi.write('\t  ' + cstring2 + '_upper_bounds =\n')
        vector_write(fidi, '\t    ', pupper, 6, 76)

    if not isempty(ptarget):
        fidi.write('\t  ' + cstring2 + '_targets =\n')
        vector_write(fidi, '\t    ', ptarget, 6, 76)

    # because qmu in files need '' for strings
    for i in range(len(rdesc)):
        if type(rdesc[i]) in [list, np.ndarray]:
            for j in range(len(rdesc[i])):
                rdesc[i][j] = "'" + rdesc[i][j] + "'"
        else:
            rdesc[i] = "'" + rdesc[i] + "'"

    return rdesc
