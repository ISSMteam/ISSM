import numpy as np
#move this later
from helpers import *

from vector_write import *

from linear_equality_constraint import *
from linear_inequality_constraint import *


def lclist_write(fidi, cstring, cstring2, dvar):
    '''
function to write linear constraint list
'''
    if dvar is None:
        return

    #  put linear constraints into lists for writing

    nvar = 0
    pmatrix = []
    plower = []
    pupper = []
    ptarget = []
    pstype = []
    pscale = []

    fnames = fieldnames(dvar)
    for i in range(np.size(fnames)):
        nvar = nvar + np.size(vars(dvar)[fnames[i]])
        pmatrix = [pmatrix, prop_matrix(vars(dvar)[fnames[i]])]
        plower = [plower, prop_lower(vars(dvar)[fnames[i]])]
        pupper = [pupper, prop_upper(vars(dvar)[fnames[i]])]
        ptarget = [ptarget, prop_target(vars(dvar)[fnames[i]])]
        pstype = [pstype, prop_stype(vars(dvar)[fnames[i]])]
        pscale = [pscale, prop_scale(vars(dvar)[fnames[i]])]

    #  write linear constraints
    print('  Writing ' + str(nvar) + ' ' + cstring + ' linear constraints.')

    if len(pmatrix) != 0:
        fidi.write('\t  ' + cstring2 + '_matrix =\n')
        vector_write(fidi, '\t    ', pmatrix, 6, 76)

    if len(plower) != 0:
        fidi.write('\t  ' + cstring2 + '_lower_bounds =\n')
        vector_write(fidi, '\t    ', plower, 6, 76)

    if len(pupper) != 0:
        fidi.write('\t  ' + cstring2 + '_upper_bounds =\n')
        vector_write(fidi, '\t    ', pupper, 6, 76)

    if len(ptarget) != 0:
        fidi.write('\t  ' + cstring2 + '_targets =\n')
        vector_write(fidi, '\t    ', ptarget, 6, 76)

    if len(pstype) != 0:
        fidi.write('\t  ' + cstring2 + '_scale_types =\n')
        vector_write(fidi, '\t    ', pstype, 6, 76)

    if len(pscale) != 0:
        fidi.write('\t  ' + cstring2 + '_scales =\n')
        vector_write(fidi, '\t    ', pscale, 6, 76)
