import numpy as np
#move this later
from helpers import *
from vector_write import *
from uniform_uncertain import *
from normal_uncertain import *


def check(a, l, p):
    '''in the event that a and b are equal, return a
    in the event that a and b are not equal, return their concatenation

    This is used for when both the input dvar and the 'cstring' variables have non - 1 length
    '''
    if np.size(a) == l:
        if p == 0:
            if type(a) in [list, np.ndarray]:
                return a
            else:
                return [a]
        else:
            return []
    elif np.size(a) == 1:
        if type(a) in [list, np.ndarray]:
            return a
        else:
            return [a]
    elif np.size(a) == 0:
        return []
    else:
        raise RuntimeError('ERROR vlist_write: input field had size ' + str(np.size(a)) + '; must have size of 0, 1, or match size of provided dvar (' + str(l) + ').')

    return


def vlist_write(fidi, cstring, cstring2, dvar):
    '''
    function to write variable list
    '''
    if dvar is None:
        return

    func = eval(cstring)

    # if type(dvar) not in [list, np.ndarray]:
    #     dvar = [dvar]
    
    fnames = fieldnames(dvar)
    
    # put variables into lists for writing
    nvar = 0
    pinitpt = []
    plower = []
    pupper = []
    pmean = []
    pstddev = []
    pinitst = []
    pstype = []
    pscale = []
    ppairs_per_variable = []
    pabscissas = []
    pcounts = []
    pdesc = []

    for i in range(len(fnames)):
        field = getattr(dvar, format(fnames[i]))
        nvar = nvar + len(field)
        pinitpt.extend(func.prop_initpt(field))
        plower.extend(func.prop_lower(field))
        pupper.extend(func.prop_upper(field))
        pmean.extend(func.prop_mean(field))
        pstddev.extend(func.prop_stddev(field))
        pinitst.extend(func.prop_initst(field))
        pstype.extend(func.prop_stype(field))
        pscale.extend(func.prop_scale(field))
        ppairs_per_variable.extend(func.prop_pairs_per_variable(field))
        pabscissas.extend(func.prop_abscissas(field))
        pcounts.extend(func.prop_counts(field))
        pdesc.extend(func.prop_desc(field, field[0].descriptor))

    pinitpt = allempty(pinitpt)
    plower = allempty(plower)
    pupper = allempty(pupper)
    pmean = allempty(pmean)
    pstddev = allempty(pstddev)
    pinitst = allempty(pinitst)
    pstype = allempty(pstype)
    pscale = allempty(pscale)
    ppairs_per_variable = allempty(ppairs_per_variable)
    pabscissas = allempty(pabscissas)
    pcounts = allempty(pcounts)
    pdesc = allempty(pdesc)

    # write variables
    print('  Writing ' + str(nvar) + ' ' + cstring + ' variables.')

    fidi.write('\t' + cstring + ' = ' + str(nvar) + '\n')
    if not isempty(pinitpt):
        fidi.write('\t  ' + cstring2 + '_initial_point =\n')
        vector_write(fidi, '\t    ', pinitpt, 6, 76)

    if not isempty(plower):
        fidi.write('\t  ' + cstring2 + '_lower_bounds =\n')
        vector_write(fidi, '\t    ', plower, 6, 76)

    if not isempty(pupper):
        fidi.write('\t  ' + cstring2 + '_upper_bounds =\n')
        vector_write(fidi, '\t    ', pupper, 6, 76)

    if not isempty(pmean):
        fidi.write('\t  ' + cstring2 + '_means =\n')
        vector_write(fidi, '\t    ', pmean, 6, 76)

    if not isempty(pstddev):
        fidi.write('\t  ' + cstring2 + '_std_deviations =\n')
        vector_write(fidi, '\t    ', pstddev, 6, 76)

    if not isempty(pinitst):
        fidi.write('\t  ' + cstring2 + '_initial_state =\n')
        vector_write(fidi, '\t    ', pinitst, 6, 76)

    if not isempty(pstype):
        fidi.write('\t  ' + cstring2 + '_scale_types =\n')
        vector_write(fidi, '\t    ', pstype, 6, 76)

    if not isempty(pscale):
        fidi.write('\t  ' + cstring2 + '_scales =\n')
        vector_write(fidi, '\t    ', pscale, 6, 76)

    if not isempty(ppairs_per_variable):
        fidi.write('\t  ' + cstring2 + '_pairs_per_variable =\n')
        vector_write(fidi, '\t    ', ppairs_per_variable, 6, 76)

    if not isempty(pabscissas):
        fidi.write('\t  ' + cstring2 + '_abscissas =\n')
        vector_write(fidi, '\t    ', pabscissas, 6, 76)

    if not isempty(pcounts):
        fidi.write('\t  ' + cstring2 + '_counts =\n')
        vector_write(fidi, '\t    ', pcounts, 6, 76)

    if not isempty(pdesc):
        fidi.write('\t  ' + 'descriptors =\n')
        vector_write(fidi, '\t    ', pdesc, 6, 76)

    return
