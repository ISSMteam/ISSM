import numpy as np


def vector_write(fidi, sbeg, vec, nmax, cmax):
    '''
function to write a vector on multiple lines
'''
    if nmax is None:
        nmax = np.inf

    if cmax is None:
        cmax = np.inf

    # set up first iteration
    svec = []
    nitem = nmax
    lsvec = cmax

    # transpose vector from column - wise to row - wise
    vec = np.array(vec).conj().T

    # assemble each line, flushing when necessary
    for i in range(np.size(vec, 0)):

        # [[1], [1], [1]...] should be [1, 1, 1, ...]
        if type(vec[i]) in [list, np.ndarray] and len(vec[i]) == 1:
            sitem = str(vec[i][0])
        else:
            sitem = str(vec[i])

        nitem = nitem + 1
        lsvec = lsvec + 1 + len(sitem)

        if (nitem <= nmax) and (lsvec <= cmax):
            svec = str(svec) + ' ' + str(sitem)
        else:
            if len(svec) > 0:
                fidi.write(str(svec) + '\n')

            svec = str(sbeg) + str(sitem)
            nitem = 1
            lsvec = len(svec)

    # flush buffer at , if necessary
    if len(svec) > 0:
        fidi.write(str(svec) + '\n')

    return
