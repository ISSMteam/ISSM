import numpy as np

from MatlabFuncs import *


def importancefactors(md, variablename, responsename, partition):
    '''
    IMPORTANCEFACTORS - compute importance factors for a certain variable and response.

    Usage:
        factors = importancefactors(md, variablename, responsename)

    Example: factors = importancefactors(md, 'drag', 'max_vel')
    '''

    variablenamelength = len(variablename)

    #go through all response functions and find the one corresponding to the correct responsename
    responsefunctions = md.qmu.results.dresp_out
    found = -1
    for i in range(len(responsefunctions)):
        if responsefunctions[i].descriptor == responsename:
            found = i
            break
    if found < 0:
        raise RuntimeError('importancefactors error message: could not find correct response function')

    responsefunctions = responsefunctions[found]
    nfun = np.size(responsefunctions.var)

    #Now recover response to the correct design variable
    importancefactors = []
    count = 0
    for i in range(nfun):
        desvar = responsefunctions.var[i]
        if strncmpi(desvar, variablename, variablenamelength):
            importancefactors.append(responsefunctions.impfac[i])
            count = count + 1

    if count == 0:
        raise RuntimeError('importancefactors error message: either response does not exist, or importancefactors are empty')

    importancefactors = np.array(importancefactors)

    if count == 1:  #we have scalar
        factors = importancefactors
        return factors
    elif count == np.max(partition + 1):
        #distribute importance factor
        factors = importancefactors[(partition.conj().T).flatten().astype(int)]
    #md.qmu.partition was created to index "c" style
    else:
        #distribute importance factor
        factors = importancefactors[(partition.conj().T).flatten().astype(int)]
    #md.qmu.partition was created to index "c" style

    #weight importancefactors by area
    #if numel(factors) == md.mesh.numberofvertices,
    #  #get areas for each vertex.
    #    aire = GetAreas(md.mesh.elements, md.mesh.x, md.mesh.y)
    #    num_elements_by_node = md.nodeconnectivity(:, )
    #    grid_aire = zeros(md.mesh.numberofvertices, 1)
    #    for i = 1:md.mesh.numberofvertices,
    #        for j = 1:num_elements_by_node(i),
    #            grid_aire(i)=grid_aire(i) + aire(md.nodeconnectivity(i, j))
    #
    #
    #    factors = factors. / grid_aire
    #
    return factors
