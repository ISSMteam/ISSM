import copy
import os.path
import numpy as np

from ContourToMesh import ContourToMesh
from ElementConnectivity import ElementConnectivity
import MatlabFuncs as m
from NodeConnectivity import NodeConnectivity


def contourenvelope(mh, *args):
    """CONTOURENVELOPE - build a set of segments enveloping a contour .exp

    Usage:
        segments = contourenvelope(mh, *args)

    Example:
        segments = contourenvelope(mh, 'Stream.exp')
        segments = contourenvelope(mh)
    """

    # Some checks
    nargs = len(args)

    if nargs > 1:
        print(contourenvelope.__doc__)
        raise Exception("contourenvelope error message: bad usage")

    if nargs == 1:
        flags = args[0]

        if isinstance(flags, str):
            file = flags
            if not os.path.exists(file):
                raise IOError("contourenvelope error message: file '%s' not found" % file)
            isfile = 1
        elif isinstance(flags, (bool, int, float)):
            #do nothing for now
            isfile = 0
        else:
            raise TypeError("contourenvelope error message: second argument should be a file or an elements flag")

    # Now, build the connectivity tables for this mesh
    # Computing connectivity
    if np.size(mh.vertexconnectivity, axis=0) != mh.numberofvertices and np.size(mh.vertexconnectivity, axis=0) != mh.numberofvertices2d:
        mh.vertexconnectivity = NodeConnectivity(mh.elements, mh.numberofvertices)
    if np.size(mh.elementconnectivity, axis=0) != mh.numberofelements and np.size(mh.elementconnectivity, axis=0) != mh.numberofelements2d:
        mh.elementconnectivity = ElementConnectivity(mh.elements, mh.vertexconnectivity)

    #get nodes inside profile
    elementconnectivity = copy.deepcopy(mh.elementconnectivity)
    if mh.dimension() == 2:
        elements = copy.deepcopy(mh.elements)
        x = copy.deepcopy(mh.x)
        y = copy.deepcopy(mh.y)
        numberofvertices = copy.deepcopy(mh.numberofvertices)
        numberofelements = copy.deepcopy(mh.numberofelements)
    else:
        elements = copy.deepcopy(mh.elements2d)
        x = copy.deepcopy(mh.x2d)
        y = copy.deepcopy(mh.y2d)
        numberofvertices = copy.deepcopy(mh.numberofvertices2d)
        numberofelements = copy.deepcopy(mh.numberofelements2d)

    if len(args) == 1:
        if isfile:
            # Get flag list of elements and nodes inside the contour
            nodein = ContourToMesh(elements, x, y, file, 'node', 1)
            elemin = (np.sum(nodein(elements), axis=1) == np.size(elements, axis=1))
            # Modify element connectivity
            elemout = np.nonzero(np.logical_not(elemin))[0]
            elementconnectivity[elemout, :] = 0
            elementconnectivity[np.nonzero(m.ismember(elementconnectivity, elemout + 1))] = 0
        else:
            # Get flag list of elements and nodes inside the contour
            nodein = np.zeros(numberofvertices)
            elemin = np.zeros(numberofelements)

            pos = np.nonzero(flags)
            elemin[pos] = 1
            nodein[elements[pos, :] - 1] = 1

            # Modify element connectivity
            elemout = np.nonzero(np.logical_not(elemin))[0]
            elementconnectivity[elemout, :] = 0
            elementconnectivity[np.nonzero(m.ismember(elementconnectivity, elemout + 1))] = 0

    # Find element on boundary
    # First: find elements on the boundary of the domain
    flag = copy.deepcopy(elementconnectivity)
    if len(args) == 1:
        flag[np.nonzero(flag)] = elemin[flag[np.nonzero(flag)]]
    elementonboundary = np.logical_and(np.prod(flag, axis=1) == 0, np.sum(flag, axis=1) > 0)

    # Find segments on boundary
    pos = np.nonzero(elementonboundary)[0]
    num_segments = np.size(pos)
    segments = np.zeros((num_segments * 3, 3), int)
    count = 0

    for el1 in pos:
        els2 = elementconnectivity[el1, np.nonzero(elementconnectivity[el1, :])[0]] - 1
        if np.size(els2) > 1:
            flag = np.intersect1d(np.intersect1d(elements[els2[0], :], elements[els2[1], :]), elements[el1, :])
            nods1 = elements[el1, :]
            nods1 = np.delete(nods1, np.nonzero(nods1 == flag))
            segments[count, :] = [nods1[0], nods1[1], el1 + 1]

            ord1 = np.nonzero(nods1[0] == elements[el1, :])[0][0]
            ord2 = np.nonzero(nods1[1] == elements[el1, :])[0][0]

    #swap segment nodes if necessary
            if ((ord1 == 0 and ord2 == 1) or (ord1 == 1 and ord2 == 2) or (ord1 == 2 and ord2 == 0)):
                temp = segments[count, 0]
                segments[count, 0] = segments[count, 1]
                segments[count, 1] = temp
            segments[count, 0:2] = np.flipud(segments[count, 0:2])
            count += 1
        else:
            nods1 = elements[el1, :]
            flag = np.setdiff1d(nods1, elements[els2, :])
            for j in range(0, 3):
                nods = np.delete(nods1, j)
                if np.any(m.ismember(flag, nods)):
                    segments[count, :] = [nods[0], nods[1], el1 + 1]
                    ord1 = np.nonzero(nods[0] == elements[el1, :])[0][0]
                    ord2 = np.nonzero(nods[1] == elements[el1, :])[0][0]
                    if ((ord1 == 0 and ord2 == 1) or (ord1 == 1 and ord2 == 2) or (ord1 == 2 and ord2 == 0)):
                        temp = segments[count, 0]
                        segments[count, 0] = segments[count, 1]
                        segments[count, 1] = temp
                    segments[count, 0:2] = np.flipud(segments[count, 0:2])
                    count += 1
    segments = segments[0:count, :]

    return segments
