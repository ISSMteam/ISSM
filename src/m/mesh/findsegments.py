import numpy as np

from ElementConnectivity import ElementConnectivity
from helpers import struct
from intersect import intersect
from pairoptions import pairoptions

def findsegments(md, *args):  # {{{
    """FINDSEGMENTS - build segments model field

    Usage:
        segments = findsegments(md, args)

    Optional inputs:
        'mesh.elementconnectivity'
    """

    # Get options
    options = pairoptions(*args)

    # Get connectivity
    mesh = struct()
    mesh.elementconnectivity = options.getfieldvalue('mesh.elementconnectivity', md.mesh.elementconnectivity)

    # Now, build the connectivity tables for this mesh if not correctly done
    if md.mesh.elementconnectivity.shape[0] != md.mesh.numberofelements:
        if options.exist('mesh.elementconnectivity'):
            raise Exception('\'mesh.elementconnectivity\' option does not have the right size.')
        else:
            mesh.elementconnectivity = ElementConnectivity(md.mesh.elements, md.mesh.vertexconnectivity)

    # Recreate the segments
    elementonboundary = np.zeros((md.mesh.numberofelements, ))
    elementonboundary[np.where(mesh.elementconnectivity[:, 2] == 0)[0]] = 1
    pos = np.nonzero(elementonboundary)[0]
    num_segments = len(pos)
    segments = np.zeros((num_segments, 3)).astype(int)
    count = 0

    # Loop over the segments
    for i in range(num_segments):
        # Get current element on boundary
        el1 = pos[i]

        # Get elements connected to 'el1'
        els2 = mesh.elementconnectivity[el1, np.nonzero(mesh.elementconnectivity[el1, :])[0]]

        # Get nodes of 'el1'
        nods1 = md.mesh.elements[el1, :]

        # 'el1' is connected to 2 other elements
        if len(els2) > 1:

            # Find the common vertices to the two elements connected to 'el1' (1 or 2)
            flag = intersect(md.mesh.elements[els2[0] - 1, :], md.mesh.elements[els2[1] - 1, :])[0] # NOTE: Throwing away second- and third- position values returned from call

            # Get the vertices on the boundary and build segment
            if hasattr(np, 'isin'): #Numpy 2017+
                tmp = np.isin(nods1, flag, assume_unique=True)
            else: #For backward compatibility
                tmp = np.in1d(nods1, flag, assume_unique=True)
            nods1 = np.delete(nods1, np.where())
            segments[count, :] = np.append(nods1, el1 + 1)

            # Swap segment nodes if necessary
            ord1 = np.where(nods1[0] == md.mesh.elements[el1, :])[0][0]
            ord2 = np.where(nods1[1] == md.mesh.elements[el1, :])[0][0]

            if ((ord1 == 0 and ord2 == 1) or (ord1 == 1 and ord2 == 2) or (ord1 == 2 and ord2 == 0)):
                temp = segments[count, 0]
                segments[count, 0] = segments[count, 1]
                segments[count, 1] = temp

            segments[count, 0:2] = np.flip(segments[count, 0:2]) # NOTE: Upper bound of index range is non-inclusive
            count = count + 1
        # 'el1' is connected to only one element
        else:
            # NOTE: This block is untested as it does not get touched by test2004 (remove this note once it has been tested)

            # Find the vertex that 'el1' does not share with 'els2'
            flag = np.setdiff1d(nods, md.mesh.elements[els2, :])

            for j in range(3):
                nods = nods1
                nods = np.delete(nods, j)
                if hasattr(np, 'isin'): #Numpy 2017+
                    tmp = np.isin(flag, nods)
                else: #For backward compatibility
                    tmp = np.in1d(flag, nods)
                if np.any(tmp):
                    segments[count, :] = np.append(nods, el1 + 1)

                    # Swap segment nodes if necessary
                    ord1 = np.where(nods1[0] == md.mesh.elements[el1, :])[0][0]
                    ord2 = np.where(nods1[1] == md.mesh.elements[el1, :])[0][0]

                    if ((ord1 == 0 and ord2 == 1) or (ord1 == 1 and ord2 == 2) or (ord1 == 2 and ord2 == 0)):
                        temp = segments[count, 0]
                        segments[count, 0] = segments[count, 1]
                        segments[count, 1] = temp

                    segments[count, 0:2] = np.flip(segments[count, 0:2]) # NOTE: Upper bound of index range is non-inclusive
                    count = count + 1

    return segments
# }}}
