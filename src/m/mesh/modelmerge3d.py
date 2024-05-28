from copy import copy, deepcopy

import numpy as np

from ElementConnectivity import ElementConnectivity
from findsegments import findsegments
from model import model
from NodeConnectivity import NodeConnectivity
from pairoptions import pairoptions


def modelmerge3d(md1, md2, *args):
    """MODELMERGE - Merge two models by merging their meshes.

    Usage:
        md = modelmerge(md1, md2)
    """

    # Process options
    options = pairoptions(*args)

    tolerance = options.getfieldvalue('tolerance', 1e-4)

    md = model()

    # First, copy md1 mesh into md.mesh to initialize, and additional classes
    md.mesh = deepcopy(md1.mesh)
    md.private = deepcopy(md1.private)

    # Some initialization
    elements1 = copy(md1.mesh.elements)
    x1 = copy(md1.mesh.x)
    y1 = copy(md1.mesh.y)
    z1 = copy(md1.mesh.z)
    nods1 = copy(md1.mesh.numberofvertices)
    nel1 = copy(md1.mesh.numberofelements)

    elements2 = copy(md2.mesh.elements)
    x2 = copy(md2.mesh.x)
    y2 = copy(md2.mesh.y)
    z2 = copy(md2.mesh.z)
    nods2 = copy(md2.mesh.numberofvertices)
    nel2 = copy(md2.mesh.numberofelements)

    # Offset elements2 by nods1
    elements2 = elements2 + nods1

    # Go into the vertices on boundary of mesh 1 and figure out which ones are common with mesh 2
    verticesonboundary = np.nonzero(md1.mesh.vertexonboundary)[0]

    # Do not display "RuntimeWarning: invalid value encountered in less" warning when comparing against np.nan values in 'x2', 'y2', and 'z2'
    #
    # TODO: Investigate if we should use some other value to represent null elements
    np.warnings.filterwarnings('ignore')

    for i in range(len(verticesonboundary)):
        node1 = verticesonboundary[i]
        xnode1 = x1[node1]
        ynode1 = y1[node1]
        znode1 = z1[node1]

        # Is there another node with these coordinates in mesh 2?
        ind = np.where(np.logical_and.reduce((~np.isnan(x2), ~np.isnan(y2), ~np.isnan(z2), ((x2 - xnode1) ** 2 + (y2 - ynode1) ** 2 + (z2 - znode1) ** 2) ** 0.5 < tolerance)))[0] # NOTE: math.sqrt cannot be applied element-wise to a list/numpy.array
        if len(ind) > 1:
            print('should reduce the tolerance, several vertices picked up!')
        if len(ind):
            x2[ind] = np.nan
            y2[ind] = np.nan
            z2[ind] = np.nan
            pos = np.where(elements2 == ((ind + 1) + nods1))
            elements2[pos] = (node1 + 1) # NOTE: 'elements2' is 2D array, so 'numpy.where' returns two lists of indices

    # Go through elements2 and drop counter on each vertex that is above the x2 and y2 vertices being dropped
    indices_nan = np.nonzero(np.isnan(x2).astype(int))[0]
    while indices_nan.size:
        # Use the index of the first instance of 'nan' value to remove that element from 'x2', 'y2', and 'z2'
        index_nan = indices_nan[0]
        pos = np.where(elements2 > ((index_nan + 1) + nods1))
        elements2[pos] = elements2[pos] - 1
        # TODO: Maybe set to None, then pop all None after this loop
        x2 = np.delete(x2, index_nan)
        y2 = np.delete(y2, index_nan)
        z2 = np.delete(z2, index_nan)

        # Check again in 'x2' for instances of 'nan'
        indices_nan = np.nonzero(np.isnan(x2).astype(int))[0]

    # Merge elements
    elements = np.concatenate((elements1, elements2))

    # Merge vertices
    x = np.concatenate((x1, x2))
    y = np.concatenate((y1, y2))
    z = np.concatenate((z1, z2))

    # Output
    md.mesh.x = x
    md.mesh.y = y
    md.mesh.z = z
    md.mesh.elements = elements
    md.mesh.numberofvertices = len(x)
    md.mesh.numberofelements = elements.shape[0]

    # Connectivities
    md.mesh.vertexconnectivity = NodeConnectivity(md.mesh.elements, md.mesh.numberofvertices)
    # print(md.mesh.vertexconnectivity)
    md.mesh.elementconnectivity = ElementConnectivity(md.mesh.elements, md.mesh.vertexconnectivity)
    # print(md.mesh.elementconnectivity)

    # Find segments
    md.mesh.segments = findsegments(md)

    # Vertex on boundary
    md.mesh.vertexonboundary = np.zeros((md.mesh.numberofvertices, ))
    md.mesh.vertexonboundary[md.mesh.segments[:, 0:2] - 1] = 1

    # Some checks
    if np.max(md.mesh.elements) > md.mesh.numberofvertices:
        raise Exception('issue in modelmerge, one of the element ids is > number of vertices!')

    return md


