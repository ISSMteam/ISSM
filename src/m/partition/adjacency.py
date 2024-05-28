import numpy as np

from GetAreas import *
import MatlabFuncs as m
from NodeConnectivity import *


def adjacency(md):
    """ADJACENCY - compute adjacency matrix, list of vertices and list of weights.

    function to create the adjacency matrix from the connectivity table.
    
    the required output is:
        md.adj_mat     (double [sparse nv x nv], vertex adjacency matrix)
        md.qmu.vertex_weight        (double [nv], vertex weights)
    """

    indi = np.array([md.mesh.elements[:, 0], md.mesh.elements[:, 1], md.mesh.elements[:, 2]])
    indj = np.array([md.mesh.elements[:, 1], md.mesh.elements[:, 2], md.mesh.elements[:, 0]])
    values = np.ones(np.shape(indi))

    md.qmu.adjacency = m.sparse(indi, indj, values, md.mesh.numberofvertices, md.mesh.numberofvertices)
    md.qmu.adjacency = np.logical_or(md.qmu.adjacency, md.qmu.adjacency.T).astype(float)  #change to reshape(-1, 1) if needed

    #now, build vwgt:
    areas = GetAreas(md.mesh.elements, md.mesh.x, md.mesh.y)

    #get node connectivity
    md.mesh.vertexconnectivity = NodeConnectivity(md.mesh.elements, md.mesh.numberofvertices)

    connectivity = md.mesh.vertexconnectivity[:, 0:-1]
    pos = np.where(connectivity)
    connectivity[pos] = areas[connectivity[pos] - 1] / 3.
    md.qmu.vertex_weight = np.sum(connectivity, 1)

    return md
