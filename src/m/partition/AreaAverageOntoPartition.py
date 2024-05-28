import copy

import numpy as np

from adjacency import adjacency
from project2d import project2d
from qmupart2npart import qmupart2npart


def AreaAverageOntoPartition(md, vector, partition, layer=None):
    '''
    AREAAVERAGEONTOPARTITION - compute partition values for a certain vector expressed on the vertices of the mesh.
    Use area weighted average.

    Usage:
        average = AreaAverageOntoPartition(md, vector)
        average = AreaAverageOntoPartition(md, vector, layer) # If in 3D, chose which layer is partitioned
    '''

    #some checks
    if(md.mesh.dimension() == 3):
        if layer is None:
            raise RuntimeError('AreaAverageOntoPartition: layer should be provided onto which Area Averaging occurs')

        #save 3D model
        md3d = copy.deepcopy(md)

        md.mesh.elements = md.mesh.elements2d
        md.mesh.x = md.mesh.x2d
        md.mesh.y = md.mesh.y2d
        md.mesh.numberofvertices = md.mesh.numberofvertices2d
        md.mesh.numberofelements = md.mesh.numberofelements2d
        md.qmu.vertex_weight = []
        md.mesh.vertexconnectivity = []

        #run connectivity routine
        md = adjacency(md)

    #finally, project vector:
        vector = project2d(md3d, vector, layer)
        partition = project2d(md3d, partition, layer)

    #ok, first check that part is Matlab indexed
    part = partition.copy()
    part = part.flatten() + 1

    #some check:
    npart = qmupart2npart(partition)
    if npart != max(part):
        raise RuntimeError('AreaAverageOntoPartition error message: ''npart'' should be equal to max(partition)')

    #initialize output
    partvector = np.zeros((max(part)))

    #start weight average
    weightedvector = vector.flatten() * md.qmu.vertex_weight
    for i in range(max(part)):
        pos = np.where((part - 1) == i)
        partvector[i] = sum(weightedvector[pos]) / sum(md.qmu.vertex_weight[pos])

    #in 3D, restore 3D model:
    if(md.mesh.dimension() == 3):
        md = copy.deepcopy(md3d)

    return partvector
