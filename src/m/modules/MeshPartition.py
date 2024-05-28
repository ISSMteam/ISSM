from MeshPartition_python import MeshPartition_python
from mesh3dprisms import *
from mesh2d import *
from mesh2dvertical import *


def MeshPartition(md, numpartitions):
    '''MESHPARTITION - Partition mesh according to the number of areas, using Metis library.

       Usage:
            [element_partitioning, node_partitioning] = MeshPartition(md.mesh, numpartitions)

       element_partitioning: Vector of partitioning area numbers, for every element.
       node_partitioning: Vector of partitioning area numbers, for every node.
'''
    if md is None or numpartitions is None:
        print((MeshPartition.__doc__))
        raise RuntimeError('Wrong usage (see above)')

    #Get mesh info from md.mesh
    numberofvertices = md.mesh.numberofvertices
    elements = md.mesh.elements
    numberofvertices2d = 0
    numberoflayers = 1
    elements2d = []
    if isinstance(md.mesh, mesh3dprisms):
        elementtype = 'Penta'
        numberofvertices2d = md.mesh.numberofvertices2d
        numberoflayers = md.mesh.numberoflayers
        elements2d = md.mesh.elements2d
    elif isinstance(md.mesh, mesh2d):
        elementtype = 'Tria'
    elif isinstance(md.mesh, mesh2dvertical):
        elementtype = 'Tria'

    #Call module
    [element_partitioning, node_partitioning] = MeshPartition_python(numberofvertices, elements, numberofvertices2d, elements2d, numberoflayers, elementtype, numpartitions)

    return [element_partitioning, node_partitioning]
