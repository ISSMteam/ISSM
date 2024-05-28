import copy

import numpy as np

from adjacency import *
from Chaco import *
from mesh2d import *
from MeshPartition import *
from pairoptions import *
from project3d import *


def partitioner(md, *args):
    '''
    PARTITIONER - partition mesh

        List of options to partitioner:
            package: 'chaco', 'metis', or 'scotch'
            npart: number of partitions
            weighting: 'on' or 'off': default off
            section: 1 by defaults(1=bisection, 2=quadrisection, 3=octasection)
            recomputeadjacency: 'on' by default (set to 'off' to compute existing one)
            type: 'node' or 'element' partition vector (default to 'node')
            Output: partitionvector: the partition vector

        Usage:
            partitionvector, md = partitioner(md, 'package', 'chaco', 'npart', 100, 'weighting', 'on')
    '''

    #get options:
    options = pairoptions(*args)

    #set defaults
    package = options.getfieldvalue('package', 'chaco')
    npart = options.getfieldvalue('npart', 10)
    weighting = options.getfieldvalue('weighting', 'on')
    section = options.getfieldvalue('section', 1)
    recomputeadjacency = options.getfieldvalue('recomputeadjacency', 'on')
    vectortype = options.getfieldvalue('type', 'node')

    # Python only: short-circuit
    if vectortype == 'element' and not package == 'linear':
        raise RuntimeError('partitioner error message: package {} does not allow element partitions.'.format(package))

    if md.mesh.dimension() == 3:
        #partitioning essentially happens in 2D. So partition in 2D, then
        #extrude the partition vector vertically.
        md3d = copy.deepcopy(md) # save for later
        md.mesh.elements = md.mesh.elements2d
        md.mesh.x = md.mesh.x2d
        md.mesh.y = md.mesh.y2d
        md.mesh.numberofvertices = md.mesh.numberofvertices2d
        md.mesh.numberofelements = md.mesh.numberofelements2d
        md.qmu.vertex_weight = []
        md.mesh.vertexconnectivity = []
        recomputeadjacency = 'on'

    #adjacency matrix if needed:
    if recomputeadjacency == 'on':
        md = adjacency(md)
    else:
        print('skipping adjacency matrix computation as requested in the options')

    if package == 'chaco':
        if vectortype == 'element':
            raise RuntimeError('partitioner error message: package {} does not allow element partitions.'.format(package))
        else:
            # default method (from chaco.m)
            method = np.array([1, 1, 0, 0, 1, 1, 50, 0, 0.001, 7654321]).conj().transpose()
            method[0] = 3  #  global method (3 = inertial (geometric))
            method[2] = 0  #  vertex weights (0 = off, 1 = on)

            #specify bisection
            method[5] = section  #  ndims (1 = bisection, 2 = quadrisection, 3 = octasection)

            #are we using weights?
            if weighting == 'on':
                weights = np.floor(md.qmu.vertex_weight / min(md.qmu.vertex_weight))
                method[2] = 1
            else:
                weights = []

            #method = method.reshape(-1, 1)  # transpose to 1x10 instead of 10

            #  partition into nparts
            if isinstance(md.mesh, mesh2d):
                part = Chaco(md.qmu.adjacency, weights, np.array([]), md.mesh.x, md.mesh.y, np.zeros((md.mesh.numberofvertices, 1)), method, npart, np.array([]))[0].conj().transpose() + 1 #index partitions from 1 up. like metis.
            else:
                part = Chaco(md.qmu.adjacency, weights, np.array([]), md.mesh.x, md.mesh.y, md.mesh.z, method, npart, np.array([]))[0].conj().transpose() + 1 #index partitions from 1 up. like metis.
    elif package == 'scotch':
        if vectortype == 'element':
            raise RuntimeError('partitioner error message: package %s does not allow element partitions.' % package)
        else:
            #are we using weights?
            if m.strcmpi(options.getfieldvalue('weighting'), 'on'):
                weights = np.floor(md.qmu.vertex_weight / min(md.qmu.vertex_weight))
            else:
                weights = []
            maptab = Scotch(md.qmu.adjacency, [], weights, [], 'cmplt', [npart])

            part = maptab[:, 1] + 1 #index partitions from 1 up. like metis.

    elif package == 'linear':
        if vectortype == 'element':
            part = np.arange(1, 1 + md.mesh.numberofelements, 1)
            print('Linear partitioner requesting partitions on elements')
        else:
            part = np.arange(1, 1 + md.mesh.numberofvertices, 1)

    elif package == 'metis':
        if vectortype == 'element':
            raise RuntimeError('partitioner error message: package %s does not allow element partitions.' % package)
        else:
            [element_partitioning, part] = MeshPartition(md, md.qmu.numberofpartitions)

    else:
        raise RuntimeError('partitioner error message: could not find {} partitioner'.format(package))

    #extrude if we are in 3D:
    if md.mesh.dimension() == 3:
        md3d.qmu.vertex_weight = md.qmu.vertex_weight
        md3d.qmu.adjacency = md.qmu.adjacency
        md = copy.deepcopy(md3d)
        if vectortype == 'element':
            part = project3d(md, 'vector', part.conj().transpose(), 'type', 'element')
        else:
            part = project3d(md, 'vector', part.conj().transpose(), 'type', 'node')

    if part.shape[0] == 1:
        part = part.conj().transpose()

    # Output
    return part, md
