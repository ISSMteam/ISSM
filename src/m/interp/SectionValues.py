import os

import numpy as np

from expread import expread
#from InterpFromMesh2d import InterpFromMesh2d
from InterpFromMeshToMesh2d import InterpFromMeshToMesh2d
from InterpFromMeshToMesh3d import InterpFromMeshToMesh3d
from project2d import project2d


def SectionValues(md, data, infile, resolution):
    """SECTIONVALUES - compute the value of a field on a section

    This routine gets the value of a given field of the model on points
    given in the file infile (Argus type file). Resolution must be a list
    [horizontal_resolution, vertical_resolution]

    Usage:
        [elements, x, y, z, s, data] = SectionValues(md, data, filename, resolution)
        [elements, x, y, z, s, data] = SectionValues(md, data, profile_structure, resolution)
    """

    if os.path.isfile(infile):
        profile = expread(infile)[0]
        nods = profile['nods']
        x = profile['x']
        y = profile['y']
    else:
        raise IOError('file %s not found' % infile)

    #get the specified resolution
    if len(resolution) != 2:
        raise ValueError('SectionValues error message: Resolution must be a list [horizontal_resolution, vertical_resolution]')
    else:
        res_h = resolution[0]

    if md.mesh.domaintype().lower() == '3d':
        if isinstance(resolution[1], int) or isinstance(resolution[1], float):
            res_v = resolution[1]
        else:
            raise ValueError('SectionValues error: resolution must be a length - 2 list of integers or floats')

    #initialization
    X = np.array([])  #X - coordinate
    Y = np.array([])  #Y - coordinate
    S = np.array([0.])  #curvilinear coordinate

    for i in range(nods - 1):

        x_start = x[i]
        x_end = x[i + 1]
        y_start = y[i]
        y_end = y[i + 1]
        s_start = S[-1]

        length_segment = np.sqrt((x_end - x_start)**2 + (y_end - y_start)**2)
        portion = int(np.ceil(length_segment / res_h))

        x_segment = np.zeros(portion)
        y_segment = np.zeros(portion)
        s_segment = np.zeros(portion)

        for j in range(int(portion)):
            x_segment[j] = x_start + (j) * (x_end - x_start) / portion
            y_segment[j] = y_start + (j) * (y_end - y_start) / portion
            s_segment[j] = s_start + j * length_segment / portion

    #plug into X and Y
        X = np.append(X, x_segment)
        Y = np.append(Y, y_segment)
        S = np.append(S, s_segment)

    X = np.append(X, x[nods - 1])
    Y = np.append(Y, y[nods - 1])

    #Number of nodes:
    numberofnodes = X.shape[0]

    #Compute Z
    Z = np.zeros(numberofnodes)

    #New mesh and Data interpolation
    if '2d' in md.mesh.domaintype().lower():

        #Interpolation of data on specified points
        #data_interp = InterpFromMesh2d(md.mesh.elements, md.mesh.x, md.mesh.y, data, X, Y)
        data_interp = InterpFromMeshToMesh2d(md.mesh.elements, md.mesh.x, md.mesh.y, data, X, Y)
    #data_interp = griddata(md.mesh.x, md.mesh.y, data, X, Y)

    #Compute index
        index = np.array([list(range(1, numberofnodes)), list(range(2, numberofnodes + 1))]).T

    else:
        #vertically extrude mesh
        #Get base and surface for each 2d point, offset to make sure that it is inside the glacier system
        offset = 1.e-3
        base = InterpFromMeshToMesh2d(md.mesh.elements2d, md.mesh.x2d, md.mesh.y2d, project2d(md, md.geometry.base, 1), X, Y) + offset
        base = base.reshape(-1, )
        surface = InterpFromMeshToMesh2d(md.mesh.elements2d, md.mesh.x2d, md.mesh.y2d, project2d(md, md.geometry.surface, 1), X, Y) - offset
        surface = surface.reshape(-1, )

    #Some useful parameters
        layers = int(np.ceil(np.mean(md.geometry.thickness) / res_v))
        nodesperlayer = int(numberofnodes)
        nodestot = int(nodesperlayer * layers)
        elementsperlayer = int(nodesperlayer - 1)
        elementstot = int((nodesperlayer - 1) * (layers - 1))

    #initialization
        X3 = np.zeros(nodesperlayer * layers)
        Y3 = np.zeros(nodesperlayer * layers)
        Z3 = np.zeros(nodesperlayer * layers)
        S3 = np.zeros(nodesperlayer * layers)
        index3 = np.zeros((elementstot, 4))

    #Get new coordinates in 3d
        for i in range(1, layers + 1):
            X3[i - 1::layers] = X
            Y3[i - 1::layers] = Y
            Z3[i - 1::layers] = base + (i - 1) * (surface - base) / (layers - 1)
            S3[i - 1::layers] = S

            if i < layers - 1:  #Build index3 with quads
                ids = np.vstack((np.arange(i, nodestot - layers, layers), np.arange(i + 1, nodestot - layers, layers), np.arange(i + layers + 1, nodestot, layers), np.arange(i + layers, nodestot, layers))).T
                index3[(i - 1) * elementsperlayer:i * elementsperlayer, :] = ids

    #Interpolation of data on specified points
        data_interp = InterpFromMeshToMesh3d(md.mesh.elements, md.mesh.x, md.mesh.y, md.mesh.z, data, X3, Y3, Z3, np.nan)

    #build outputs
        X = X3
        Y = Y3
        Z = Z3
        S = S3

        index = index3

    return index, X, Y, Z, S, data_interp
