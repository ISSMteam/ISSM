import warnings
from collections import OrderedDict

import numpy as np

from pairoptions import pairoptions


def isoline(md, field, *args):
    """isoline - construct isovalue lines based on field provided

    Support options:
        'value':    isoline value, default is 0
        'output':   'struct' exp structure with individual contours (default)
                    'matrix' contours are concatenated and separated by NaNs
                    'longest' only the longest contour is returned
                    'filename.exp' saved as exp file

    Example:
        contours = isoline(md, md.results.TransientSolution[-1].MaskOceanLevelset, 'value', 0)
        contours = isoline(md, md.results.TransientSolution[-1].MaskOceanLevelset, 'output', 'matrix')
        contours = isoline(md, md.results.TransientSolution[-1].MaskIceLevelset, 'output', 'Icefront.exp')
    """

    # Process options
    options = pairoptions(*args)

    # Process data
    if md.mesh.dimension() == 3:
        x = md.mesh.x2d
        y = md.mesh.y2d
        index = md.mesh.elements2d
    else:
        x = md.mesh.x
        y = md.mesh.y
        index = md.mesh.elements

    if options.exist('amr'):
        amr = options.getfieldvalue('amr')
        x = amr.MeshX
        y = amr.MeshY
        index = amr.MeshElements

    # ISSM arrays are often (n, 1) columns: flatten so indexing behaves
    x = np.asarray(x, dtype=float).flatten()
    y = np.asarray(y, dtype=float).flatten()
    index = np.asarray(index).astype(int) - 1  # 1-based -> 0-based
    field = np.asarray(field, dtype=float).flatten()

    # Deal with z coordinate
    if hasattr(md.mesh, 'z') and not options.exist('amr'):
        z = np.asarray(md.mesh.z, dtype=float).flatten()
    else:
        z = np.zeros(len(x))

    if len(field) == 0:
        raise ValueError('field provided is empty')

    if md.mesh.dimension() == 3:
        if len(field) != md.mesh.numberofvertices2d:
            raise ValueError('field provided should be of size md.mesh.numberofvertices2d')
    else:
        if len(field) != len(x):
            raise ValueError('field provided should be of size md.mesh.numberofvertices')

    # the value we are trying to track
    level = options.getfieldvalue('value', 0.)

    # Initialization of some variables
    numberofelements = np.size(index, 0)
    elementslist = np.arange(numberofelements)

    # Get unique edges in mesh
    # 1: list of edges
    edges = np.vstack((index[:, [0, 1]], index[:, [1, 2]], index[:, [2, 0]]))
    # 2: find unique edges
    edges, J = np.unique(np.sort(edges, 1), axis=0, return_inverse=True)
    # 3: unique edge numbers
    vec = np.ravel(J)
    # 4: unique edges numbers in each triangle (2 triangles sharing the same
    #    edge will have the same edge number)
    edges_tria = np.column_stack((vec[elementslist], vec[elementslist + numberofelements], vec[elementslist + 2 * numberofelements]))

    # Segments [nodes1 nodes2]
    Seg1 = index[:, [0, 1]]
    Seg2 = index[:, [1, 2]]
    Seg3 = index[:, [2, 0]]

    # Segment numbers [1;4;6;...]
    Seg1_num = edges_tria[:, 0]
    Seg2_num = edges_tria[:, 1]
    Seg3_num = edges_tria[:, 2]

    # Value of data on each tips of the segments
    Data1 = field[Seg1]
    Data2 = field[Seg2]
    Data3 = field[Seg3]

    # Get the ranges for each segment
    Range1 = np.sort(Data1, 1)
    Range2 = np.sort(Data2, 1)
    Range3 = np.sort(Data3, 1)

    # Find the segments that contain this value
    pos1 = (Range1[:, 0] < level) & (Range1[:, 1] >= level)
    pos2 = (Range2[:, 0] < level) & (Range2[:, 1] >= level)
    pos3 = (Range3[:, 0] < level) & (Range3[:, 1] >= level)

    # Get elements
    poselem12 = (pos1 & pos2)
    poselem13 = (pos1 & pos3)
    poselem23 = (pos2 & pos3)
    poselem = np.flatnonzero(poselem12 | poselem13 | poselem23)
    numelems = len(poselem)

    # If no element has been flagged, skip to the next level
    if numelems == 0:
        warnings.warn('isoline warning message: no elements found with corresponding value')
        return []

    # Go through the elements and build the coordinates for each segment (1 by element)
    x1 = np.zeros(numelems)
    x2 = np.zeros(numelems)
    y1 = np.zeros(numelems)
    y2 = np.zeros(numelems)
    z1 = np.zeros(numelems)
    z2 = np.zeros(numelems)

    edge_l = np.zeros((numelems, 2), dtype=int)

    for j in range(numelems):

        with np.errstate(divide='ignore', invalid='ignore'):
            weight1 = (level - Data1[poselem[j], 0]) / (Data1[poselem[j], 1] - Data1[poselem[j], 0])
            weight2 = (level - Data2[poselem[j], 0]) / (Data2[poselem[j], 1] - Data2[poselem[j], 0])
            weight3 = (level - Data3[poselem[j], 0]) / (Data3[poselem[j], 1] - Data3[poselem[j], 0])

        if poselem12[poselem[j]]:

            x1[j] = x[Seg1[poselem[j], 0]] + weight1 * (x[Seg1[poselem[j], 1]] - x[Seg1[poselem[j], 0]])
            x2[j] = x[Seg2[poselem[j], 0]] + weight2 * (x[Seg2[poselem[j], 1]] - x[Seg2[poselem[j], 0]])
            y1[j] = y[Seg1[poselem[j], 0]] + weight1 * (y[Seg1[poselem[j], 1]] - y[Seg1[poselem[j], 0]])
            y2[j] = y[Seg2[poselem[j], 0]] + weight2 * (y[Seg2[poselem[j], 1]] - y[Seg2[poselem[j], 0]])
            z1[j] = z[Seg1[poselem[j], 0]] + weight1 * (z[Seg1[poselem[j], 1]] - z[Seg1[poselem[j], 0]])
            z2[j] = z[Seg2[poselem[j], 0]] + weight2 * (z[Seg2[poselem[j], 1]] - z[Seg2[poselem[j], 0]])

            edge_l[j, 0] = Seg1_num[poselem[j]]
            edge_l[j, 1] = Seg2_num[poselem[j]]

        elif poselem13[poselem[j]]:

            x1[j] = x[Seg1[poselem[j], 0]] + weight1 * (x[Seg1[poselem[j], 1]] - x[Seg1[poselem[j], 0]])
            x2[j] = x[Seg3[poselem[j], 0]] + weight3 * (x[Seg3[poselem[j], 1]] - x[Seg3[poselem[j], 0]])
            y1[j] = y[Seg1[poselem[j], 0]] + weight1 * (y[Seg1[poselem[j], 1]] - y[Seg1[poselem[j], 0]])
            y2[j] = y[Seg3[poselem[j], 0]] + weight3 * (y[Seg3[poselem[j], 1]] - y[Seg3[poselem[j], 0]])
            z1[j] = z[Seg1[poselem[j], 0]] + weight1 * (z[Seg1[poselem[j], 1]] - z[Seg1[poselem[j], 0]])
            z2[j] = z[Seg3[poselem[j], 0]] + weight3 * (z[Seg3[poselem[j], 1]] - z[Seg3[poselem[j], 0]])

            edge_l[j, 0] = Seg1_num[poselem[j]]
            edge_l[j, 1] = Seg3_num[poselem[j]]

        elif poselem23[poselem[j]]:

            x1[j] = x[Seg2[poselem[j], 0]] + weight2 * (x[Seg2[poselem[j], 1]] - x[Seg2[poselem[j], 0]])
            x2[j] = x[Seg3[poselem[j], 0]] + weight3 * (x[Seg3[poselem[j], 1]] - x[Seg3[poselem[j], 0]])
            y1[j] = y[Seg2[poselem[j], 0]] + weight2 * (y[Seg2[poselem[j], 1]] - y[Seg2[poselem[j], 0]])
            y2[j] = y[Seg3[poselem[j], 0]] + weight3 * (y[Seg3[poselem[j], 1]] - y[Seg3[poselem[j], 0]])
            z1[j] = z[Seg2[poselem[j], 0]] + weight2 * (z[Seg2[poselem[j], 1]] - z[Seg2[poselem[j], 0]])
            z2[j] = z[Seg3[poselem[j], 0]] + weight3 * (z[Seg3[poselem[j], 1]] - z[Seg3[poselem[j], 0]])

            edge_l[j, 0] = Seg2_num[poselem[j]]
            edge_l[j, 1] = Seg3_num[poselem[j]]

        else:
            pass  # should never get here

    # Now that we have the segments, we must try to connect them...

    # Loop over the subcontours
    contours = []

    while len(edge_l) > 0:

        # Take the right edge of the second segment and connect it to the next segments if any
        e1 = edge_l[0, 0]
        e2 = edge_l[0, 1]
        xc = np.array([x1[0], x2[0]])
        yc = np.array([y1[0], y2[0]])
        zc = np.array([z1[0], z2[0]])

        # Erase the lines corresponding to this edge
        edge_l = np.delete(edge_l, 0, axis=0)
        x1 = np.delete(x1, 0)
        x2 = np.delete(x2, 0)
        y1 = np.delete(y1, 0)
        y2 = np.delete(y2, 0)
        z1 = np.delete(z1, 0)
        z2 = np.delete(z2, 0)

        ro1, co1 = np.where(edge_l == e1)

        while len(ro1) > 0:

            if np.all(co1 == 0):
                xc = np.concatenate((x2[ro1], xc))
                yc = np.concatenate((y2[ro1], yc))
                zc = np.concatenate((z2[ro1], zc))

                # Next edge:
                e1 = edge_l[ro1, 1]

            else:
                xc = np.concatenate((x1[ro1], xc))
                yc = np.concatenate((y1[ro1], yc))
                zc = np.concatenate((z1[ro1], zc))

                # Next edge:
                e1 = edge_l[ro1, 0]

            # Erase the lines of this
            edge_l = np.delete(edge_l, ro1, axis=0)
            x1 = np.delete(x1, ro1)
            x2 = np.delete(x2, ro1)
            y1 = np.delete(y1, ro1)
            y2 = np.delete(y2, ro1)
            z1 = np.delete(z1, ro1)
            z2 = np.delete(z2, ro1)

            # Next connection
            ro1, co1 = np.where(edge_l == e1)

        # Same thing the other way (to the right)
        ro2, co2 = np.where(edge_l == e2)

        while len(ro2) > 0:

            if np.all(co2 == 0):
                xc = np.concatenate((xc, x2[ro2]))
                yc = np.concatenate((yc, y2[ro2]))
                zc = np.concatenate((zc, z2[ro2]))

                # Next edge:
                e2 = edge_l[ro2, 1]

            else:
                xc = np.concatenate((xc, x1[ro2]))
                yc = np.concatenate((yc, y1[ro2]))
                zc = np.concatenate((zc, z1[ro2]))

                # Next edge:
                e2 = edge_l[ro2, 0]

            # Erase the lines of this
            edge_l = np.delete(edge_l, ro2, axis=0)
            x1 = np.delete(x1, ro2)
            x2 = np.delete(x2, ro2)
            y1 = np.delete(y1, ro2)
            y2 = np.delete(y2, ro2)
            z1 = np.delete(z1, ro2)
            z2 = np.delete(z2, ro2)

            # Next connection
            ro2, co2 = np.where(edge_l == e2)

        # Save xc, yc contour
        newcontour = OrderedDict()
        newcontour['nods'] = len(xc)
        newcontour['density'] = 1
        newcontour['closed'] = 0
        newcontour['x'] = xc
        newcontour['y'] = yc
        newcontour['z'] = zc
        newcontour['name'] = ''
        contours.append(newcontour)

    # Process output
    outputformat = options.getfieldvalue('output', 'struct')
    if outputformat == 'matrix':
        xm = np.concatenate([np.append(c['x'], np.nan) for c in contours])
        ym = np.concatenate([np.append(c['y'], np.nan) for c in contours])
        contours = np.column_stack((xm, ym))
    elif outputformat.endswith('.exp'):
        print('Saving output as ' + outputformat)
        from expwrite import expwrite
        expwrite(contours, outputformat)
    elif outputformat == 'struct':
        pass  # nothing to do, this is the default
    elif outputformat == 'longest':
        contours = contours[int(np.argmax([c['nods'] for c in contours]))]
    else:
        print('output format not supported, returning struct')

    return contours
