import os.path

import numpy as np

from ElementConnectivity import ElementConnectivity
from mesh2d import mesh2d
from NodeConnectivity import NodeConnectivity
from Triangle_python import Triangle_python


def triangle(md, domainname, *args):
    """triangle - create model mesh using the triangle package

    This routine creates a model mesh using Triangle and a domain outline, to 
    within a certain resolution where md is a @model object, domainname is the 
    name of an Argus domain outline file, and resolution is a characteristic 
    length for the mesh (same unit as the domain outline unit). Riftname is an optional argument (Argus domain outline) describing rifts.

    Usage:
        md = triangle(md, domainname, resolution)
        OR
        md = triangle(md, domainname, resolution, riftname)

    Examples:
        md = triangle(md, 'DomainOutline.exp', 1000)
        md = triangle(md, 'DomainOutline.exp', 1000, 'Rifts.exp')
    """

    # Figure out a characteristic area. Resolution is a node oriented concept 
    # (ex a 1000m resolution node would be made of 1000 * 1000 area squares).

    if len(args) == 1:
        resolution = args[0]
        riftname = ''
    if len(args) == 2:
        riftname = args[0]
        resolution = args[1]

    # Check that mesh was not already run, and warn user
    if md.mesh.numberofelements:
        choice = input('This model already has a mesh. Are you sure you want to go ahead? (y / n)')
        if choice not in ['y', 'n']:
            print('bad answer try you should use \'y\' or \'n\' ... exiting')
            return None
        if choice == 'n':
            print('no meshing done ... exiting')
            return None

    area = resolution ** 2

    # Check that file exists (this is a very common mistake)
    if not os.path.exists(domainname):
        raise IOError('file {} not found'.format(domainname))

    # Mesh using Triangle
    elements, x, y, segments, segmentmarkers = Triangle_python(domainname, riftname, area)

    # Check that all the created nodes belong to at least one element
    removeorphans = 1
    if removeorphans:
        uniqueelements = np.sort(np.unique(elements))
        orphans = np.nonzero((~np.isin(range(1, len(x)), uniqueelements)).astype(int))[0]
        for i in range(0, len(orphans)):
            print('WARNING: removing orphans')
            # Get rid of the orphan node i
            # Update x and y
            x = np.concatenate((x[0:(orphans[i] - i)], x[(orphans[i] - i + 1):]))
            y = np.concatenate((y[0:(orphans[i] - i)], y[(orphans[i] - i + 1):]))
            # Update elements
            pos = np.nonzero((elements > (orphans[i] - i)).flatten(order='F'))[0]
            elementstmp = elements.flatten(order='F')
            elementstmp[pos] -= 1
            elements = elementstmp.reshape(np.shape(elements), order='F')
            # Update segments
            pos1 = np.nonzero(segments[:,0] > (orphans[i] - i))[0]
            pos2 = np.nonzero(segments[:,1] > (orphans[i] - i))[0]
            segments[pos1, 0] -= 1
            segments[pos2, 1] -= 1

    # Plug into md
    md.mesh = mesh2d()
    md.mesh.x = x
    md.mesh.y = y
    md.mesh.elements = elements.astype(int)
    md.mesh.segments = segments.astype(int)
    md.mesh.segmentmarkers = segmentmarkers.astype(int)

    # Fill in rest of fields
    md.mesh.numberofelements = np.size(md.mesh.elements, axis=0)
    md.mesh.numberofvertices = np.size(md.mesh.x)
    md.mesh.vertexonboundary = np.zeros(md.mesh.numberofvertices, int)
    md.mesh.vertexonboundary[md.mesh.segments[:, 0:2] - 1] = 1

    # Now, build the connectivity tables for this mesh
    md.mesh.vertexconnectivity = NodeConnectivity(md.mesh.elements, md.mesh.numberofvertices)
    md.mesh.elementconnectivity = ElementConnectivity(md.mesh.elements, md.mesh.vertexconnectivity)

    return md
