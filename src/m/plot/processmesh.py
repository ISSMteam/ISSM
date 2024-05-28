#import matplotlib.delaunay as delaunay
import numpy as np


def processmesh(md, data, options):
    """PROCESSMESH - process mesh to be plotted

    Usage:
        x, y, z, elements, is2d = processmech(md, data, options)

    See also: PLOTMODEL, PROCESSDATA

    TODO:
    - Test application of matlplotlib.delaunay
    - Test that output of delaunay matches output of delaunay in MATLAB (see
    src/m/plot/processmesh.m)
    """

    #some checks
    if md.mesh.numberofvertices == 0:
        raise ValueError('processmesh error: mesh is empty')
    if md.mesh.numberofvertices == md.mesh.numberofelements:
        raise ValueError('processmesh error: the number of elements is the same as the number of nodes')

    #special case for mesh 2dvertical
    if md.mesh.domaintype() == '2Dvertical':
        x, y, z, elements, is2d, isplanet = processmesh(md.mesh, options)
        return x, y, z, elements, is2d, isplanet

    # # special case for mesh 3dsurface
    # if md.mesh.domaintype() == '3Dsurface':
    #     x, y, z, elements, is2d, isplanet = processmesh(md.mesh, options)
    #     if options.getfieldvalue('coord', 'xy') == 'latlon':
    #         x = md.mesh.long
    #         y = md.mesh.lat
    #         elements = delaunay(x, y)
    #         z = md.mesh.lat
    #         z[:] = 0
    #     return x, y, z, elements, is2d, isplanet

    if hasattr(md.mesh, 'elements2d'):
        elements2d = md.mesh.elements2d
        numofvertices2d = md.mesh.numberofvertices2d
        numofelements2d = md.mesh.numberofelements2d
    else:
        numofvertices2d = np.nan
        numofelements2d = np.nan

    if options.exist('amr'):
        step = options.getfieldvalue('amr')
        nonan = np.nonzero(~np.isnan(md.results.TransientSolution[step].MeshX))
        x = md.results.TransientSolution[step].MeshX[nonan]
        y = md.results.TransientSolution[step].MeshY[nonan]
        nonan = np.nonzero(~np.isnan(md.results.TransientSolution[step].MeshElements))
        elements = md.results.TransientSolution[step].MeshElements[nonan] - 1
        eldim = np.shape(md.results.TransientSolution[step].MeshElements)[1]
        elements = np.reshape(elements, ((int(len(elements) / eldim), eldim)))

    else:
        elements = md.mesh.elements - 1
        if options.getfieldvalue('coord', 'xy') != 'latlon':
            x = md.mesh.x
            if hasattr(md.mesh, 'x2d'):
                x2d = md.mesh.x2d
            y = md.mesh.y
            if hasattr(md.mesh, 'y2d'):
                y2d = md.mesh.y2d

    if hasattr(md.mesh, 'z'):
        z = md.mesh.z
    else:
        z = np.zeros(np.shape(md.mesh.x))
    z = options.getfieldvalue('z', z)
    if isinstance(z, str):
        z = getattr(md, z)

    force2D = numofelements2d in np.shape(data) or numofvertices2d in np.shape(data)
    #is it a 2D plot?
    if md.mesh.dimension() == 2 or options.getfieldvalue('layer', 0) >= 1 or force2D:
        is2d = 1
    else:
        is2d = 0

    #layer projection?
    if options.getfieldvalue('layer', 0) >= 1 or force2D:
        if options.getfieldvalue('coord', 'xy') == 'latlon':
            raise Exception('processmesh error message: cannot work with 3D meshes for now')
        #we modify the mesh temporarily to a 2d mesh from which the 3d mesh was extruded
        x = x2d
        y = y2d
        z = np.zeros(np.shape(x2d))
        elements = md.mesh.elements2d - 1

    #units
    if options.exist('unit'):
        unit = options.getfieldvalue('unit')
        x = x * unit
        y = y * unit
        z = z * unit

    #is model a member of planet class?
    #
    # TODO: Change this when planet class defined (see src/m/plot/processmesh.m)
    #
    if md.__class__.__name__ != 'model':
        isplanet = 1
    else:
        isplanet = 0

    return x, y, z, elements, is2d, isplanet
