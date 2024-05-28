import numpy as np
try:
    from scipy.sparse import csc_matrix
except ImportError:
    print("could not import scipy, no averaging capabilities enabled")

from GetAreas import GetAreas


def averaging(md, data, iterations, layer=0):
    """AVERAGING - smooths the input over the mesh

    This routine takes a list of the elements or the nodes in input and return 
    a list over the nodes.
    For each iterations it computes the average over each element (average of 
    the vertices values) and then computes the average over each node by taking 
    the average of the element around a node weighted by the elements volume.
    For 3d mesh, a last argument can be added to specify the layer to be 
    averaged on.

    Usage:
        smoothdata = averaging(md, data, iterations)
        smoothdata = averaging(md, data, iterations, layer)

    Examples:
        velsmoothed = averaging(md, md.initialization.vel, 4)
        pressure = averaging(md, md.initialization.pressure, 0)
        temperature = averaging(md, md.initialization.temperature, 1, 1)
    """

    if (len(data) != md.mesh.numberofelements) & (len(data) != md.mesh.numberofvertices):
        raise Exception('averaging error message: data not supported yet')
    if (md.mesh.dimension() == 3) & (layer != 0):
        if (layer <= 0) | (layer > md.mesh.numberoflayers):
            raise ValueError('layer should be between 1 and md.mesh.numberoflayers')
    else:
        layer = 0

    # Initialization
    if layer == 0:
        weights = np.zeros((md.mesh.numberofvertices, ))
        data = np.asarray(data).flatten()
    else:
        weights = np.zeros((md.mesh.numberofvertices2d, ))
        data = data[(layer - 1) * md.mesh.numberofvertices2d + 1:layer * md.mesh.numberofvertices2d - 1, :]

    # Load some variables (it is much faster if the variables are loaded from md once for all)
    if layer == 0:
        index = md.mesh.elements
        numberofnodes = md.mesh.numberofvertices
        numberofelements = md.mesh.numberofelements
    else:
        index = md.mesh.elements2d
        numberofnodes = md.mesh.numberofvertices2d
        numberofelements = md.mesh.numberofelements2d

    # Build some variables
    if (md.mesh.dimension() == 3) & (layer == 0):
        rep = 6
        areas = GetAreas(index, md.mesh.x, md.mesh.y, md.mesh.z)
    elif md.mesh.dimension() == 2:
        rep = 3
        areas = GetAreas(index, md.mesh.x, md.mesh.y)
    else:
        rep = 3
        areas = GetAreas(index, md.mesh.x2d, md.mesh.y2d)

    index = index - 1  # Python indexes from zero
    line=index.T.flatten()
    areas = np.vstack(areas).reshape(-1, )
    summation = 1. / rep * np.ones((rep,1) )
    linesize = rep * numberofelements

    # Update weights that hold the volume of all the element holding the node i
    weights = csc_matrix((np.tile(areas, (1, rep)).reshape(-1,), (line, np.zeros(linesize, ))), shape=(numberofnodes, 1))

    # Initialization
    if len(data) == numberofelements:
        average_node = csc_matrix((np.tile(np.multiply(areas,data), (1, rep)).reshape(-1, ), (line, np.zeros(linesize, ))), shape=(numberofnodes, 1))
        average_node = np.divide(average_node,weights)
        average_node = csc_matrix(average_node)
    else:
        average_node = csc_matrix(data.reshape(-1, 1))

    # Loop over iteration
    for i in np.arange(1, iterations + 1):
        average_el = np.asarray(average_node.todense()[index].reshape(numberofelements, rep)*summation).reshape(-1, )
        average_node = csc_matrix((np.tile(np.multiply(areas,average_el.reshape(-1)), (1, rep)).reshape(-1, ), (line, np.zeros(linesize, ))), shape=(numberofnodes, 1))
        average_node = np.divide(average_node,weights)
        average_node = csc_matrix(average_node)

    # Return output as a full matrix (C code does not like sparse matrices)
    average = np.expand_dims(np.asarray(average_node.todense()).reshape(-1, ),axis=1)

    return average
