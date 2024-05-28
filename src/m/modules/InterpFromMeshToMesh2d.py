from InterpFromMeshToMesh2d_python import InterpFromMeshToMesh2d_python


def InterpFromMeshToMesh2d(*args):
    """INTERPFROMMESHTOMESH2D - Interpolation from a 2d triangular mesh onto a list of points

    Usage:
        data_interp = InterpFromMeshToMesh2d(index, x, y, data, x_interp, y_interp)
        OR
        data_interp = InterpFromMeshToMesh2d(index, x, y, data, x_interp, y_interp, OPTIONS)

    index:              index of the mesh where data is defined (e.g. md.mesh.elements)
    x, y:               coordinates of the nodes where data is defined
    data:               matrix holding the data to be interpolated onto the mesh (one column per field)
    x_interp, y_interp: coordinates of the points onto which we interpolate
    data_interp:        vector of mesh interpolated data
    Available options:
            default:    default value if point is outsite of triangulation (instead of linear interpolation)

    All fields should be arrays and not lists as this would crash the code.

    Example:
        load('temperature.mat')
        md.initialization.temperature = InterpFromMeshToMesh2d(index, x, y, temperature, md.mesh.x, md.mesh.y)
        md.initialization.temperature = InterpFromMeshToMesh2d(index, x, y, temperature, md.mesh.x, md.mesh.y, 'default', 253)
    """

    # Check usage
    nargs = len(args)
    if nargs != 6 and nargs != 8:
        print(InterpFromMeshToMesh2d.__doc__)
        raise Exception('Wrong usage (see above)')

    # Call Python module
    data_interp = InterpFromMeshToMesh2d_python(*args)
    return data_interp[0] # NOTE: Value returned from wrapper function is a tuple, the first element of which being the result we actually want
