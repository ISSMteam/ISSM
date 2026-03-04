from InterpFromMeshToMesh3d_python import InterpFromMeshToMesh3d_python


def InterpFromMeshToMesh3d(index, x, y, z, data, x_prime, y_prime, z_prime, default_value):
    """INTERPFROMMESHTOMESH3D - Interpolation from a 3d hexahedron mesh onto a list of points

    Usage:
        index:                      index of the mesh where data is defined
        x, y, z:                    coordinates of the nodes where data is defined
        data:                       matrix holding the data to be interpolated onto the mesh
        x_prime, y_prime, z_prime:  coordinates of the points onto which we interpolate
        default_value:              default value if no data is found (holes)
        data_prime:                 vector of mesh interpolated data

    Example:
        load('temperature.mat')
        md.initialization.temperature = InterpFromMeshToMesh3d(index, x, y, z, temperature, md.mesh.x, md.mesh.y, md.mesh.z, 253)
    """

    # Check usage
    #nargs = len(args)
    #if nargs != 9:
    #    print(InterpFromMeshToMesh3d.__doc__)
    #    raise Exception('Wrong usage (see above)')

    # Call Python module
    data_prime = InterpFromMeshToMesh3d_python(index, x, y, z, data, x_prime, y_prime, z_prime, default_value)

    return data_prime
