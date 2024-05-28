from InterpFromMesh2d_python import InterpFromMesh2d_python


def InterpFromMesh2d(*args):  # {{{
    """INTERPFROMMESH2D

    Usage:
        data_prime = InterpFromMesh2d(index, x, y, data, x_prime, y_prime)
        OR
        data_prime = InterpFromMesh2d(index, x, y, data, x_prime, y_prime, default_value)
        OR
        data_prime = InterpFromMesh2d(index, x, y, data, x_prime, y_prime, defualt_value, contourname)

    index:              index of the mesh where data is defined
    x, y:               coordinates of the nodes where data is defined
    data:               vector holding the data to be intepolated onto the points
    x_prime, y_prime:   coordinates of the mesh vertices onto which we interpolate
    default_value:      a scalar or vector of size len(x_prime)
    contourname:        linear interpolation will happen on all x_interp, y_interp inside the contour, default vlaue will be adopted on the rest of the mesh
    data_prime:         vector of prime interpolated data
    """

    # Check usage
    nargs = len(args)
    if nargs != 6 and nargs != 7 and nargs != 8:
        print(InterpFromMesh2d.__doc__)
        raise Exception('Wrong usage (see above)')

    # Call Python module
    if nargs == 6:
        data_prime = InterpFromMesh2d_python(args[0], args[1], args[2], args[3], args[4], args[5])
    elif nargs == 7:
        data_prime = InterpFromMesh2d_python(args[0], args[1], args[2], args[3], args[4], args[5], args[6])
    elif nargs == 8:
        data_prime = InterpFromMesh2d_python(args[0], args[1], args[2], args[3], args[4], args[5], args[6], args[7])
    else:
        # NOTE: Should never get here because of previous check
        raise Exception('InterpFromMesh2d not supported')

    return data_prime[0] # NOTE: Value returned from wrapper function is a tuple, the first element of which being the result we actually want
# }}}
