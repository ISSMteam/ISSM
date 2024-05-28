from InterpFromMeshToGrid_python import InterpFromMeshToGrid_python


def InterpFromMeshToGrid(index, x, y, data, xgrid, ygrid, default_value):
    """
    INTERPFROMMESHTOGRID - Interpolation of a data defined on a mesh onto a grid

        This function is a multi - threaded mex file that interpolates a field defined
        on a triangular mesh onto a regular grid

        index, x, y:    delaunay triangulation defining the mesh
        meshdata:    vertex values of data to be interpolated

        xgrid, ygrid, :    parameters that define the grid
        default_value:    value of points located out of the mesh
    """
    # Call mex module
    grid = InterpFromMeshToGrid_python(index, x, y, data, xgrid, ygrid, default_value)
    # Return
    return grid
