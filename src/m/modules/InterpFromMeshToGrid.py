import numpy as np
from InterpFromMeshToGrid_python import InterpFromMeshToGrid_python

def InterpFromMeshToGrid(index, x, y, data, xgrid, ygrid, default_value):
    """InterpFromMeshToGrid - Interpolation of a data defined on a mesh onto a 
    grid

    This function is a wrapper to a multi-threaded Python module that 
    interpolates a field defined on a triangular mesh onto a regular grid.

    index, x, y:    delaunay triangulation defining the mesh
    meshdata:       vertex values of data to be interpolated

    xgrid, ygrid:   parameters that define the grid
    default_value:  value of points located out of the mesh
    """

    # Call Python module
    return np.squeeze(InterpFromMeshToGrid_python(index, x, y, data, xgrid, ygrid, default_value))
