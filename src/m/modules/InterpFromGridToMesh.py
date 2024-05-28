import numpy as np
from InterpFromGridToMesh_python import InterpFromGridToMesh_python


def InterpFromGridToMesh(x, y, data, x_mesh, y_mesh, default_value):
    """
    INTERPFROMGRIDTOMESH - Interpolation from a grid onto a list of points

   Usage:
      data_mesh = InterpFromGridToMesh(x, y, data, x_mesh, y_mesh, default_value)

   data:        matrix holding the data to be interpolated onto the mesh
   x, y:        coordinates of matrix data (x and y must be in increasing order)
   x_mesh, y_mesh:    coordinates of the points onto which we interpolate
    default_value:    vector of mesh interpolated data

    Example:
        load('velocities.mat')
        md.inversion.vx_obs = InterpFromGridToMesh(x_n, y_m, vx, md.mesh.x, md.mesh.y, 0)
    """
    # Call mex module
    data_mesh = InterpFromGridToMesh_python(x, y, data, x_mesh, y_mesh, default_value)
    # Return
    return np.squeeze(data_mesh)
