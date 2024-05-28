import numpy as np


def project2d(md3d, value, layer):
    '''
        returns the value of a field for a given layer of the mesh

    returns the value of a vector for a given layer from extruded mesh onto the 2d mesh
    used to do the extrusion. This function is used to compare values between different
    layers of a 3d mesh.

    Usage:
      projection_value = project2d(md3d, value, layer)

    Example:
      vel2 = project2d(md3d, md3d.initialization.vel, 2)
      returns the velocity of the second layer (1 is the base)
        '''

    if md3d.mesh.domaintype().lower() != '3d':
        raise Exception("model passed to project2d function should be 3D")

    if layer < 1 or layer > md3d.mesh.numberoflayers:
        raise ValueError("layer must be between 0 and {}".format(md3d.mesh.numberoflayers))

    # coerce to array in case float is passed
    if type(value) not in [np.ndarray, np.ma.core.MaskedArray]:
        print('coercing array')
        value = np.array(value)

    vec2d = False
    if value.ndim == 2 and value.shape[1] == 1:
        value = value.reshape(-1, )
        vec2d = True

    if value.size == 1:
        projection_value = value[(layer - 1) * md3d.mesh.numberofelements2d:layer * md3d.mesh.numberofelements2d]
    elif value.shape[0] == md3d.mesh.numberofvertices:
        #print 'indices: ', (layer - 1) * md3d.mesh.numberofvertices2d, layer * md3d.mesh.numberofvertices2d
        projection_value = value[(layer - 1) * md3d.mesh.numberofvertices2d:layer * md3d.mesh.numberofvertices2d]
    elif value.shape[0] == md3d.mesh.numberofvertices + 1:
        if np.ndim(value) == 1:
            projection_value = np.hstack((value[(layer - 1) * md3d.mesh.numberofvertices2d:layer * md3d.mesh.numberofvertices2d], value[-1]))
        else:
            projection_value = np.vstack((value[(layer - 1) * md3d.mesh.numberofvertices2d:layer * md3d.mesh.numberofvertices2d], value[-1]))
    else:
        projection_value = value[(layer - 1) * md3d.mesh.numberofelements2d:layer * md3d.mesh.numberofelements2d]

    if vec2d:
        projection_value = projection_value.reshape(-1, )

    return projection_value
