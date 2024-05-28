import numpy as np
from project2d import project2d


def DepthAverage(md, vector):
    '''
    computes depth average of 3d vector using the trapezoidal rule, and returns
    the value on the 2d mesh.

    Usage:
            vector_average = DepthAverage(md, vector)

    Example:
            vel_bar = DepthAverage(md, md.initialization.vel)
    '''

    #check that the model given in input is 3d
    if md.mesh.elementtype() != 'Penta':
        raise TypeError('DepthAverage error message: the model given in input must be 3d')

    # coerce to array in case float is passed
    if type(vector) not in [np.ndarray, np.ma.core.MaskedArray]:
        print('coercing array')
        vector = np.array(vector)

    vec2d = False
    if vector.ndim == 2:
        vec2d = True
        vector = vector.reshape(-1, )

    #nods data
    if vector.shape[0] == md.mesh.numberofvertices:
        vector_average = np.zeros(md.mesh.numberofvertices2d)
        for i in range(1, md.mesh.numberoflayers):
            vector_average = vector_average + (project2d(md, vector, i) + project2d(md, vector, i + 1)) / 2. * (project2d(md, md.mesh.z, i + 1) - project2d(md, md.mesh.z, i))
        vector_average = vector_average / project2d(md, md.geometry.thickness, 1)

    #element data
    elif vector.shape[0] == md.mesh.numberofelements:
        vector_average = np.zeros(md.mesh.numberofelements2d)
        for i in range(1, md.mesh.numberoflayers):
            vertices_dz = (project2d(md, md.mesh.z, i + 1) - project2d(md, md.mesh.z, i))
            elements_dz = vertices_dz.mean(1)
            vector_average = vector_average + project2d(md, vector, i) * elements_dz
    #vector_average = vector_average + project2d(md, vector, i) * (project2d(md, md.mesh.z, i + 1) - project2d(md, md.mesh.z, i))
        vertices_thickness = project2d(md, md.geometry.thickness, 1)
        elements_thickness = vertices_thickness.mean(1)
        vector_average = vector_average / elements_thickness
    #vector_average = vector_average / project2d(md, md.geometry.thickness, 1)

    else:
        raise ValueError('vector size not supported yet')

    if vec2d:
        vector_average = vector_average.reshape(-1, )

    return vector_average
