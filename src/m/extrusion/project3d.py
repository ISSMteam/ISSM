import numpy as np
from pairoptions import pairoptions


def project3d(md, *args):
    '''
    PROJECT3D - vertically project a vector from 2d mesh

        vertically project a vector from 2d mesh (split in noncoll and coll
        areas) into a 3d mesh.
        This vector can be a node vector of size (md.mesh.numberofvertices2d,
        N/A) or an element vector of size (md.mesh.numberofelements2d, N/A).

        arguments:
            'vector': 2d vector
            'type': 'element' or 'node' or 'poly'

        options:
            'layer'     a layer number where vector should keep its values. If
                        not specified, all layers adopt the value of the 2d
                        vector.
            'padding':  default to 0 (value adopted by other 3d layers not
                        being projected.
            'degree':   degree of polynomials when extrude from bottom to the top

        Examples:
            extruded_vector = project3d(md, 'vector', vector2d, 'type', 'node', 'layer', 1, 'padding', NaN)
            extruded_vector = project3d(md, 'vector', vector2d, 'type', 'element', 'padding', 0)
            extruded_vector = project3d(md, 'vector', vector2d, 'type', 'node')
    '''

    #some regular checks
    if not md:
        raise TypeError("bad usage")
    if md.mesh.domaintype().lower() != '3d':
        raise TypeError("input model is not 3d")

    #retrieve parameters from options.
    options = pairoptions(*args)
    vector2d = options.getfieldvalue('vector')  #mandatory
    vectype = options.getfieldvalue('type')  #mandatory
    layer = options.getfieldvalue('layer', 0)  #optional (do all layers otherwise)
    paddingvalue = options.getfieldvalue('padding', 0)  #0 by default
    polyexponent = options.getfieldvalue('degree', 0)  #0 by default, 0-degree polynomial

    #Handle special case where vector2d is single element (differs from representation in MATLAB)
    if isinstance(vector2d, (bool, int, float)):
        projected_vector = vector2d

    if np.size(vector2d) == 1:
        projected_vector = vector2d

    elif vectype.lower() == 'node':
        #Initialize 3d vector
        if np.ndim(vector2d) == 1:
            if vector2d.shape[0] == md.mesh.numberofvertices2d:
                projected_vector = (paddingvalue * np.ones((md.mesh.numberofvertices))).astype(vector2d.dtype)
            elif vector2d.shape[0] == md.mesh.numberofvertices2d + 1:
                projected_vector = (paddingvalue * np.ones((md.mesh.numberofvertices + 1))).astype(vector2d.dtype)
                projected_vector[-1] = vector2d[-1]
                vector2d = vector2d[:-1]
            else:
                raise TypeError("vector length not supported")
            #Fill in
            if layer == 0:
                for i in range(md.mesh.numberoflayers):
                    projected_vector[(i * md.mesh.numberofvertices2d):((i + 1) * md.mesh.numberofvertices2d)] = vector2d
            else:
                projected_vector[((layer - 1) * md.mesh.numberofvertices2d):(layer * md.mesh.numberofvertices2d)] = vector2d
        else:
            if vector2d.shape[0] == md.mesh.numberofvertices2d:
                projected_vector = (paddingvalue * np.ones((md.mesh.numberofvertices, np.size(vector2d, axis=1)))).astype(vector2d.dtype)
            elif vector2d.shape[0] == md.mesh.numberofvertices2d + 1:
                projected_vector = (paddingvalue * np.ones((md.mesh.numberofvertices + 1, np.size(vector2d, axis=1)))).astype(vector2d.dtype)
                projected_vector[-1, :] = vector2d[-1, :]
                vector2d = vector2d[:-1, :]
            else:
                raise TypeError("vector length not supported")
            #Fill in
            if layer == 0:
                for i in range(md.mesh.numberoflayers):
                    projected_vector[(i * md.mesh.numberofvertices2d):((i + 1) * md.mesh.numberofvertices2d), :] = vector2d
            else:
                projected_vector[((layer - 1) * md.mesh.numberofvertices2d):(layer * md.mesh.numberofvertices2d), :] = vector2d

    elif vectype.lower() == 'element':
        #Initialize 3d vector
        if np.ndim(vector2d) == 1:
            if vector2d.shape[0] == md.mesh.numberofelements2d:
                projected_vector = (paddingvalue * np.ones((md.mesh.numberofelements))).astype(vector2d.dtype)
            elif vector2d.shape[0] == md.mesh.numberofelements2d + 1:
                projected_vector = (paddingvalue * np.ones((md.mesh.numberofelements + 1))).astype(vector2d.dtype)
                projected_vector[-1] = vector2d[-1]
                vector2d = vector2d[:-1]
            else:
                raise TypeError("vector length not supported")
            #Fill in
            if layer == 0:
                for i in range(md.mesh.numberoflayers - 1):
                    projected_vector[(i * md.mesh.numberofelements2d):((i + 1) * md.mesh.numberofelements2d)] = vector2d
            else:
                projected_vector[((layer - 1) * md.mesh.numberofelements2d):(layer * md.mesh.numberofelements2d)] = vector2d
        else:
            if vector2d.shape[0] == md.mesh.numberofelements2d:
                projected_vector = (paddingvalue * np.ones((md.mesh.numberofelements, np.size(vector2d, axis=1)))).astype(vector2d.dtype)
            elif vector2d.shape[0] == md.mesh.numberofelements2d + 1:
                projected_vector = (paddingvalue * np.ones((md.mesh.numberofelements + 1, np.size(vector2d, axis=1)))).astype(vector2d.dtype)
                projected_vector[-1, :] = vector2d[-1, :]
                vector2d = vector2d[:-1, :]
            else:
                raise TypeError("vector length not supported")
            #Fill in
            if layer == 0:
                for i in range(md.mesh.numberoflayers - 1):
                    projected_vector[(i * md.mesh.numberofelements2d):((i + 1) * md.mesh.numberofelements2d), :] = vector2d
            else:
                projected_vector[((layer - 1) * md.mesh.numberofelements2d):(layer * md.mesh.numberofelements2d), :] = vector2d
    elif vectype.lower() == 'poly':
        #Initialize 3d vector
        if np.ndim(vector2d) == 1:
            if vector2d.shape[0] == md.mesh.numberofvertices2d:
                projected_vector = (paddingvalue * np.ones((md.mesh.numberofvertices))).astype(vector2d.dtype)
            elif vector2d.shape[0] == md.mesh.numberofvertices2d + 1:
                projected_vector = (paddingvalue * np.ones((md.mesh.numberofvertices + 1))).astype(vector2d.dtype)
                projected_vector[-1] = vector2d[-1]
                vector2d = vector2d[:-1]
            else:
                raise TypeError("vector length not supported")
            #Fill in
            if layer == 0:
                for i in range(md.mesh.numberoflayers - 1):
                    projected_vector[(i * md.mesh.numberofvertices2d):((i + 1) * md.mesh.numberofvertices2d)] = vector2d * (1.0 - (1.0 - i / (md.mesh.numberoflayers - 1.0))**polyexponent)
            else:
                projected_vector[((layer - 1) * md.mesh.numberofvertices2d):(layer * md.mesh.numberofvertices2d)] = vector2d * (1.0 - (1.0 - layer / (md.mesh.numberoflayers - 1.0))**polyexponent)
        else:
            if vector2d.shape[0] == md.mesh.numberofvertices2d:
                projected_vector = (paddingvalue * np.ones((md.mesh.numberofvertices, np.size(vector2d, axis=1)))).astype(vector2d.dtype)
            elif vector2d.shape[0] == md.mesh.numberofvertices2d + 1:
                projected_vector = (paddingvalue * np.ones((md.mesh.numberofvertices + 1, np.size(vector2d, axis=1)))).astype(vector2d.dtype)
                projected_vector[-1, :] = vector2d[-1, :]
                vector2d = vector2d[:-1, :]
            else:
                raise TypeError("vector length not supported")
            #Fill in
            if layer == 0:
                for i in range(md.mesh.numberoflayers - 1):
                    projected_vector[(i * md.mesh.numberofvertices2d):((i + 1) * md.mesh.numberofvertices2d), :] = vector2d * (1.0 - (1.0 - i / (md.mesh.numberoflayers - 1.0))**polyexponent)
            else:
                projected_vector[((layer - 1) * md.mesh.numberofvertices2d):(layer * md.mesh.numberofvertices2d), :] = vector2d * (1.0 - (1.0 - layer / (md.mesh.numberoflayers - 1.0))**polyexponent)
    else:
        raise TypeError("project3d error message: unknown projection type")

    return projected_vector
