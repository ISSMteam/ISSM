from Chaco_python import Chaco_python


def Chaco(A, vwgts, ewgts, x, y, z, options, nparts, goal):
    '''CHACO

   Usage:
      assgn = Chaco(A, vwgts, ewgts, x, y, z, options, nparts, goal)

   A:            Input adjacency matrix
   vwgts:        weights for all vertices
   ewgts:        weights for all edges
   x, y, z:        coordinates for inertial method
   options:        architecture and partitioning options
   nparts:        number of parts options
   goal:        desired set sizes
'''
    # Call mex module
    assgn = Chaco_python(A, vwgts, ewgts, x, y, z, options, nparts, goal)
    return assgn
