import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

from pairoptions import pairoptions


def meshintersect3d(x, y, z, xs, ys, zs, *args):  # {{{
    """MESHINTERSECT - returns indices (into x, y, and z) of common values 
    between (x, y, z) and (xs, ys, zs) (i.e. x(index) = xs; y(index) = ys).
    """

    # Process options
    options = pairoptions(*args)

    # Retrieve tolerance
    maxtol = options.getfieldvalue('maxtol', 100000) # 100 km
    tolincrement = options.getfieldvalue('tolincrement', 10)
    force = options.getfieldvalue('force', 0)

    # Go through lats, longs and find, within tolerance, the index of the corresponding value in lat, long
    indices = np.zeros((len(xs), ))

    for i in range(len(xs)):
        tolerance = 0
        distance = ((x - xs[i]) ** 2 + (y - ys[i]) ** 2 + (z - zs[i]) ** 2) ** 0.5 # NOTE: math.sqrt cannot be applied element-wise to a list/numpy.array
        s = np.where(distance == 0)[0]

        if len(s):
            if len(s) > 1:
                # We have two vertices that are coincident! Not good
                #
                # TODO: Reconfigure the following in the process of bringing plotting online
                for j in range(len(s)):
                    plot(x[s[j]], y[s[j]], z[s[j]], c='cyan', s=40)
                print('Vertex %i of input mesh coincides with the following ouput mesh vertices ' % i)
                print(s)
                if force:
                    indices[i] = s[0]
                else:
                    raise RuntimeError('')
            else:
                indices[i] = s
        else:
            # We could not find a 0 distance, find the lowest tolerance that generates a find
            count = 1
            while not len(s):
                if count > 1000:
                    raise Exception('could not find a vertex matching vertex {} of input mesh!\nMight think anbout changing tolerance increment'.format(i))
                tolerance = tolerance + tolincrement
                s = np.where(distance < tolerance)[0]
                count = count + 1
            if tolerance > maxtol:
                raise Exception('found matching vertices {} in output mesh for input mesh vertex {}\nhowever, these vertices are farther than the max tolerance allowed!'.format(s, i))

            # Recover minimum distance
            sf = distance[s]
            pos = np.where(sf == np.min(sf))[0]
            s = s[pos]
            indices[i] = s

    if len(np.where(indices == 0)[0]) > 1: # NOTE: This check is different than the corresponding one under MATLAB as one index may indeed be '0'
        raise RuntimeError('issue with transition vector having one empty slot')

    # Convert results to type 'int' to avoid modifying structures to which 
    # results are assigned
    indices = indices.astype(int)

    return indices
# }}}
