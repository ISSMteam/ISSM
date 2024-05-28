import numpy as np

from ElementConnectivity import ElementConnectivity
from mesh2d import mesh2d
from NodeConnectivity import NodeConnectivity


def squaremesh(md, Lx, Ly, nx, ny):
    """SQUAREMESH - create a structured square mesh

    This script will generate a structured square mesh
    Lx and Ly are the dimension of the domain (in meters)
    nx anx ny are the number of nodes in the x and y direction
    The coordinates x and y returned are in meters.

    Usage:
        [md] = squaremesh(md, Lx, Ly, nx, ny)
    """

    #get number of elements and number of nodes
    nel = (nx - 1) * (ny - 1) * 2
    nods = nx * ny

    #initialization
    index = np.zeros((nel, 3), int)
    x = np.zeros((nx * ny, ))
    y = np.zeros((nx * ny, ))

    #create coordinates
    for n in range(0, nx):
        for m in range(0, ny):
            x[n * ny + m] = float(n)
            y[n * ny + m] = float(m)

    #create index
    for n in range(0, nx - 1):
        for m in range(0, ny - 1):
            A = n * ny + (m + 1)
            B = A + 1
            C = (n + 1) * ny + (m + 1)
            D = C + 1
            index[n * (ny - 1) * 2 + 2 * m, :] = [A, C, B]
            index[n * (ny - 1) * 2 + 2 * (m + 1) - 1, :] = [B, C, D]

    #Scale  x and y
    x = x / np.max(x) * Lx
    y = y / np.max(y) * Ly

    #create segments
    segments = np.zeros((2 * (nx - 1) + 2 * (ny - 1), 3), int)
    #left edge:
    segments[0:ny - 1, :] = np.vstack((np.arange(2, ny + 1), np.arange(1, ny), (2 * np.arange(1, ny) - 1))).T
    #right edge:
    segments[ny - 1:2 * (ny - 1), :] = np.vstack((np.arange(ny * (nx - 1) + 1, nx * ny), np.arange(ny * (nx - 1) + 2, nx * ny + 1), 2 * np.arange((ny - 1) * (nx - 2) + 1, (nx - 1) * (ny - 1) + 1))).T
    #front edge:
    segments[2 * (ny - 1):2 * (ny - 1) + (nx - 1), :] = np.vstack((np.arange(2 * ny, ny * nx + 1, ny), np.arange(ny, ny * (nx - 1) + 1, ny), np.arange(2 * (ny - 1), 2 * (nx - 1) * (ny - 1) + 1, 2 * (ny - 1)))).T
    #back edge
    segments[2 * (ny - 1) + (nx - 1):2 * (nx - 1) + 2 * (ny - 1), :] = np.vstack((np.arange(1, (nx - 2) * ny + 2, ny), np.arange(ny + 1, ny * (nx - 1) + 2, ny), np.arange(1, 2 * (nx - 2) * (ny - 1) + 2, 2 * (ny - 1)))).T

    #plug coordinates and nodes
    md.mesh = mesh2d()
    md.mesh.x = x
    md.mesh.y = y
    md.mesh.numberofvertices = nods
    md.mesh.vertexonboundary = np.zeros(nods, int)
    md.mesh.vertexonboundary[segments[:, 0:2] - 1] = 1

    #plug elements
    md.mesh.elements = index
    md.mesh.segments = segments
    md.mesh.numberofelements = nel

    #Now, build the connectivity tables for this mesh.
    md.mesh.vertexconnectivity = NodeConnectivity(md.mesh.elements, md.mesh.numberofvertices)
    md.mesh.elementconnectivity = ElementConnectivity(md.mesh.elements, md.mesh.vertexconnectivity)

    return md
