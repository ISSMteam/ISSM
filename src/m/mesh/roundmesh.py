import numpy as np
import os
from collections import OrderedDict
from expwrite import expwrite
from MatlabFuncs import tempname
from triangle import triangle


def roundmesh(md, radius, resolution, *args):
    """roundmesh - create an unstructured round mesh

    This script will generate an unstructured round mesh
    - radius     : specifies the radius of the circle in meters
    - resolution : specifies the resolution in meters

    Usage:
        md = roundmesh(md, radius, resolution)
        md = roundmesh(md, radius, resolution, 'domain.exp')
    """
    # First we have to create the domain outline
    if len(args):
        expname = args[0]
    else:
        expname = tempname() + '.exp'

    # Get number of points on the circle
    pointsonedge = int(np.floor((2. * np.pi * radius) / resolution) + 1)  # +1 to close the outline

    # Calculate the Cartesian coordinates of the points
    theta = np.linspace(0., 2. * np.pi, pointsonedge)
    x_list = roundsigfig(radius * np.cos(theta), 12)
    y_list = roundsigfig(radius * np.sin(theta), 12)
    A = OrderedDict()
    A['x'] = x_list
    A['y'] = y_list
    A['density'] = 1.
    expwrite(A, expname)

    # Call mesher
    md = triangle(md, expname, resolution)
    # md = bamg(md, 'domain', 'RoundDomainOutline.exp', 'hmin', resolution)

    # Move the closest node to the center
    pos = np.argmin(np.add(np.power(md.mesh.x, 2), np.power(md.mesh.y, 2)))
    md.mesh.x[pos] = 0.
    md.mesh.y[pos] = 0.

    # Delete domain
    if not len(args):
        os.remove(expname)

    return md


def roundsigfig(x, n):
    nonzeros = np.where(x != 0)
    digits = np.ceil(np.log10(np.abs(x[nonzeros])))
    x[nonzeros] = x[nonzeros] / 10.**digits
    x[nonzeros] = np.round(x[nonzeros], decimals=n)
    x[nonzeros] = x[nonzeros] * 10.**digits
    return x
